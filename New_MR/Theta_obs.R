library(MASS)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gt)


expit <- function(x) {
  1 / (1 + exp(-x))
}

logit <- function(p) {
  log(p / (1 - p))
}

solve_psd <- function(A, tol = 1e-10) {
  A <- (A + t(A)) / 2
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  
  if (min(ev) < tol) {
    A <- A + diag(tol - min(ev) + tol, nrow(A))
  }
  
  solve(A)
}

trim_prob <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}

as_binary <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  as.numeric(x)
}


build_mr_data <- function(dat, z, x, y) {
  vars <- unique(c(z, x, y))
  dat <- dat[, vars, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  if (nrow(dat) == 0) {
    stop("No complete cases remain after filtering.")
  }
  
  for (v in c(x, y, z)) {
    dat[[v]] <- as_binary(dat[[v]])
  }
  
  if (!all(dat[[x]] %in% c(0, 1))) {
    stop(sprintf("%s must be coded 0/1.", x))
  }
  
  if (!all(dat[[y]] %in% c(0, 1))) {
    stop(sprintf("%s must be coded 0/1.", y))
  }
  
  for (v in z) {
    if (!all(dat[[v]] %in% c(0, 1))) {
      stop(sprintf("IV %s must be coded 0/1.", v))
    }
    if (stats::var(dat[[v]]) == 0) {
      stop(sprintf("IV %s has no variation.", v))
    }
  }
  
  dat
}


fit_first_stage <- function(dat, z, x) {
  fml <- as.formula(
    paste(x, "~", paste(z, collapse = " + "))
  )
  
  glm(
    fml,
    data = dat,
    family = binomial()
  )
}


get_F_diag <- function(dat, z, x, eps = 1e-10) {
  out_fail <- tibble::tibble(
    R2_mcf = NA_real_,
    R2_cs = NA_real_,
    R2_nk = NA_real_,
    F_mcf = NA_real_,
    F_cs = NA_real_,
    F_nk = NA_real_,
    F_anova = NA_real_
  )
  
  n <- nrow(dat)
  n_iv <- length(z)
  
  if (n <= n_iv + 1 || n_iv < 1) return(out_fail)
  
  fml_null <- as.formula(paste(x, "~ 1"))
  fml_full <- as.formula(paste(x, "~", paste(z, collapse = " + ")))
  
  fit_null <- tryCatch(
    glm(fml_null, data = dat, family = binomial()),
    error = function(e) NULL
  )
  
  fit_full <- tryCatch(
    glm(fml_full, data = dat, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit_null) || is.null(fit_full)) return(out_fail)
  
  ll_null <- tryCatch(as.numeric(logLik(fit_null)), error = function(e) NA_real_)
  ll_full <- tryCatch(as.numeric(logLik(fit_full)), error = function(e) NA_real_)
  
  if (!is.finite(ll_null) || !is.finite(ll_full) || abs(ll_null) < eps) {
    return(out_fail)
  }
  
  R2_mcf <- 1 - (ll_full / ll_null)
  R2_cs <- 1 - exp((2 / n) * (ll_null - ll_full))
  R2_nk_denom <- 1 - exp((2 / n) * ll_null)
  
  R2_nk <- if (is.finite(R2_nk_denom) && R2_nk_denom > eps) {
    R2_cs / R2_nk_denom
  } else {
    NA_real_
  }
  
  calc_F <- function(R2) {
    if (!is.finite(R2) || R2 <= 0 || R2 >= 1) return(NA_real_)
    (R2 / (1 - R2)) * ((n - n_iv - 1) / n_iv)
  }
  
  F_mcf <- calc_F(R2_mcf)
  F_cs  <- calc_F(R2_cs)
  F_nk  <- calc_F(R2_nk)
  
  fit_lm <- tryCatch(
    lm(fml_full, data = dat),
    error = function(e) NULL
  )
  
  F_anova <- if (!is.null(fit_lm)) {
    fs <- tryCatch(summary(fit_lm)$fstatistic, error = function(e) NULL)
    if (is.null(fs)) NA_real_ else unname(fs[1])
  } else {
    NA_real_
  }
  
  tibble::tibble(
    R2_mcf = as.numeric(R2_mcf),
    R2_cs = as.numeric(R2_cs),
    R2_nk = as.numeric(R2_nk),
    F_mcf = as.numeric(F_mcf),
    F_cs = as.numeric(F_cs),
    F_nk = as.numeric(F_nk),
    F_anova = as.numeric(F_anova)
  )
}

## ------------------------------------------------------------
## 2SPS
## ------------------------------------------------------------

est_2sps <- function(dat, z, x, y, alpha = 0.05, eps = 1e-8) {
  
  F_fail <- c(
    R2_mcf = NA_real_,
    R2_cs = NA_real_,
    R2_nk = NA_real_,
    F_mcf = NA_real_,
    F_cs = NA_real_,
    F_nk = NA_real_,
    F_anova = NA_real_
  )
  
  out_fail <- c(
    method = "2SPS",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0,
    F_fail
  )
  
  fit_1 <- tryCatch(
    fit_first_stage(dat, z, x),
    error = function(e) NULL
  )
  
  if (is.null(fit_1)) return(out_fail)
  
  ## ----------------------------------------------------------
  ## First-stage R2 and F diagnostics calculated during 2SPS
  ## ----------------------------------------------------------
  
  F_diag_2sps <- tryCatch(
    get_F_diag(dat, z, x),
    error = function(e) tibble::as_tibble_row(as.list(F_fail))
  )
  
  F_diag_2sps <- as.numeric(F_diag_2sps[1, names(F_fail), drop = TRUE])
  names(F_diag_2sps) <- names(F_fail)
  
  dat$x_hat <- tryCatch(
    predict(fit_1, type = "response"),
    error = function(e) rep(NA_real_, nrow(dat))
  )
  
  if (any(!is.finite(dat$x_hat))) {
    out_fail[names(F_fail)] <- F_diag_2sps
    return(out_fail)
  }
  
  fit_2 <- tryCatch(
    glm(
      as.formula(paste(y, "~ x_hat")),
      data = dat,
      family = binomial()
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_2)) {
    out_fail[names(F_fail)] <- F_diag_2sps
    return(out_fail)
  }
  
  ct <- summary(fit_2)$coefficients
  
  if (!("x_hat" %in% rownames(ct))) {
    out_fail[names(F_fail)] <- F_diag_2sps
    return(out_fail)
  }
  
  est <- unname(ct["x_hat", "Estimate"])
  se  <- unname(ct["x_hat", "Std. Error"])
  z   <- unname(ct["x_hat", "z value"])
  p   <- unname(ct["x_hat", "Pr(>|z|)"])
  
  if (!is.finite(est) || !is.finite(se) || se <= eps) {
    out_fail[names(F_fail)] <- F_diag_2sps
    return(out_fail)
  }
  
  zcrit <- qnorm(1 - alpha / 2)
  lcl <- est - zcrit * se
  ucl <- est + zcrit * se
  
  c(
    method = "2SPS",
    est = est,
    se = se,
    lcl = lcl,
    ucl = ucl,
    z = z,
    p = p,
    conv = 1,
    F_diag_2sps
  )
}


## ------------------------------------------------------------
## 2SRI
## ------------------------------------------------------------

est_2sri <- function(dat, z, x, y, alpha = 0.05, eps = 1e-8) {
  out_fail <- c(
    method = "2SRI",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0
  )
  
  fit_1 <- tryCatch(
    fit_first_stage(dat, z, x),
    error = function(e) NULL
  )
  if (is.null(fit_1)) return(out_fail)
  
  dat$x_hat <- tryCatch(
    predict(fit_1, type = "response"),
    error = function(e) rep(NA_real_, nrow(dat))
  )
  if (any(!is.finite(dat$x_hat))) return(out_fail)
  
  dat$r_1 <- dat[[x]] - dat$x_hat
  
  fit_2 <- tryCatch(
    glm(
      as.formula(paste(y, "~", x, "+ r_1")),
      data = dat,
      family = binomial()
    ),
    error = function(e) NULL
  )
  if (is.null(fit_2)) return(out_fail)
  
  ct <- summary(fit_2)$coefficients
  if (!(x %in% rownames(ct))) return(out_fail)
  
  est <- unname(ct[x, "Estimate"])
  se  <- unname(ct[x, "Std. Error"])
  z   <- unname(ct[x, "z value"])
  p   <- unname(ct[x, "Pr(>|z|)"])
  
  if (!is.finite(est) || !is.finite(se) || se <= eps) {
    return(out_fail)
  }
  
  zcrit <- qnorm(1 - alpha / 2)
  lcl <- est - zcrit * se
  ucl <- est + zcrit * se
  
  c(
    method = "2SRI",
    est = est,
    se = se,
    lcl = lcl,
    ucl = ucl,
    z = z,
    p = p,
    conv = 1
  )
}

## ------------------------------------------------------------
## GMM
## ------------------------------------------------------------

est_gmm <- function(dat, z, x, y,
                    start = c(0, 0),
                    alpha = 0.05,
                    maxit = 1000,
                    eps = 1e-8) {
  y_vec <- dat[[y]]
  X_mat <- cbind(1, dat[[x]])
  Z_mat <- cbind(1, as.matrix(dat[, z, drop = FALSE]))
  n <- nrow(dat)
  
  if (length(y_vec) != n || nrow(X_mat) != n || nrow(Z_mat) != n) {
    stop("Input dimensions do not match.")
  }
  
  p <- ncol(X_mat)
  
  if (length(start) != p) {
    stop("Length of 'start' must equal the number of regression parameters.")
  }
  
  out_fail <- c(
    method = "GMM",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0,
    obj = NA_real_
  )
  
  W <- tryCatch(
    solve_psd(crossprod(Z_mat) / n),
    error = function(e) NULL
  )
  if (is.null(W) || any(!is.finite(W))) return(out_fail)
  
  g_bar <- function(b) {
    mu <- expit(drop(X_mat %*% b))
    r <- y_vec - mu
    drop(crossprod(Z_mat, r) / n)
  }
  
  obj_gmm <- function(b) {
    g <- g_bar(b)
    if (any(!is.finite(g))) return(1e12)
    
    val <- as.numeric(t(g) %*% W %*% g)
    if (!is.finite(val)) 1e12 else val
  }
  
  opt <- tryCatch(
    optim(
      par = start,
      fn = obj_gmm,
      method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-10)
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || is.null(opt$par) || any(!is.finite(opt$par))) {
    return(out_fail)
  }
  
  b_hat <- opt$par
  
  if (opt$convergence != 0) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  mu_hat <- expit(drop(X_mat %*% b_hat))
  r_hat <- y_vec - mu_hat
  
  if (any(!is.finite(mu_hat)) || any(mu_hat <= 0) || any(mu_hat >= 1)) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  G_hat <- Z_mat * as.numeric(r_hat)
  S_hat <- crossprod(G_hat) / n
  
  w_hat <- as.numeric(mu_hat * (1 - mu_hat))
  D_hat <- -crossprod(Z_mat, X_mat * w_hat) / n
  
  B_mat <- t(D_hat) %*% W %*% D_hat
  B_inv <- tryCatch(solve_psd(B_mat), error = function(e) NULL)
  
  if (is.null(B_inv) || any(!is.finite(B_inv))) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  M_mat <- t(D_hat) %*% W %*% S_hat %*% W %*% D_hat
  
  if (any(!is.finite(M_mat))) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  vcov_b <- (B_inv %*% M_mat %*% B_inv) / n
  
  if (any(!is.finite(vcov_b))) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  var_b <- diag(vcov_b)
  var_b[var_b < 0] <- NA_real_
  se_b <- sqrt(var_b)
  
  est <- b_hat[2]
  se <- se_b[2]
  
  if (!is.finite(est) || !is.finite(se) || is.na(se) || se <= eps) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  z_stat <- est / se
  p_val <- 2 * (1 - pnorm(abs(z_stat)))
  zcrit <- qnorm(1 - alpha / 2)
  lcl <- est - zcrit * se
  ucl <- est + zcrit * se
  
  c(
    method = "GMM",
    est = est,
    se = se,
    lcl = lcl,
    ucl = ucl,
    z = z_stat,
    p = p_val,
    conv = 1,
    obj = opt$value
  )
}

## ------------------------------------------------------------
## IV-MVB
## ------------------------------------------------------------

est_ivmvb <- function(dat, z, x, y,
                      start = NULL,
                      alpha = 0.05,
                      maxit = 1000,
                      eps = 1e-8) {
  
  x_obs <- dat[[x]]
  y_obs <- dat[[y]]
  Z_mat <- as.matrix(dat[, z, drop = FALSE])
  k <- ncol(Z_mat)
  
  p_x_obs  <- mean(x_obs, na.rm = TRUE)
  p_y_obs  <- mean(y_obs, na.rm = TRUE)
  p11_obs  <- mean(x_obs * y_obs, na.rm = TRUE)
  
  denom_theta_obs <- p_x_obs * p_y_obs * (1 - p_x_obs) * (1 - p_y_obs)
  denom_theta_obs <- max(denom_theta_obs, eps)
  
  theta_obs <- (p11_obs - p_x_obs * p_y_obs) / denom_theta_obs
  
  if (!is.finite(theta_obs)) {
    theta_obs <- NA_real_
  }
  
  if (is.null(start)) {
    start <- rep(0, k + 3)
  }
  
  out_fail <- c(
    method = "IV-MVB",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0,
    obj = NA_real_,
    theta_obs = theta_obs
  )
  
  if (!is.finite(theta_obs)) {
    return(out_fail)
  }
  
  nll <- function(par) {
    a0 <- par[1]
    a  <- par[2:(k + 1)]
    b0 <- par[k + 2]
    b1 <- par[k + 3]
    
    pi_x <- trim_prob(expit(as.numeric(a0 + Z_mat %*% a)), eps)
    pi_y <- trim_prob(expit(b0 + b1 * pi_x), eps)
    
    p11 <- pi_x * pi_y * (1 + theta_obs * (1 - pi_x) * (1 - pi_y))
    p10 <- pi_x * (1 - pi_y) * (1 - theta_obs * (1 - pi_x) * pi_y)
    p01 <- (1 - pi_x) * pi_y * (1 - theta_obs * pi_x * (1 - pi_y))
    p00 <- (1 - pi_x) * (1 - pi_y) * (1 + theta_obs * pi_x * pi_y)
    
    if (any(!is.finite(p11)) || any(!is.finite(p10)) ||
        any(!is.finite(p01)) || any(!is.finite(p00))) {
      return(1e12)
    }
    
    if (any(p11 <= 0) || any(p10 <= 0) ||
        any(p01 <= 0) || any(p00 <= 0)) {
      return(1e12)
    }
    
    p_sum <- p11 + p10 + p01 + p00
    if (any(abs(p_sum - 1) > 1e-6)) {
      return(1e12)
    }
    
    ll <- x_obs * y_obs * log(p11) +
      x_obs * (1 - y_obs) * log(p10) +
      (1 - x_obs) * y_obs * log(p01) +
      (1 - x_obs) * (1 - y_obs) * log(p00)
    
    if (any(!is.finite(ll))) {
      return(1e12)
    }
    
    -sum(ll)
  }
  
  opt <- tryCatch(
    optim(
      par = start,
      fn = nll,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = maxit, reltol = 1e-10)
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || is.null(opt$par)) {
    return(out_fail)
  }
  
  conv <- as.integer(opt$convergence == 0)
  b_hat <- opt$par
  
  if (conv == 0 || any(!is.finite(b_hat))) {
    out <- out_fail
    out["obj"] <- if (!is.null(opt$value)) opt$value else NA_real_
    return(out)
  }
  
  se_b <- rep(NA_real_, length(b_hat))
  if (!is.null(opt$hessian)) {
    V <- tryCatch(solve_psd(opt$hessian), error = function(e) NULL)
    if (!is.null(V)) {
      var_b <- diag(V)
      var_b[var_b < 0] <- NA_real_
      se_b <- sqrt(var_b)
    }
  }
  
  est <- b_hat[k + 3]
  se  <- se_b[k + 3]
  
  if (!is.finite(est) || !is.finite(se) || is.na(se) || se <= eps) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  z_stat <- est / se
  p_val  <- 2 * (1 - pnorm(abs(z_stat)))
  zcrit  <- qnorm(1 - alpha / 2)
  lcl    <- est - zcrit * se
  ucl    <- est + zcrit * se
  
  c(
    method = "IV-MVB",
    est = est,
    se = se,
    lcl = lcl,
    ucl = ucl,
    z = z_stat,
    p = p_val,
    conv = 1,
    obj = opt$value,
    theta_obs = theta_obs
  )
}

MRbinary <- function(dat,
                     z,
                     x,
                     y,
                     alpha = 0.05,
                     maxit = 1000) {
  
  dat <- build_mr_data(dat, z = z, x = x, y = y)
  
  safe_est <- function(expr, method) {
    tryCatch(
      expr,
      error = function(e) c(
        method = method,
        est = NA_real_,
        se = NA_real_,
        lcl = NA_real_,
        ucl = NA_real_,
        z = NA_real_,
        p = NA_real_,
        conv = 0
      )
    )
  }
  
  res <- list()
  
  res$stage1 <- tryCatch(
    fit_first_stage(dat, z, x),
    error = function(e) NULL
  )
  
  res$`2SPS` <- safe_est(
    est_2sps(dat, z, x, y, alpha = alpha),
    "2SPS"
  )
  
  res$`2SRI` <- safe_est(
    est_2sri(dat, z, x, y, alpha = alpha),
    "2SRI"
  )
  
  res$GMM <- safe_est(
    est_gmm(dat, z, x, y, alpha = alpha, maxit = maxit),
    "GMM"
  )
  
  res$`IV-MVB` <- safe_est(
    est_ivmvb(dat, z, x, y, alpha = alpha, maxit = maxit),
    "IV-MVB"
  )
  
  tab <- do.call(
    rbind,
    lapply(res[c("2SPS", "2SRI", "GMM", "IV-MVB")], function(o) {
      
      get_num <- function(name) {
        if (name %in% names(o)) {
          as.numeric(unname(o[name]))
        } else {
          NA_real_
        }
      }
      
      data.frame(
        method = as.character(o["method"]),
        est    = get_num("est"),
        se     = get_num("se"),
        lcl    = get_num("lcl"),
        ucl    = get_num("ucl"),
        z      = get_num("z"),
        p      = get_num("p"),
        conv   = get_num("conv"),
        
        R2_mcf_2sps = get_num("R2_mcf"),
        R2_cs_2sps  = get_num("R2_cs"),
        R2_nk_2sps  = get_num("R2_nk"),
        F_mcf_2sps  = get_num("F_mcf"),
        F_cs_2sps   = get_num("F_cs"),
        F_nk_2sps   = get_num("F_nk"),
        F_anova_2sps = get_num("F_anova"),
        
        stringsAsFactors = FALSE
      )
    })
  )
  
  zcrit <- qnorm(1 - alpha / 2)
  ok <- tab$conv == 1 & is.finite(tab$est) & is.finite(tab$se) & tab$se > 0
  
  tab$lcl[ok] <- tab$est[ok] - zcrit * tab$se[ok]
  tab$ucl[ok] <- tab$est[ok] + zcrit * tab$se[ok]
  tab$z[ok]   <- tab$est[ok] / tab$se[ok]
  tab$p[ok]   <- 2 * (1 - pnorm(abs(tab$z[ok])))
  
  tab$lcl[!ok] <- NA_real_
  tab$ucl[!ok] <- NA_real_
  tab$z[!ok]   <- NA_real_
  tab$p[!ok]   <- NA_real_
  
  list(
    dat = dat,
    stage1 = res$stage1,
    summary = tab,
    raw = res
  )
}