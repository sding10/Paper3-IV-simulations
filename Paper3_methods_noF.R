## Packages
library(MASS)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gt)

## Data preparation
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


## Instrumental strength
fit_first_stage <- function(dat, z, x) {
  fml <- as.formula(
    paste(x, "~", paste(z, collapse = " + "))
  )
  
  fit <- glm(
    fml,
    data = dat,
    family = binomial()
  )
  
  fit
}


## 2SPS
est_2sps <- function(dat, z, x, y, alpha = 0.05, eps = 1e-8) {
  out_fail <- c(
    method = "2SPS",
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
  
  fit_2 <- tryCatch(
    glm(
      as.formula(paste(y, "~ x_hat")),
      data = dat,
      family = binomial()
    ),
    error = function(e) NULL
  )
  if (is.null(fit_2)) return(out_fail)
  
  ct <- summary(fit_2)$coefficients
  if (!("x_hat" %in% rownames(ct))) return(out_fail)
  
  est <- unname(ct["x_hat", "Estimate"])
  se  <- unname(ct["x_hat", "Std. Error"])
  z   <- unname(ct["x_hat", "z value"])
  p   <- unname(ct["x_hat", "Pr(>|z|)"])
  
  if (!is.finite(est) || !is.finite(se) || se <= eps) {
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
    conv = 1
  )
}


## 2SRI
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


## GMM
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


## IV-MVB
est_ivmvb <- function(dat, z, x, y,
                      start = NULL,
                      alpha = 0.05,
                      maxit = 1000,
                      eps = 1e-8) {
  x_obs <- dat[[x]]
  y_obs <- dat[[y]]
  Z_mat <- as.matrix(dat[, z, drop = FALSE])
  k <- ncol(Z_mat)
  
  rho_obs <- suppressWarnings(cor(x_obs, y_obs, use = "complete.obs"))
  
  if (!is.finite(rho_obs)) {
    rho_obs <- NA_real_
  }
  
  if (is.finite(rho_obs)) {
    rho_obs <- max(min(rho_obs, 1 - eps), -1 + eps)
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
    rho_obs = rho_obs
  )
  
  if (!is.finite(rho_obs)) {
    return(out_fail)
  }
  
  nll <- function(par) {
    a0 <- par[1]
    a  <- par[2:(k + 1)]
    b0 <- par[k + 2]
    b1 <- par[k + 3]
    
    pi_x <- trim_prob(expit(as.numeric(a0 + Z_mat %*% a)), eps)
    pi_y <- trim_prob(expit(b0 + b1 * pi_x), eps)
    
    denom <- sqrt(pi_x * (1 - pi_x) * pi_y * (1 - pi_y))
    denom <- pmax(denom, eps)
    theta <- rho_obs / denom
    
    p11 <- pi_x * pi_y * (1 + theta * (1 - pi_x) * (1 - pi_y))
    p10 <- pi_x * (1 - pi_y) * (1 - theta * (1 - pi_x) * pi_y)
    p01 <- (1 - pi_x) * pi_y * (1 - theta * pi_x * (1 - pi_y))
    p00 <- (1 - pi_x) * (1 - pi_y) * (1 + theta * pi_x * pi_y)
    
    if (any(!is.finite(p11)) || any(!is.finite(p10)) ||
        any(!is.finite(p01)) || any(!is.finite(p00))) {
      return(1e12)
    }
    
    if (any(p11 <= 0) || any(p10 <= 0) || any(p01 <= 0) || any(p00 <= 0)) {
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
    rho_obs = rho_obs
  )
}


## Warpper
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
      data.frame(
        method = as.character(o["method"]),
        est    = as.numeric(o["est"]),
        se     = as.numeric(o["se"]),
        lcl    = as.numeric(o["lcl"]),
        ucl    = as.numeric(o["ucl"]),
        z      = as.numeric(o["z"]),
        p      = as.numeric(o["p"]),
        conv   = if ("conv" %in% names(o)) as.numeric(unname(o["conv"])) else 1,
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


