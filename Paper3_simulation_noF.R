sim_binary_iv <- function(n, p_z) {
  k <- length(p_z)
  
  Z_mat <- sapply(seq_len(k), function(j) {
    rbinom(n, size = 1, prob = p_z[j])
  })
  
  Z_mat <- matrix(Z_mat, nrow = n, ncol = k)
  colnames(Z_mat) <- paste0("Z", seq_len(k))
  
  Z_mat
}


sim_once <- function(n,
                     a1,
                     b1,
                     c_x,
                     c_y,
                     p_z = 0.2,
                     sigma_u = 0.5,
                     a0 = 0,
                     b0 = 0,
                     scen = NA_character_,
                     strength_id = NA_character_,
                     strength_label = NA_character_,
                     a1_value = NA_real_,
                     maxit = 1000) {
  
  if (length(a1) < 1) {
    stop("a1 must contain at least one instrument-strength value.")
  }
  
  k <- length(a1)
  
  if (length(p_z) == 1) {
    p_z <- rep(p_z, k)
  }
  
  if (length(p_z) != k) {
    stop("Length of p_z must be 1 or equal to length(a1).")
  }
  
  Z_mat <- sim_binary_iv(
    n = n,
    p_z = p_z
  )
  
  U <- rnorm(n, mean = 0, sd = sigma_u)
  
  lin_x <- a0 + as.numeric(Z_mat %*% a1) + c_x * U
  p_x <- trim_prob(expit(lin_x))
  X <- rbinom(n, size = 1, prob = p_x)
  
  lin_y <- b0 + b1 * X + c_y * U
  p_y <- trim_prob(expit(lin_y))
  Y <- rbinom(n, size = 1, prob = p_y)
  
  dat <- data.frame(
    X = X,
    Y = Y,
    U = U,
    Z_mat,
    check.names = FALSE
  )
  
  fit <- tryCatch(
    MRbinary(
      dat = dat[, c(colnames(Z_mat), "X", "Y")],
      z = colnames(Z_mat),
      x = "X",
      y = "Y",
      maxit = maxit
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || is.null(fit$summary)) {
    return(data.frame(
      method = c("2SPS", "2SRI", "GMM", "IV-MVB"),
      est = NA_real_,
      se = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      z = NA_real_,
      p = NA_real_,
      conv = 0,
      b1 = b1,
      c_x = c_x,
      c_y = c_y,
      strength_id = strength_id,
      strength_label = strength_label,
      a1_value = a1_value,
      stringsAsFactors = FALSE
    ))
  }
  
  tab <- fit$summary
  
  tab$b1 <- b1
  tab$c_x <- c_x
  tab$c_y <- c_y
  tab$strength_id <- strength_id
  tab$strength_label <- strength_label
  tab$a1_value <- a1_value
  
  tab
}


run_sim <- function(n_sim,
                    n,
                    a1,
                    b1,
                    c_x,
                    c_y,
                    p_z,
                    sigma_u = 0.5,
                    scen = NA_character_,
                    strength_id = NA_character_,
                    strength_label = NA_character_,
                    a1_value = NA_real_,
                    methods = c("2SPS", "2SRI", "GMM", "IV-MVB"),
                    seed = NULL,
                    maxit = 1000) {
  
  if (!is.null(seed)) set.seed(seed)
  
  res_list <- vector("list", n_sim)
  
  make_fail <- function(i) {
    data.frame(
      sim = i,
      scen = scen,
      n = n,
      b1 = b1,
      c_x = c_x,
      c_y = c_y,
      strength_id = strength_id,
      strength_label = strength_label,
      a1_value = a1_value,
      method = methods,
      est = NA_real_,
      se = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      z = NA_real_,
      p = NA_real_,
      conv = 0,
      stringsAsFactors = FALSE
    )
  }
  
  for (i in seq_len(n_sim)) {
    
    out <- tryCatch(
      sim_once(
        n = n,
        a1 = a1,
        b1 = b1,
        c_x = c_x,
        c_y = c_y,
        p_z = p_z,
        sigma_u = sigma_u,
        scen = scen,
        strength_id = strength_id,
        strength_label = strength_label,
        a1_value = a1_value,
        maxit = maxit
      ),
      error = function(e) NULL
    )
    
    if (is.null(out)) {
      res_list[[i]] <- make_fail(i)
      next
    }
    
    tab <- out
    
    if (!"method" %in% names(tab) && "Method" %in% names(tab)) {
      names(tab)[names(tab) == "Method"] <- "method"
    }
    if (!"est" %in% names(tab) && "Estimate" %in% names(tab)) {
      names(tab)[names(tab) == "Estimate"] <- "est"
    }
    if (!"se" %in% names(tab) && "SE" %in% names(tab)) {
      names(tab)[names(tab) == "SE"] <- "se"
    }
    if (!"lcl" %in% names(tab) && "LCL" %in% names(tab)) {
      names(tab)[names(tab) == "LCL"] <- "lcl"
    }
    if (!"ucl" %in% names(tab) && "UCL" %in% names(tab)) {
      names(tab)[names(tab) == "UCL"] <- "ucl"
    }
    
    if (!"z" %in% names(tab)) tab$z <- NA_real_
    if (!"p" %in% names(tab)) tab$p <- NA_real_
    if (!"conv" %in% names(tab)) tab$conv <- NA_real_
    
    tab$sim <- i
    tab$scen <- scen
    tab$n <- n
    tab$b1 <- b1
    tab$c_x <- c_x
    tab$c_y <- c_y
    tab$strength_id <- strength_id
    tab$strength_label <- strength_label
    tab$a1_value <- a1_value
    
    miss <- setdiff(methods, tab$method)
    
    if (length(miss) > 0) {
      tab_miss <- data.frame(
        sim = i,
        scen = scen,
        n = n,
        b1 = b1,
        c_x = c_x,
        c_y = c_y,
        strength_id = strength_id,
        strength_label = strength_label,
        a1_value = a1_value,
        method = miss,
        est = NA_real_,
        se = NA_real_,
        lcl = NA_real_,
        ucl = NA_real_,
        z = NA_real_,
        p = NA_real_,
        conv = 0,
        stringsAsFactors = FALSE
      )
      
      tab <- dplyr::bind_rows(tab, tab_miss)
    }
    
    keep <- c(
      "sim",
      "scen",
      "n",
      "b1",
      "c_x",
      "c_y",
      "strength_id",
      "strength_label",
      "a1_value",
      "method",
      "est",
      "se",
      "lcl",
      "ucl",
      "z",
      "p",
      "conv"
    )
    
    for (v in keep) {
      if (!v %in% names(tab)) tab[[v]] <- NA
    }
    
    tab <- tab[, keep, drop = FALSE]
    tab <- tab[match(methods, tab$method), , drop = FALSE]
    
    res_list[[i]] <- tab
  }
  
  res <- dplyr::bind_rows(res_list)
  rownames(res) <- NULL
  
  res
}