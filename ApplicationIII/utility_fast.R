#### get initial value
solvek_fdr <- function(pi0, pvals, pvals_cutoff = 1e-15) {
  pvals[pvals <= pvals_cutoff] <- pvals_cutoff
  pvals[pvals >= 1 - pvals_cutoff] <- 1 - pvals_cutoff

  pi1 <- 1 - pi0

  fun <- function(k) {
    pk <- pvals^(-k)
    sum(-pi1 * pk * (1 + (1 - k) * log(pvals)) / (pi0 + pi1 * (1 - k) * pk))
  }

  err <- try(
    {
      k <- uniroot(fun, c(0.01, 0.99))$root
    },
    silent = TRUE
  )

  if (class(err) == "try-error") {
    k <- 0.75
  }
  return(k)
}

#### calculate weight
calculate_weights <- function(pvals, pi0_var, f1_var, EM_paras) {
  ## initialization
  n <- length(pvals)
  len_pi0 <- ncol(pi0_var)
  len_f1 <- ncol(f1_var)
  pi0_est <- EM_paras$pi0_init
  k_init <- EM_paras$k_init
  tol <- EM_paras$tol
  iterlim <- EM_paras$iterlim
  nlm_iter <- EM_paras$nlm_iter
  # No covariate effect - initialization
  theta1 <- c(binomial()$linkfun(pi0_est), rep(0, len_pi0 - 1))
  beta1 <- c(binomial()$linkfun(k_init), rep(0, len_f1 - 1))

  iter <- 0
  lpval <- log(pvals)

  func2 <- function(x, len1, len2, q0) {
    theta <- x[1:len1]
    beta <- x[(len1 + 1):(len1 + len2)]
    q1 <- 1 - q0
    exp_eta_theta <- exp(as.vector(pi0_var %*% theta))
    pi0_temp <- (exp_eta_theta / (1 + exp_eta_theta))
    exp_eta_beta <- exp(as.vector(f1_var %*% beta))
    k_temp <- (exp_eta_beta / (1 + exp_eta_beta))
    res <- -sum(q0 * log(pi0_temp) + q1 * log(1 - pi0_temp) + q1 * (-k_temp * lpval + log(1 - k_temp)))
    return(res)
  }

  loglik0 <- 0
  while (iter < iterlim) {
    # E step
    exp_eta_theta <- exp(as.vector(pi0_var %*% theta1))
    exp_eta_beta <- exp(as.vector(f1_var %*% beta1))
    pi0 <- exp_eta_theta / (1 + exp_eta_theta)
    k <- exp_eta_beta / (1 + exp_eta_beta)
    f1 <- (1 - k) * pvals^(-k)
    f01 <- (1 - pi0) * f1 + pi0
    q0 <- pi0 / f01

    # M step
    st <- c(theta1, beta1)
    obj <- suppressWarnings(nlm(func2, st, len1 = len_pi0, len2 = len_f1, q0 = q0, iterlim = nlm_iter))

    # Update new likelihood
    theta2 <- obj$estimate[1:len_pi0]
    beta2 <- obj$estimate[(len_pi0 + 1):(len_pi0 + len_f1)]

    loglik1 <- sum(log(f01))

    if (abs((loglik1 - loglik0) / loglik0) < tol) {
      break
    } else {
      theta1 <- theta2
      beta1 <- beta2
      iter <- iter + 1
      loglik0 <- loglik1
    }
  }

  result <- list(
    pi0 = pi0, k = k, pi0_coef = theta1, k_coef = beta1,
    EM_paras = EM_paras, loglik = loglik1, pvals = pvals,
    EM_iter = iter
  )
  return(result)
}

gbh_fun <- function(n_test, lfdr_fit, invlfdr, alpha_bh) {
  ## initialization
  evalue_vec <- numeric(n_test)
  hat_invlfdr <- numeric(n_test)
  ## Do GBH
  ordered_lfdr <- sort(lfdr_fit, decreasing = TRUE)
  for (iter_lfdr in seq_len(n_test + 1)) {
    if (iter_lfdr == n_test + 1) {
      index_select <- NULL
      sum_fit <- 0
    } else {
      lfdr_tmp <- ordered_lfdr[iter_lfdr]
      n_reject <- n_test + 1 - iter_lfdr
      hat_invlfdr <- invlfdr(lfdr_tmp)
      sum_fit <- sum(hat_invlfdr)
      thres <- sum_fit / n_reject
      if (thres <= alpha_bh) {
        index_select <- which(lfdr_fit <= lfdr_tmp)
        e_tmp <- 1 / sum_fit
        evalue_vec[index_select] <- e_tmp
        break
      }
    }
  }
  return(list(
    evalue_vec = evalue_vec, hat_invlfdr = hat_invlfdr,
    index_select = index_select, sum_fit = sum_fit
  ))
}

modelfit_fun <- function(pvalue_train, pi0_var_train, f1_var_train, EM_paras,
                         pvalue_test, pi0_var_test, f1_var_test, alpha_bh) {
  n_test <- length(pvalue_test)
  ## fit the model
  res_fit <- calculate_weights(
    pvals = pvalue_train, pi0_var = pi0_var_train,
    f1_var = f1_var_train, EM_paras = EM_paras
  )
  ## pi0
  pi0_coef <- res_fit$pi0_coef
  exp_eta_theta <- exp(as.vector(pi0_var_test %*% pi0_coef))
  pi0_fit <- exp_eta_theta / (1 + exp_eta_theta)
  pi0_fit[pi0_fit < 0.1] <- 0.1
  pi0_fit[pi0_fit > 1 - 1e-5] <- 1 - 1e-5
  ## kappa
  k_coef <- res_fit$k_coef
  exp_eta_beta <- exp(as.vector(f1_var_test %*% k_coef))
  kappa_fit <- exp_eta_beta / (1 + exp_eta_beta)
  ## calculate lfdr
  f1_fit <- (1 - kappa_fit) * pvalue_test^(-kappa_fit)
  lfdr_fit <- pi0_fit / (pi0_fit + (1 - pi0_fit) * f1_fit)
  ## define inv local FDR function
  invlfdr <- function(t) {
    num_tmp <- pi0_fit * (1 - t)
    dom_tmp <- (1 - pi0_fit) * t * (1 - kappa_fit)
    res <- (num_tmp / dom_tmp)^(-1 / kappa_fit)
    return(res)
  }
  gbh_fit <- gbh_fun(n_test, lfdr_fit, invlfdr, alpha_bh)
  evalue_vec <- gbh_fit$evalue_vec
  sum_fit <- gbh_fit$sum_fit
  index_select <- gbh_fit$index_select + 1

  # ## save inverse lfdr for weight
  hat_invlfdr <- numeric(n_test)
  ## calculate weight vector
  for (iter_lfdr_tmp in seq_len(n_test)) {
    lfdr_fit_tmp <- lfdr_fit
    lfdr_fit_tmp[iter_lfdr_tmp] <- 0
    gbh_fit_tmp <- gbh_fun(n_test, lfdr_fit_tmp, invlfdr, alpha_bh)
    hat_invlfdr[iter_lfdr_tmp] <- gbh_fit_tmp$hat_invlfdr[iter_lfdr_tmp]
  }
  ## return results
  return(list(
    hat_invlfdr = hat_invlfdr, sum_fit = sum_fit,
    index_select = index_select, evalue_vec = evalue_vec
  ))
}

## cross fit function
crossfit_fun <- function(pvalue_vec, pi0_var, f1_var, EM_paras, split_groups,
                         n_split = 10, alpha_bh = 0.05) {
  ## Shuffle the index vector randomly
  n <- length(pvalue_vec)
  ## cross fit
  evalue_vec <- numeric(n)
  evalue_vec_list <- list()
  sum_fit_list <- list()
  hat_invlfdr_list <- list()
  for (iter_split in seq_len(n_split)) {
    ## split data
    # fitted data
    index_test <- split_groups[[iter_split]]
    pvalue_test <- pvalue_vec[index_test]
    pi0_var_test <- pi0_var[index_test, , drop = FALSE]
    f1_var_test <- f1_var[index_test, , drop = FALSE]
    # trainning data
    index_train <- unlist(split_groups[-iter_split])
    pvalue_train <- pvalue_vec[index_train]
    pi0_var_train <- pi0_var[index_train, , drop = FALSE]
    f1_var_train <- f1_var[index_train, , drop = FALSE]
    ## fit the model
    model_fit <- modelfit_fun(
      pvalue_train, pi0_var_train, f1_var_train, EM_paras,
      pvalue_test, pi0_var_test, f1_var_test, alpha_bh
    )
    evalue_vec_list[[iter_split]] <- model_fit$evalue_vec
    sum_fit_list[[iter_split]] <- model_fit$sum_fit
    hat_invlfdr_list[[iter_split]] <- model_fit$hat_invlfdr
  }
  for (iter_split in seq_len(n_split)) {
    if (sum(evalue_vec_list[[iter_split]]) > 0) {
      evalue_vec[split_groups[[iter_split]]] <- n * evalue_vec_list[[iter_split]] *
        sum_fit_list[[iter_split]] / (sum_fit_list[[iter_split]] + sum(unlist(hat_invlfdr_list[-iter_split])))
    }
  }
  return(evalue_vec)
}

## main function
ebh_gbh <- function(pvalue_vec, pi0_var = NULL, f1_var = NULL, n_sample = 10, n_split = 10,
                    alpha_ebh = 0.05, alpha_bh = 0.025, pvals_cutoff = 1e-15,
                    EM_paras = list(iterlim = 50, tol = 1e-5, k_init = NULL, pi0_init = NULL, nlm_iter = 5)) {
  #### initialization
  n <- length(pvalue_vec)
  ## EM parameters
  if (is.null(EM_paras$pi0_init)) {
    pi0_est <- qvalue::qvalue(pvalue_vec, pi0.method = "bootstrap")$pi0
    if (pi0_est >= 0.99) pi0_est <- 0.99
    EM_paras$pi0_init <- pi0_est
  } else {
    pi0_est <- EM_paras$pi0_init
  }
  if (is.null(EM_paras$k_init)) {
    k_init <- solvek_fdr(pi0_est, pvalue_vec, pvals_cutoff)
    EM_paras$k_init <- k_init
  } else {
    k_init <- EM_paras$k_init
  }
  if (is.null(EM_paras$iterlim)) {
    iterlim <- 50
    EM_paras$iterlim <- iterlim
  } else {
    iterlim <- EM_paras$iterlim
  }
  if (is.null(EM_paras$nlm_iter)) {
    nlm_iter <- 5
    EM_paras$nlm_iter <- nlm_iter
  } else {
    nlm_iter <- EM_paras$nlm_iter
  }
  if (is.null(EM_paras$tol)) {
    tol <- 1e-5
    EM_paras$tol <- tol
  } else {
    tol <- EM_paras$tol
  }
  ## covariate information
  if (is.null(pi0_var)) {
    pi0_var <- matrix(rep(1, n))
  } else {
    pi0_var <- cbind(rep(1, n), pi0_var)
  }
  if (is.null(f1_var)) {
    f1_var <- matrix(rep(1, n))
  } else {
    f1_var <- cbind(rep(1, n), f1_var)
  }

  ## multiple try
  index_list <- list()
  for (iter_try in seq_len(n_sample)) {
    shuffled_index <- sample(seq_len(n))
    index_list[[iter_try]] <- split(shuffled_index, ceiling(seq_len(n) * n_split / n))
  }

  evalue_mat <- matrix(NA, nrow = n_sample, ncol = n)
  for (iter_try in seq_len(n_sample)) {
    split_groups <- index_list[[iter_try]]
    evalue_mat[iter_try, ] <- crossfit_fun(pvalue_vec, pi0_var, f1_var, EM_paras, split_groups,
      n_split = n_split, alpha_bh = alpha_bh
    )
  }
  evalue_ave <- colMeans(evalue_mat)

  ## e-BH
  orderd_evalue <- sort(evalue_ave, decreasing = FALSE)
  for (iter in seq_len(n + 1)) {
    if (iter == n + 1) {
      index_select <- NULL
    } else {
      e_tmp <- orderd_evalue[iter]
      thres <- n / (alpha_ebh * (n + 1 - iter))
      if (e_tmp >= thres) {
        index_select <- which(evalue_ave >= e_tmp)
        break
      }
    }
  }
  return(list(
    index_select = index_select, evalue_ave = evalue_ave,
    evalue_mat = evalue_mat
  ))
}
