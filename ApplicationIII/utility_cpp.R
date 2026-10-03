library(Rcpp)
library(RcppArmadillo)
sourceCpp("ApplicationIII/utilities.cpp")
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

## main function
ebh_gbc <- function(pvalue_vec, pi0_var = NULL, f1_var = NULL, n_try = 10, n_split = 10,
                    alpha_ebh = 0.05, alpha_bc = 0.025, pvals_cutoff = 1e-15,
                    EM_paras = list(iterlim = 50, tol = 1e-5, k_init = NULL, pi0_init = NULL, nlm_iter = 5),
                    weight_bool = TRUE, parallel = FALSE, n_cores = 1) {
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
  ## gbh-ebh
  if (parallel) {
    library(parallel)
    out_res <- mclapply(seq_len(n_try), function(iter_try) {
      split_groups <- shuffle_fun(n, n_split)
      evalue_tmp <- crossfit_fun(
        pvalue_vec, pi0_var, f1_var, alpha_bc, pi0_est,
        k_init, tol, iterlim, nlm_iter, n_split,
        split_groups, weight_bool
      )
      return(evalue_tmp)
    }, mc.cores = min(n_cores, n_try))
    evalue_mat <- matrix(unlist(out_res), nrow = n_try, byrow = TRUE)
    evalue_ave <- colMeans(evalue_mat)
    index_select <- ebh_fun_cpp(evalue_ave, alpha_ebh)
    outres <- list(
      index_select = index_select,
      evalue_mat = evalue_mat,
      evalue_ave = evalue_ave
    )
    outres$index_select <- outres$index_select + 1
  } else {
    outres <- suppressWarnings(ebh_gbc_cpp(
      pvalue_vec, pi0_var, f1_var, n_try, n_split,
      alpha_ebh, alpha_bc, pi0_est, k_init, tol, iterlim, nlm_iter,
      weight_bool
    ))
    outres$index_select <- outres$index_select + 1
  }
  return(outres)
}
