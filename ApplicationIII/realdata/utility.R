library(parallel)
library(Rcpp)
library(RcppArmadillo)
source("ApplicationIII/utility_cpp.R")
source("ApplicationIII/utility_sabha.R")
#### simu_fun
simu_fun <- function(alpha_use) {
  alpha_bc <- alpha_use / (1 + alpha_use)

  ################## BH method ##################
  p_adj <- p.adjust(pvalue_vec, method = "BH")
  index_select <- which(p_adj <= alpha_use)
  bh_res <- length(index_select)

  # ################## BC method ##################
  # bc_fit <- bc_fun_cpp(pvalue_vec, alpha_use)
  # index_select <- bc_fit$index_select + 1
  # bc_res <- length(index_select)
  
  ################## IHW ##################
  ihw_res <- try(IHW::ihw(pvalue_vec, x_covariate, alpha_use
  ), silent = TRUE)
  if (class(ihw_res) == "try-error") {
    ihw_res <- 0
  } else {
    ihw_res <- sum(IHW::adj_pvalues(ihw_res) <= alpha_use)
  }
  # ################## IHW Storey ##################
  # ihw_storey <- try(IHWStatsPaper::ihw_bh(pvalue_vec, x_covariate, alpha_use,
  #                                         Storey = TRUE
  # ), silent = TRUE)
  # if (class(ihw_storey) == "try-error") {
  #   ihw_storey_res <- 0
  # } else {
  #   ihw_storey_res <- sum(ihw_storey)
  # }

  # ################## IHW betamix ##################
  # ihw_betamix <- try(IHWStatsPaper::ihw_betamix_censored(pvalue_vec, x_covariate, alpha_use,
  #                                                        Storey = TRUE
  # ), silent = TRUE)
  # if (class(ihw_betamix) == "try-error") {
  #   ihw_betamix_res <- 0
  # } else {
  #   ihw_betamix_res <- sum(ihw_betamix)
  # }

  ################## AdaPT ##################
  ## load data
  x_mat <- as.matrix(x_covariate)
  feature_no <- length(pvalue_vec)
  adapt_padj <- rep(1, feature_no)
  x_df <- data.frame(x_mat)
  colnames(x_df) <- paste0("cov", seq_len(ncol(x_mat)))
  formulas <- paste0(".")
  adapt_fit <- try(adaptMT::adapt_glm(
    x = x_df, pvals = pvalue_vec, pi_formulas = formulas,
    mu_formulas = formulas, alphas = alpha_use,
    verbose = list(print = FALSE, fit = FALSE, ms = FALSE)
  ), silent = TRUE)
  if (class(adapt_fit)[1] == "try-error") {
    adapt_res <- 0
  } else {
    adapt_res <- adapt_fit$nrejs
  }

  ################## SABHA ##################
  n <- length(pvalue_vec)
  index0 <- order(x_covariate)
  index1 <- order(index0)
  fdr <- rep(1, n)
  pvals <- pvalue_vec[index0]
  qhat <- Solve_q_step(pvals, 0.5, 0.1)
  fdr[SABHA_method(pvals, qhat, alpha_use, 0.5)] <- alpha_use / 2
  fdr <- fdr[index1]
  index_select <- which(fdr < alpha_use)
  sabha_res <- length(index_select)

  ################## cross_fit average ##################
  #### try 1
  ebhgbh_res <- ebh_gbc(
    pvalue_vec = pvalue_vec, pi0_var = x_covariate, f1_var = x_covariate,
    n_try = 1, n_split = 2, alpha_ebh = alpha_use, alpha_bc = alpha_bc,
    weight_bool = TRUE
  )
  index_select <- ebhgbh_res$index_select
  ebh_res <- length(index_select)

  ################## save results ##################
  # return(c(
  #   bh_res, ihw_storey_res, ihw_betamix_res,
  #   adapt_res, sabha_res, ebh_res
  # ))
  return(c(
    bh_res, ihw_res, adapt_res, sabha_res, ebh_res
  ))
}
