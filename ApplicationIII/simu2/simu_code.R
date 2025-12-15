library(parallel)
source("ApplicationIII/utility_sabha.R")
source("ApplicationIII/utility_cpp.R")
#### define parameters
n <- 3e3
alpha <- 0.1
alpha_ebh <- 0.1
alpha_bc <- alpha_ebh / (1 + alpha_ebh)
a0_vec <- c(3.5, 2.5, 1.5)
a1_vec <- c(1.5, 2, 2.5)
af <- 0.5
signa_sten_vec <- seq(2.5, 3.4, by = 0.15)
n_simu <- 500
#### simulation function ####
simu_fun <- function(a0, a1, af, signa_sten) {
  ## generate datasets
  # pi0
  x_covariate <- rnorm(n)
  a_use <- a0 + a1 * x_covariate
  pi0_true <- exp(a_use) / (1 + exp(a_use))
  indicator_true <- sapply(pi0_true, function(pi0) {
    sample(c(0, 1), size = 1, prob = c(pi0, 1 - pi0))
  })
  index_alter <- which(indicator_true == 1)
  # fi1
  x_covariate_f <- rnorm(n)
  signa_sten_scale <- 2 * exp(af * x_covariate_f) / (1 + exp(af * x_covariate_f))
  signa_sten <- signa_sten_scale * signa_sten
  # p-value
  z_vec <- rnorm(n, mean = signa_sten * indicator_true, sd = 1)
  pvalue_vec <- 1 - pnorm(z_vec)

  ################## BH method ##################
  p_adj <- p.adjust(pvalue_vec, method = "BH")
  index_select <- which(p_adj <= alpha)
  bh_res <- numeric(2)
  if (length(index_select) > 0) {
    bh_res <- c(
      length(intersect(index_select, index_alter)) / length(index_alter),
      length(setdiff(index_select, index_alter)) / length(index_select)
    )
  }

  ################## IHW Storey ##################
  x_mat <- cbind(x_covariate, x_covariate_f)
  ihw_storey <- try(IHWStatsPaper::ihw_bh(pvalue_vec, x_mat, alpha,
    Storey = TRUE
  ), silent = TRUE)
  if (class(ihw_storey) == "try-error") {
    ihw_storey_res <- numeric(2)
  } else {
    ihw_storey_res <- numeric(2)
    index_select <- which(ihw_storey == TRUE)
    if (length(index_select) > 0) {
      ihw_storey_res <- c(
        length(intersect(index_select, index_alter)) / length(index_alter),
        length(setdiff(index_select, index_alter)) / length(index_select)
      )
    }
  }

  ################## IHW betamix ##################
  ihw_betamix <- try(IHWStatsPaper::ihw_betamix_censored(pvalue_vec, x_mat, alpha,
    Storey = TRUE
  ), silent = TRUE)
  if (class(ihw_betamix) == "try-error") {
    ihw_betamix_res <- numeric(2)
  } else {
    ihw_betamix_res <- numeric(2)
    index_select <- which(ihw_betamix == TRUE)
    if (length(index_select) > 0) {
      ihw_betamix_res <- c(
        length(intersect(index_select, index_alter)) / length(index_alter),
        length(setdiff(index_select, index_alter)) / length(index_select)
      )
    }
  }

  ################## AdaPT ##################
  ## load data
  feature_no <- length(pvalue_vec)
  adapt_padj <- rep(1, feature_no)
  x_df <- data.frame(x1 = x_covariate, x2 = x_covariate_f)
  pi_formulas <- paste0("x1")
  mu_formulas <- paste0("x2")
  adapt_fit <- try(adaptMT::adapt_glm(
    x = x_df, pvals = pvalue_vec, pi_formulas = pi_formulas,
    mu_formulas = mu_formulas, alphas = alpha,
    verbose = list(print = FALSE, fit = FALSE, ms = FALSE)
  ), silent = TRUE)
  if (class(adapt_fit)[1] == "try-error") {
    adapt_res <- numeric(2)
  } else {
    adapt_res <- numeric(2)
    index_select <- adapt_fit$rejs
    if (length(index_select) > 0) {
      adapt_res <- c(
        length(intersect(index_select, index_alter)) / length(index_alter),
        length(setdiff(index_select, index_alter)) / length(index_select)
      )
    }
  }

  ################## SABHA ##################
  index0 <- order(x_covariate)
  index1 <- order(index0)
  fdr <- rep(1, n)
  pvals <- pvalue_vec[index0]
  qhat <- Solve_q_step(pvals, 0.5, 0.1)
  fdr[SABHA_method(pvals, qhat, alpha, 0.5)] <- alpha / 2
  fdr <- fdr[index1]
  index_select <- which(fdr < alpha)
  sabha_res <- numeric(2)
  if (length(index_select) > 0) {
    sabha_res <- c(
      length(intersect(index_select, index_alter)) / length(index_alter),
      length(setdiff(index_select, index_alter)) / length(index_select)
    )
  }

  ################## cross_fit average ##################
  #### try 1
  ebhgbh_res <- ebh_gbc(
    pvalue_vec = pvalue_vec, pi0_var = x_covariate, f1_var = x_covariate_f,
    n_try = 1, n_split = 2, alpha_ebh = alpha_ebh, alpha_bc = alpha_bc,
    weight_bool = TRUE
  )
  index_select <- ebhgbh_res$index_select
  ebh_res <- numeric(2)
  if (length(index_select) > 0) {
    ebh_res <- c(
      length(intersect(index_select, index_alter)) / length(index_alter),
      length(setdiff(index_select, index_alter)) / length(index_select)
    )
  }

  ################## return results ##################
  return(c(
    bh_res, bc_res, ihw_storey_res, ihw_betamix_res, adapt_res, sabha_res, ebh_res
  ))
}

#### simulations ####
for (a0 in a0_vec) {
  for (a1 in a1_vec) {
    for (signa_sten in signa_sten_vec) {
      set.seed(1)
      print(c(a0, a1, signa_sten))
      out_res <- mclapply(seq_len(n_simu), function(i) {
        return(simu_fun(a0, a1, af, signa_sten))
      }, mc.cores = 50)
      res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
      saveRDS(res_mat, paste0(
        "ApplicationIII/simu2/results/a0", a0, "a1", a1, "sign", signa_sten, ".rds"
      ))
      print(round(colMeans(res_mat), 4))
    }
  }
}
