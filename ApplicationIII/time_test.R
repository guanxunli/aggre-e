set.seed(1)
source("stru_ada_bc/utility_sabha.R")
source("stru_ada_bc/utility_cpp.R")
#### define parameters
n <- 3e3
alpha <- 0.1
alpha_ebh <- 0.1
alpha_bc <- alpha_ebh / (1 + alpha_ebh)
a0_vec <- c(3.5, 2.5, 1.5)
a1_vec <- c(1.5, 2, 2.5)
af <- 1
signa_sten_vec <- seq(2.5, 3.4, by = 0.15)
n_simu <- 1e2

#### simulation function ####
simu_fun <- function(a0, a1, af, signa_sten) {
  time_res <- numeric(7)
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
  time1 <- Sys.time()
  p_adj <- p.adjust(pvalue_vec, method = "BH")
  index_select <- which(p_adj <= alpha)
  bh_res <- numeric(2)
  if (length(index_select) > 0) {
    bh_res <- c(
      length(intersect(index_select, index_alter)) / length(index_alter),
      length(setdiff(index_select, index_alter)) / length(index_select)
    )
  }
  time_res[1] <- as.numeric(Sys.time() - time1, units = "secs")
  
  ################## BC method ##################
  time1 <- Sys.time()
  bc_fit <- bc_fun_cpp(pvalue_vec, alpha)
  index_select <- bc_fit$index_select + 1
  bc_res <- numeric(2)
  if (length(index_select) > 0) {
    bc_res <- c(
      length(intersect(index_select, index_alter)) / length(index_alter),
      length(setdiff(index_select, index_alter)) / length(index_select)
    )
  }
  time_res[2] <- as.numeric(Sys.time() - time1, units = "secs")
  
  ################## IHW Storey ##################
  tiem1 <- Sys.time()
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
  time_res[3] <- as.numeric(Sys.time() - time1, units = "secs")
  ################## IHW betamix ##################
  time1 <- Sys.time()
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
  time_res[4] <- as.numeric(Sys.time() - time1, units = "secs")
  ################## AdaPT ##################
  ## load data
  time1 <- Sys.time()
  feature_no <- length(pvalue_vec)
  adapt_padj <- rep(1, feature_no)
  x_df <- data.frame(x1 = x_covariate, x2 = x_covariate_f)
  pi_formulas <- paste0("x1")
  mu_formulas <- paste0("x2")
  adapt_fit <- try(adaptMT::adapt_glm(
    x = x_df, pvals = pvalue_vec, pi_formulas = pi_formulas,
    mu_formulas = mu_formulas, alphas = 0.1,
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
  time_res[5] <- as.numeric(Sys.time() - time1, units = "secs")
  ################## SABHA ##################
  time1 <- Sys.time()
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
  time_res[6] <- as.numeric(Sys.time() - time1, units = "secs")
  ################## cross_fit average ##################
  time1 <- Sys.time()
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
  time_res[7] <- as.numeric(Sys.time() - time1, units = "secs")
  
  ################## return results ##################
  return(time_res)
}

#### simulations ####
a0 <- 2.5
a1 <- 2
signa_sten <- 3
set.seed(1)
print(c(a0, a1, signa_sten))
time_mat <- matrix(NA, nrow = n_simu, ncol = 7)
for (iter_simu in seq_len(n_simu)) {
  if (iter_simu %% 10 == 0) print(iter_simu)
  time_mat[iter_simu, ] <- simu_fun(a0, a1, af, signa_sten)
}
print(round(colMeans(time_mat), 4))