library(parallel)
library(doParallel)
library(doRNG)
#### define parameters
n_cores <- 10
n <- 3e3
alpha <- 0.1
alpha_ebh <- 0.1
alpha_bc <- alpha_ebh / (1 + alpha_ebh)
a0_vec <- c(3.5, 2.5, 1.5)
a1_vec <- c(1.5, 2, 2.5)
signa_sten_vec <- seq(1.9, 3.85, by = 0.15)
n_simu <- 100

#### simulations ####
cl <- makeCluster(n_cores)
registerDoParallel(cl)
# Load the utility functions once in each worker.
clusterEvalQ(cl, {
  source("ApplicationIII/utility_cpp.R")
  source("ApplicationIII/utility_sabha.R")
  NULL
})
dir.create("ApplicationIII/simu1/results", recursive = TRUE, showWarnings = FALSE)
for (a0 in a0_vec) {
  for (a1 in a1_vec) {
    for (signa_sten in signa_sten_vec) {
      set.seed(1)
      print(c(a0, a1, signa_sten))
      out_res <- foreach(iter_simu = seq_len(n_simu),
                         .noexport = c("ebh_gbc", "Solve_q_step", "SABHA_method")) %dorng% {
        ## generate datasets
        x_covariate <- rnorm(n)
        u_covariate <- a0 + a1 * x_covariate
        pi0_true <- exp(u_covariate) / (1 + exp(u_covariate))
        indicator_true <- sapply(pi0_true, function(pi0) {
          sample(c(0, 1), size = 1, prob = c(pi0, 1 - pi0))
        })
        index_alter <- which(indicator_true == 1)
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

        ################## IHW ##################
        ihw_fit <- try(IHW::ihw(pvalue_vec, x_covariate, alpha),
          silent = TRUE
        )
        if (class(ihw_fit) == "try-error") {
          ihw_res <- numeric(2)
        } else {
          ihw_res <- numeric(2)
          index_select <- which(IHW::adj_pvalues(ihw_fit) <= alpha)
          if (length(index_select) > 0) {
            ihw_res <- c(
              length(intersect(index_select, index_alter)) / length(index_alter),
              length(setdiff(index_select, index_alter)) / length(index_select)
            )
          }
        }

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
          mu_formulas = formulas, alphas = 0.1,
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
          pvalue_vec = pvalue_vec, pi0_var = x_covariate, f1_var = x_covariate,
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
        c(bh_res, ihw_res, adapt_res, sabha_res, ebh_res)
      }
      res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
      saveRDS(res_mat, paste0(
        "ApplicationIII/simu1/results/a0", a0, "a1", a1, "sign", signa_sten, ".rds"
      ))
      print(round(colMeans(res_mat), 4))
    }
  }
}
stopCluster(cl)
registerDoSEQ()
