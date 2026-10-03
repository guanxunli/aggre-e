library(parallel)
library(doParallel)
library(doRNG)
## define parameters
n_cores <- 10
alpha <- alpha_ebh <- 0.05
alpha_bh <- alpha_ebh / 2
alpha_bha <- alpha_ebh / (1 + alpha_ebh)
n_simu <- 100 # Total repetitions per signal strength.
seed <- 2026092301
n <- 1000
ratio <- 0.5
sd_a <- 0.3
signa_stre_vec <- c(1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2)
n1 <- n * ratio
index_alter <- seq_len(n1)
dir.create("ApplicationII/results", recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
registerDoParallel(cl)
# Load the C++ functions once in each worker.
clusterEvalQ(cl, {
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp("ApplicationII/utilities.cpp")
  NULL
})
RNGkind("L'Ecuyer-CMRG")
for (iter_sign in seq_along(signa_stre_vec)) {
  signa_stre <- signa_stre_vec[iter_sign]
  mu <- signa_stre # Actual mean; do not multiply by log(n).
  print(signa_stre)
  set.seed(seed + 20000 + iter_sign * 100) # Derive each grid point from the single base seed.
  out_res <- foreach(iter_simu = seq_len(n_simu),
                     .noexport = c("bh_fun_cpp", "bc_fun_cpp", "bhbc_fix", "bhbc_ada",
                                   "fastbhbc_ada", "storey_fun_cpp")) %dorng% {
    ## generate p-values
    x_all <- c(rnorm(n1, mean = mu, sd = sd_a), rnorm(n - n1))
    p_all <- pnorm(x_all, lower.tail = FALSE)
    ## BH
    bh_res <- bh_fun_cpp(p_all, alpha = alpha)
    indexbh <- bh_res$index_select + 1
    bh_res <- numeric(2)
    if (length(indexbh) > 0) {
      bh_res[1] <- length(intersect(indexbh, index_alter)) / length(index_alter)
      bh_res[2] <- length(setdiff(indexbh, index_alter)) / length(indexbh)
    }
    ## BC
    bc_res <- bc_fun_cpp(p_all, alpha = alpha)
    indexbc <- bc_res$index_select + 1
    bc_res <- numeric(2)
    if (length(indexbc) > 0) {
      bc_res[1] <- length(intersect(indexbc, index_alter)) / length(index_alter)
      bc_res[2] <- length(setdiff(indexbc, index_alter)) / length(indexbc)
    }
    ## ave
    bhbc_res <- bhbc_fix(p_all, alpha_bh = alpha_bh, alpha_ebh = alpha_ebh)
    indexbhbc <- bhbc_res$index_select + 1
    bhbc_res <- numeric(2)
    if (length(indexbhbc) > 0) {
      bhbc_res[1] <- length(intersect(indexbhbc, index_alter)) / length(index_alter)
      bhbc_res[2] <- length(setdiff(indexbhbc, index_alter)) / length(indexbhbc)
    }
    ## ada
    bhbcada_res <- bhbc_ada(p_all, alpha_bh = alpha_bha, alpha_ebh = alpha_ebh)
    indexbhbcada <- bhbcada_res$index_select + 1
    bhbcada_res <- numeric(2)
    if (length(indexbhbcada) > 0) {
      bhbcada_res[1] <- length(intersect(indexbhbcada, index_alter)) / length(index_alter)
      bhbcada_res[2] <- length(setdiff(indexbhbcada, index_alter)) / length(indexbhbcada)
    }
    ## fast ada
    fastbhbcada_res <- fastbhbc_ada(p_all, alpha_bh = alpha_bha, alpha_ebh = alpha_ebh)
    indexbhbcadaf <- fastbhbcada_res$index_select + 1
    fastbhbcada_res <- numeric(2)
    if (length(indexbhbcadaf) > 0) {
      fastbhbcada_res[1] <- length(intersect(indexbhbcadaf, index_alter)) / length(index_alter)
      fastbhbcada_res[2] <- length(setdiff(indexbhbcadaf, index_alter)) / length(indexbhbcadaf)
    }
    ## Storey
    storey_res <- storey_fun_cpp(p_all, alpha = alpha)
    index_storey <- storey_res$index_select + 1
    storey_res <- numeric(2)
    if (length(index_storey) > 0) {
      storey_res[1] <- length(intersect(index_storey, index_alter)) / length(index_alter)
      storey_res[2] <- length(setdiff(index_storey, index_alter)) / length(index_storey)
    }
    c(bh_res, bc_res, bhbc_res, bhbcada_res, fastbhbcada_res, storey_res)
  }
  res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
  colnames(res_mat) <- c("BH.power", "BH.FDP", "BC.power", "BC.FDP",
                         "eBH_Ave.power", "eBH_Ave.FDP", "eBH_Ada.power", "eBH_Ada.FDP",
                         "fast_eBH_Ada.power",
                         "fast_eBH_Ada.FDP", "ST.power", "ST.FDP")
  saveRDS(res_mat, paste0(
    "ApplicationII/results/bcbh_ratio", ratio, "sign", signa_stre, ".rds"
  ))
  print(round(colMeans(res_mat), 4))
}
stopCluster(cl)
registerDoSEQ()
