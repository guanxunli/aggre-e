library(parallel)
library(doParallel)
library(doRNG)
## define parameters
n_simu <- 500
n_cores <- 50
alpha <- alpha_ebh <- 0.05
alpha_bh <- alpha_ebh / 2
alpha_bha <- alpha_ebh / (1 + alpha_ebh)
n <- 1e3
signa_stre_vec <- round(seq(0.3, 0.5, length.out = 7), 2)
ratio <- 5e-2
n1 <- n * ratio
index_alter <- seq_len(n1)
for (signa_stre in signa_stre_vec) {
  print(signa_stre)
  mu <- signa_stre * log(n)
  set.seed(2024)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  out_res <- foreach(iter_simu = seq_len(n_simu)) %dorng% {
    library(Rcpp)
    library(RcppArmadillo)
    sourceCpp("ApplicationII/utilities.cpp")
    x_all <- c(rnorm(n1, mean = mu, sd = 1), rnorm(n - n1))
    p_all <- 1 - pnorm(x_all)
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
    if (length(indexbhbcada) > 0) {
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
    round(c(bh_res, bc_res, bhbc_res, bhbcada_res, fastbhbcada_res, storey_res), 4)
  }
  stopCluster(cl)
  res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
  saveRDS(res_mat, paste0(
    "ApplicationII/results/bhbc_ratio", ratio, "sign", signa_stre, ".rds"
  ))
  print(round(colMeans(res_mat), 4))
}