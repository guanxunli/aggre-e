## define parameters
library(Rcpp)
library(RcppArmadillo)
sourceCpp("ApplicationII/utilities.cpp")
n_simu <- 500
alpha <- alpha_ebh <- 0.05
alpha_bh <- alpha_ebh / 2
alpha_bha <- alpha_ebh / (1 + alpha_ebh)
n <- 1e3
signa_stre_vec <- round(seq(0.3, 0.5, length.out = 7), 2)
#### BH > BC setting
ratio <- 5e-2
n1 <- n * ratio
index_alter <- seq_len(n1)
signa_stre <- 0.4
mu <- signa_stre * log(n)
set.seed(1)
time_mat <- matrix(NA, nrow = n_simu, ncol = 6)
for (iter_simu in seq_len(n_simu)) {
  if (iter_simu %% 10 == 0) print(iter_simu)
  x_all <- c(rnorm(n1, mean = mu, sd = 1), rnorm(n - n1))
  p_all <- 1 - pnorm(x_all)
  ## BH
  time1 <- Sys.time()
  bh_res <- bh_fun_cpp(p_all, alpha = alpha)
  indexbh <- bh_res$index_select + 1
  bh_res <- numeric(2)
  if (length(indexbh) > 0) {
    bh_res[1] <- length(intersect(indexbh, index_alter)) / length(index_alter)
    bh_res[2] <- length(setdiff(indexbh, index_alter)) / length(indexbh)
  }
  time_mat[iter_simu, 1] <- as.numeric(Sys.time() - time1, units = "secs")
  ## BC
  time1 <- Sys.time()
  bc_res <- bc_fun_cpp(p_all, alpha = alpha)
  indexbc <- bc_res$index_select + 1
  bc_res <- numeric(2)
  if (length(indexbc) > 0) {
    bc_res[1] <- length(intersect(indexbc, index_alter)) / length(index_alter)
    bc_res[2] <- length(setdiff(indexbc, index_alter)) / length(indexbc)
  }
  time_mat[iter_simu, 2] <- as.numeric(Sys.time() - time1, units = "secs")
  ## ave
  time1 <- Sys.time()
  bhbc_res <- bhbc_fix(p_all, alpha_bh = alpha_bh, alpha_ebh = alpha_ebh)
  indexbhbc <- bhbc_res$index_select + 1
  bhbc_res <- numeric(2)
  if (length(indexbhbc) > 0) {
    bhbc_res[1] <- length(intersect(indexbhbc, index_alter)) / length(index_alter)
    bhbc_res[2] <- length(setdiff(indexbhbc, index_alter)) / length(indexbhbc)
  }
  time_mat[iter_simu, 3] <- as.numeric(Sys.time() - time1, units = "secs")
  ## ada
  time1 <- Sys.time()
  bhbcada_res <- bhbc_ada(p_all, alpha_bh = alpha_bha, alpha_ebh = alpha_ebh)
  indexbhbcada <- bhbcada_res$index_select + 1
  bhbcada_res <- numeric(2)
  if (length(indexbhbcada) > 0) {
    bhbcada_res[1] <- length(intersect(indexbhbcada, index_alter)) / length(index_alter)
    bhbcada_res[2] <- length(setdiff(indexbhbcada, index_alter)) / length(indexbhbcada)
  }
  time_mat[iter_simu, 4] <- as.numeric(Sys.time() - time1, units = "secs")
  ## fast ada
  time1 <- Sys.time()
  fastbhbcada_res <- fastbhbc_ada(p_all, alpha_bh = alpha_bha, alpha_ebh = alpha_ebh)
  indexbhbcadaf <- fastbhbcada_res$index_select + 1
  fastbhbcada_res <- numeric(2)
  if (length(indexbhbcada) > 0) {
    fastbhbcada_res[1] <- length(intersect(indexbhbcadaf, index_alter)) / length(index_alter)
    fastbhbcada_res[2] <- length(setdiff(indexbhbcadaf, index_alter)) / length(indexbhbcadaf)
  }
  time_mat[iter_simu, 5] <- as.numeric(Sys.time() - time1, units = "secs")
  ## Storey
  time1 <- Sys.time()
  storey_res <- storey_fun_cpp(p_all, alpha = alpha)
  index_storey <- storey_res$index_select + 1
  storey_res <- numeric(2)
  if (length(index_storey) > 0) {
    storey_res[1] <- length(intersect(index_storey, index_alter)) / length(index_alter)
    storey_res[2] <- length(setdiff(index_storey, index_alter)) / length(index_storey)
  }
  time_mat[iter_simu, 6] <- as.numeric(Sys.time() - time1, units = "secs")
}
print(round(colMeans(time_mat), 5))
