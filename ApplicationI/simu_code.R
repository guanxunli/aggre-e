library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
n_simu <- 1000
alpha <- 0.05
alpha_ebh <- 0.05
alpha_bh <- alpha_ebh
alpha_abh <- alpha_ebh

######################### setting E1 ###########################
n_vec <- c(1e2, 1e3)
n <- sum(n_vec)
n1_vec <- c(20, 20)
index_alter_list <- list()
for (iter in seq_len(2)) {
  n_use <- n_vec[iter]
  n1_use <- n1_vec[iter]
  index_alter_list[[iter]] <- seq(n_use - n1_use + 1, n_use)
}
## results
set.seed(1)
cl <- makeCluster(10)
registerDoParallel(cl)
out_res <- foreach(iter_simu = seq_len(n_simu)) %dorng% {
  library(Rcpp)
  library(RcppArmadillo)
  source("ApplicationI/utilities.R")
  ## generate p-value
  p1 <- c(runif(n_vec[1] - n1_vec[1]), rbeta(n1_vec[1], shape1 = 4, shape2 = 500))
  p2 <- c(runif(n_vec[2] - n1_vec[2]), rbeta(n1_vec[2], shape1 = 0.1, shape2 = 500))
  pvalue_list <- list(p1, p2)
  ## joint BC
  bcj_res <- bc_joint(pvalue_list, index_alter_list, alpha)
  ## separate BC
  bcs_res <- bc_separate(pvalue_list, index_alter_list, alpha)
  ## ebh weight1
  ebh1_res <- ebh_fun(
    pvalue_list, index_alter_list, list(1, 1),
    alpha_bh, alpha_ebh
  )
  ## ebh weight2
  ebh2_res <- ebh_fun(
    pvalue_list, index_alter_list,
    list(n / (2 * n_vec[1]), n / (2 * n_vec[2])),
    alpha_bh, alpha_ebh
  )
  ## ebh weight3
  ebh3_res <- ebh_apa(
    pvalue_list, index_alter_list,
    alpha_abh, alpha_ebh
  )
  round(c(bcj_res, bcs_res, ebh1_res, ebh2_res, ebh3_res), 4)
}
stopCluster(cl)
res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
saveRDS(res_mat, paste0(
  "ApplicationI/results/setting1.rds"
))
res <- round(colMeans(res_mat), 3)
#### print results
cat(
  "BC\\_j &", res[1], "&", res[3], "&", res[5], "&",
  res[2], "&", res[4], "&", res[6], "\\\\\n",
  "BC\\_s &", res[7], "&", res[9], "&", res[11], "&",
  res[8], "&", res[10], "&", res[12], "\\\\\n",
  "eBH\\_1 &", res[13], "&", res[15], "&", res[17], "&",
  res[14], "&", res[16], "&", res[18], "\\\\\n",
  "eBH\\_2 &", res[19], "&", res[21], "&", res[23], "&",
  res[20], "&", res[22], "&", res[24], "\\\\\n",
  "eBH\\_a &", res[25], "&", res[27], "&", res[29], "&",
  res[26], "&", res[28], "&", res[30], "\\\\\n"
)

######################### setting E2 ###########################
n_vec <- c(1e2, 1e3)
n <- sum(n_vec)
n1_vec <- c(20, 20)
index_alter_list <- list()
for (iter in seq_len(2)) {
  n_use <- n_vec[iter]
  n1_use <- n1_vec[iter]
  index_alter_list[[iter]] <- seq(n_use - n1_use + 1, n_use)
}
## results
set.seed(1)
cl <- makeCluster(10)
registerDoParallel(cl)
out_res <- foreach(iter_simu = seq_len(n_simu)) %dorng% {
  library(Rcpp)
  library(RcppArmadillo)
  source("ApplicationI/utilities.R")
  ## generate p-value
  p1 <- c(runif(n_vec[1] - n1_vec[1]), rbeta(n1_vec[1], shape1 = 0.5, shape2 = 500))
  p2 <- c(runif(n_vec[2] - n1_vec[2]), rbeta(n1_vec[2], shape1 = 0.5, shape2 = 500))
  pvalue_list <- list(p1, p2)
  ## joint BC
  bcj_res <- bc_joint(pvalue_list, index_alter_list, alpha)
  ## separate BC
  bcs_res <- bc_separate(pvalue_list, index_alter_list, alpha)
  ## ebh weight1
  ebh1_res <- ebh_fun(
    pvalue_list, index_alter_list, list(1, 1),
    alpha_bh, alpha_ebh
  )
  ## ebh weight2
  ebh2_res <- ebh_fun(
    pvalue_list, index_alter_list,
    list(n / (2 * n_vec[1]), n / (2 * n_vec[2])),
    alpha_bh, alpha_ebh
  )
  ## ebh weight3
  ebh3_res <- ebh_apa(
    pvalue_list, index_alter_list,
    alpha_abh, alpha_ebh
  )
  round(c(bcj_res, bcs_res, ebh1_res, ebh2_res, ebh3_res), 4)
}
stopCluster(cl)
res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
saveRDS(res_mat, paste0(
  "ApplicationI/results/setting2.rds"
))
res <- round(colMeans(res_mat), 3)
#### print results
cat(
  "BC\\_j &", res[1], "&", res[3], "&", res[5], "&",
  res[2], "&", res[4], "&", res[6], "\\\\\n",
  "BC\\_s &", res[7], "&", res[9], "&", res[11], "&",
  res[8], "&", res[10], "&", res[12], "\\\\\n",
  "eBH\\_1 &", res[13], "&", res[15], "&", res[17], "&",
  res[14], "&", res[16], "&", res[18], "\\\\\n",
  "eBH\\_2 &", res[19], "&", res[21], "&", res[23], "&",
  res[20], "&", res[22], "&", res[24], "\\\\\n",
  "eBH\\_a &", res[25], "&", res[27], "&", res[29], "&",
  res[26], "&", res[28], "&", res[30], "\\\\\n"
)

######################### setting F1 ###########################
n_vec <- c(1e2, 1e2, 1e3, 1e3)
n <- sum(n_vec)
n1_vec <- c(20, 20, 20, 20)
index_alter_list <- list()
for (iter in seq_len(4)) {
  n_use <- n_vec[iter]
  n1_use <- n1_vec[iter]
  index_alter_list[[iter]] <- seq(n_use - n1_use + 1, n_use)
}
## results
set.seed(1)
cl <- makeCluster(10)
registerDoParallel(cl)
out_res <- foreach(iter_simu = seq_len(n_simu)) %dorng% {
  library(Rcpp)
  library(RcppArmadillo)
  source("ApplicationI/utilities.R")
  ## generate p-value
  p1 <- c(runif(n_vec[1] - n1_vec[1]), rbeta(n1_vec[1], shape1 = 0.1, shape2 = 500))
  p2 <- c(runif(n_vec[2] - n1_vec[2]), rbeta(n1_vec[2], shape1 = 0.1, shape2 = 500))
  p3 <- c(runif(n_vec[3] - n1_vec[3]), rbeta(n1_vec[3], shape1 = 0.1, shape2 = 500))
  p4 <- c(runif(n_vec[4] - n1_vec[4]), rbeta(n1_vec[4], shape1 = 0.1, shape2 = 500))
  pvalue_list <- list(p1, p2, p3, p4)
  ## joint BH
  bcj_res <- bc_joint(pvalue_list, index_alter_list, alpha)
  ## separate BH
  bcs_res <- bc_separate(pvalue_list, index_alter_list, alpha)
  ## naive BH
  bcsn_res <- bc_separate(pvalue_list, index_alter_list, alpha / 2)
  ## ebh weight1
  ebh1_res <- ebh_fun(
    pvalue_list, index_alter_list, list(1, 1, 1, 1),
    alpha_bh, alpha_ebh
  )
  ## ebh weight2
  ebh2_res <- ebh_fun(
    pvalue_list, index_alter_list,
    list(
      n / (4 * n_vec[1]), n / (4 * n_vec[2]),
      n / (4 * n_vec[3]), n / (4 * n_vec[4])
    ),
    alpha_bh, alpha_ebh
  )
  ## ebh weight3
  ebh3_res <- ebh_apa(
    pvalue_list, index_alter_list,
    alpha_abh, alpha_ebh
  )
  round(c(bcj_res, bcs_res, bcsn_res, ebh1_res, ebh2_res, ebh3_res), 4)
}
stopCluster(cl)
res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
saveRDS(res_mat, paste0(
  "ApplicationI/results/setting3.rds"
))
res <- round(colMeans(res_mat), 3)
#### print results
cat(
  "BC\\_j &", res[1], "&", res[3], "&", res[5], "&", res[7], "&", res[9], "&",
  res[2], "&", res[4], "&", res[6], "&", res[8], "&", res[10], "\\\\\n",
  "BC\\_s &", res[11], "&", res[13], "&", res[15], "&", res[17], "&", res[19], "&",
  res[12], "&", res[14], "&", res[16], "&", res[18], "&", res[20], "\\\\\n",
  "BC\\_sn &", res[21], "&", res[23], "&", res[25], "&", res[27], "&", res[29], "&",
  res[22], "&", res[24], "&", res[26], "&", res[28], "&", res[30], "\\\\\n",
  "eBH\\_1 &", res[31], "&", res[33], "&", res[35], "&", res[37], "&", res[39], "&",
  res[32], "&", res[34], "&", res[36], "&", res[38], "&", res[40], "\\\\\n",
  "eBH\\_2 &", res[41], "&", res[43], "&", res[45], "&", res[47], "&", res[49], "&",
  res[42], "&", res[44], "&", res[46], "&", res[48], "&", res[50], "\\\\\n",
  "eBH\\_a &", res[51], "&", res[53], "&", res[55], "&", res[57], "&", res[59], "&",
  res[52], "&", res[54], "&", res[56], "&", res[58], "&", res[60], "\\\\\n"
)

######################### setting F2 ###########################
## define parameters
n_vec <- c(1e2, 1e2, 1e2, 1e2)
n <- sum(n_vec)
n1_vec <- c(1, 20, 20, 20)
index_alter_list <- list()
for (iter in seq_len(4)) {
  n_use <- n_vec[iter]
  n1_use <- n1_vec[iter]
  index_alter_list[[iter]] <- seq(n_use - n1_use + 1, n_use)
}
## results
set.seed(1)
cl <- makeCluster(10)
registerDoParallel(cl)
out_res <- foreach(iter_simu = seq_len(n_simu)) %dorng% {
  library(Rcpp)
  library(RcppArmadillo)
  source("ApplicationI/utilities.R")
  ## generate p-value
  p1 <- c(runif(n_vec[1] - n1_vec[1]), rbeta(n1_vec[1], shape1 = 0.01, shape2 = 5000))
  p2 <- c(runif(n_vec[2] - n1_vec[2]), rbeta(n1_vec[2], shape1 = 0.1, shape2 = 500))
  p3 <- c(runif(n_vec[3] - n1_vec[3]), rbeta(n1_vec[3], shape1 = 0.2, shape2 = 500))
  p4 <- c(runif(n_vec[4] - n1_vec[4]), rbeta(n1_vec[4], shape1 = 0.3, shape2 = 500))
  pvalue_list <- list(p1, p2, p3, p4)
  ## joint BH
  bcj_res <- bc_joint(pvalue_list, index_alter_list, alpha)
  ## separate BH
  bcs_res <- bc_separate(pvalue_list, index_alter_list, alpha)
  ## naive BH
  bcsn_res <- bc_separate(pvalue_list, index_alter_list, alpha / 2)
  ## ebh weight1
  ebh1_res <- ebh_fun(
    pvalue_list, index_alter_list, list(1, 1, 1, 1),
    alpha_bh, alpha_ebh
  )
  ## ebh weight2
  ebh2_res <- ebh_fun(
    pvalue_list, index_alter_list,
    list(
      n / (4 * n_vec[1]), n / (4 * n_vec[2]),
      n / (4 * n_vec[3]), n / (4 * n_vec[4])
    ),
    alpha_bh, alpha_ebh
  )
  ## ebh weight3
  ebh3_res <- ebh_apa(
    pvalue_list, index_alter_list,
    alpha_abh, alpha_ebh
  )
  round(c(bcj_res, bcs_res, bcsn_res, ebh1_res, ebh2_res, ebh3_res), 4)
}
stopCluster(cl)
res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
saveRDS(res_mat, paste0(
  "ApplicationI/results/setting4.rds"
))
res <- round(colMeans(res_mat), 3)
#### print results
cat(
  "BC\\_j &", res[1], "&", res[3], "&", res[5], "&", res[7], "&", res[9], "&",
  res[2], "&", res[4], "&", res[6], "&", res[8], "&", res[10], "\\\\\n",
  "BC\\_s &", res[11], "&", res[13], "&", res[15], "&", res[17], "&", res[19], "&",
  res[12], "&", res[14], "&", res[16], "&", res[18], "&", res[20], "\\\\\n",
  "BC\\_sn &", res[21], "&", res[23], "&", res[25], "&", res[27], "&", res[29], "&",
  res[22], "&", res[24], "&", res[26], "&", res[28], "&", res[30], "\\\\\n",
  "eBH\\_1 &", res[31], "&", res[33], "&", res[35], "&", res[37], "&", res[39], "&",
  res[32], "&", res[34], "&", res[36], "&", res[38], "&", res[40], "\\\\\n",
  "eBH\\_2 &", res[41], "&", res[43], "&", res[45], "&", res[47], "&", res[49], "&",
  res[42], "&", res[44], "&", res[46], "&", res[48], "&", res[50], "\\\\\n",
  "eBH\\_a &", res[51], "&", res[53], "&", res[55], "&", res[57], "&", res[59], "&",
  res[52], "&", res[54], "&", res[56], "&", res[58], "&", res[60], "\\\\\n"
)

######################### setting F3 ###########################
alpha <- 0.2
alpha_ebh <- 0.2
alpha_bh <- alpha_ebh
alpha_abh <- alpha_ebh
## define parameters
n_vec <- c(50, 100, 50, 100)
n <- sum(n_vec)
n1_vec <- c(2, 2, 4, 4)
index_alter_list <- list()
for (iter in seq_len(4)) {
  n_use <- n_vec[iter]
  n1_use <- n1_vec[iter]
  index_alter_list[[iter]] <- seq(n_use - n1_use + 1, n_use)
}
## results
set.seed(1)
cl <- makeCluster(10)
registerDoParallel(cl)
out_res <- foreach(iter_simu = seq_len(n_simu)) %dorng% {
  library(Rcpp)
  library(RcppArmadillo)
  source("ApplicationI/utilities.R")
  ## generate p-value
  p1 <- c(runif(n_vec[1] - n1_vec[1]), rbeta(n1_vec[1], shape1 = 0.1, shape2 = 500))
  p2 <- c(runif(n_vec[2] - n1_vec[2]), rbeta(n1_vec[2], shape1 = 0.1, shape2 = 500))
  p3 <- c(runif(n_vec[3] - n1_vec[3]), rbeta(n1_vec[3], shape1 = 0.2, shape2 = 500))
  p4 <- c(runif(n_vec[4] - n1_vec[4]), rbeta(n1_vec[4], shape1 = 0.3, shape2 = 500))
  pvalue_list <- list(p1, p2, p3, p4)
  ## joint BH
  bcj_res <- bc_joint(pvalue_list, index_alter_list, alpha)
  ## separate BH
  bcs_res <- bc_separate(pvalue_list, index_alter_list, alpha)
  ## naive BH
  bcsn_res <- bc_separate(pvalue_list, index_alter_list, alpha / 2)
  ## ebh weight1
  ebh1_res <- ebh_fun(
    pvalue_list, index_alter_list, list(1, 1, 1, 1),
    alpha_bh, alpha_ebh
  )
  ## ebh weight2
  ebh2_res <- ebh_fun(
    pvalue_list, index_alter_list,
    list(
      n / (4 * n_vec[1]), n / (4 * n_vec[2]),
      n / (4 * n_vec[3]), n / (4 * n_vec[4])
    ),
    alpha_bh, alpha_ebh
  )
  ## ebh weight3
  ebh3_res <- ebh_apa(
    pvalue_list, index_alter_list,
    alpha_abh, alpha_ebh
  )
  round(c(bcj_res, bcs_res, bcsn_res, ebh1_res, ebh2_res, ebh3_res), 4)
}
stopCluster(cl)
res_mat <- matrix(unlist(out_res), nrow = n_simu, byrow = TRUE)
saveRDS(res_mat, paste0(
  "ApplicationI/results/setting5.rds"
))
res <- round(colMeans(res_mat), 3)
#### print results
cat(
  "BC\\_j &", res[1], "&", res[3], "&", res[5], "&", res[7], "&", res[9], "&",
  res[2], "&", res[4], "&", res[6], "&", res[8], "&", res[10], "\\\\\n",
  "BC\\_s &", res[11], "&", res[13], "&", res[15], "&", res[17], "&", res[19], "&",
  res[12], "&", res[14], "&", res[16], "&", res[18], "&", res[20], "\\\\\n",
  "BC\\_sn &", res[21], "&", res[23], "&", res[25], "&", res[27], "&", res[29], "&",
  res[22], "&", res[24], "&", res[26], "&", res[28], "&", res[30], "\\\\\n",
  "eBH\\_1 &", res[31], "&", res[33], "&", res[35], "&", res[37], "&", res[39], "&",
  res[32], "&", res[34], "&", res[36], "&", res[38], "&", res[40], "\\\\\n",
  "eBH\\_2 &", res[41], "&", res[43], "&", res[45], "&", res[47], "&", res[49], "&",
  res[42], "&", res[44], "&", res[46], "&", res[48], "&", res[50], "\\\\\n",
  "eBH\\_a &", res[51], "&", res[53], "&", res[55], "&", res[57], "&", res[59], "&",
  res[52], "&", res[54], "&", res[56], "&", res[58], "&", res[60], "\\\\\n"
)
