set.seed(2)
n_simu <- 1000
alpha <- 0.05
alpha_ebh <- 0.05
alpha_bh <- alpha_ebh
alpha_abh <- alpha_ebh
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
library(Rcpp)
library(RcppArmadillo)
sourceCpp("ApplicationI/utilities.cpp")
## generate p-value
p1 <- c(runif(n_vec[1] - n1_vec[1]), rbeta(n1_vec[1], shape1 = 4, shape2 = 500))
p2 <- c(runif(n_vec[2] - n1_vec[2]), rbeta(n1_vec[2], shape1 = 0.1, shape2 = 500))
pvalue_list <- list(p1, p2)

## initialization
n_group <- length(pvalue_list)
n_vec <- numeric(n_group)
## BC method
index_select_list <- list()
evec <- NULL
for (iter_group in seq_len(n_group)) {
  pvalue_tmp <- pvalue_list[[iter_group]]
  index_alter_tmp <- index_alter_list[[iter_group]]
  res_tmp <- bc_fun_cpp(pvalue_tmp, alpha)
  evec <- c(evec, res_tmp$evalue)
}
index_select <- which(evec != 0)
n_select <- length(index_select)
n1 <- sum(index_select <= 100)
n2 <- n_select - n1
evec <- evec[index_select]
## weighted version
res <- ebh_apa_cpp(pvalue_list, alpha_abh, alpha_ebh)
evec_weight <- res$evalue[index_select]

## plot
plot_df <- data.frame("evalue" = c(evec, evec_weight),
                      "index" = rep(seq_len(n_select), 2),
                      "weight" = c(rep("without", n_select), rep("with", n_select)),
                      "group" = c(rep(c(rep("group1", n1), rep("group2", n2)), 2)))

library(ggplot2)
p1 <- ggplot(plot_df, aes(x = index, y = evalue, color = group)) +
  geom_point(aes(shape = weight))
p1