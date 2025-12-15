library(Rcpp)
library(RcppArmadillo)
sourceCpp("ApplicationI/utilities.cpp")

############################ joint BH ##################################
bc_joint <- function(p1, p2, index_alter1, index_alter2, alpha = 0.05) {
  ## initialization
  n1 <- length(p1)
  n2 <- length(p2)
  index_alter <- c(index_alter1, index_alter2 + n1)
  p_all <- c(p1, p2)
  res_vec <- numeric(6)
  ## BH
  res <- bc_fun_cpp(p_all, alpha)
  index_select <- res$index_select + 1
  index_select1 <- index_select[which(index_select < n1 + 1)]
  index_select2 <- index_select[which(index_select > n1)] - n1
  
  ## save results
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
  }
  # group1
  if (length(index_select1) > 0) {
    all_pow1 <- length(intersect(index_select1, index_alter1)) / length(index_alter1)
    all_fdp1 <- length(setdiff(index_select1, index_alter1)) / length(index_select1)
    res_vec[c(3, 4)] <- c(all_pow1, all_fdp1)
  }
  # group2
  if (length(index_select2) > 0) {
    all_pow2 <- length(intersect(index_select2, index_alter2)) / length(index_alter2)
    all_fdp2 <- length(setdiff(index_select2, index_alter2)) / length(index_select2)
    res_vec[c(5, 6)] <- c(all_pow2, all_fdp2)
  }
  ## return
  return(res_vec)
}

############################ assemble BH ############################
bc_assemble <- function(p1, p2, index_alter1, index_alter2, alpha = 0.025) {
  ## initialization
  n1 <- length(p1)
  n2 <- length(p2)
  index_alter <- c(index_alter1, index_alter2 + n1)
  res_vec <- numeric(6)
  ## BH
  # group 1
  res1 <- bc_fun_cpp(p1, alpha)
  index_select1 <- res1$index_select + 1
  # group 2
  res2 <- bc_fun_cpp(p2, alpha)
  index_select2 <- res2$index_select + 1
  # overall
  index_select <- c(index_select1, index_select2 + n1)
  
  ## save results
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
  }
  # group1
  if (length(index_select1) > 0) {
    all_pow1 <- length(intersect(index_select1, index_alter1)) / length(index_alter1)
    all_fdp1 <- length(setdiff(index_select1, index_alter1)) / length(index_select1)
    res_vec[c(3, 4)] <- c(all_pow1, all_fdp1)
  }
  # group2
  if (length(index_select2) > 0) {
    all_pow2 <- length(intersect(index_select2, index_alter2)) / length(index_alter2)
    all_fdp2 <- length(setdiff(index_select2, index_alter2)) / length(index_select2)
    res_vec[c(5, 6)] <- c(all_pow2, all_fdp2)
  }
  ## return
  return(res_vec)
}

############################ evalue method #############################
ebh_fun <- function(p1, p2, index_alter1, index_alter2, w1, w2,
                    alpha_bh = 0.05, alpha_ebh = 0.05) {
  ## initialization
  n1 <- length(p1)
  n2 <- length(p2)
  n <- n1 + n2
  index_alter <- c(index_alter1, index_alter2 + n1)
  res_vec <- numeric(6)
  ## BH
  # group 1
  res1 <- bc_fun_cpp(p1, alpha_bh)
  evalue1 <- w1 * res1$evalue
  # group 2
  res2 <- bc_fun_cpp(p2, alpha_bh)
  evalue2 <- w2 * res2$evalue
  ## eBH
  evalue <- c(evalue1, evalue2)
  index_select <- ebh_fun_cpp(evalue) + 1
  index_select1 <- index_select[which(index_select < n1 + 1)]
  index_select2 <- index_select[which(index_select > n1)] - n1
  
  ## save results
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
  }
  # group1
  if (length(index_select1) > 0) {
    all_pow1 <- length(intersect(index_select1, index_alter1)) / length(index_alter1)
    all_fdp1 <- length(setdiff(index_select1, index_alter1)) / length(index_select1)
    res_vec[c(3, 4)] <- c(all_pow1, all_fdp1)
  }
  # group2
  if (length(index_select2) > 0) {
    all_pow2 <- length(intersect(index_select2, index_alter2)) / length(index_alter2)
    all_fdp2 <- length(setdiff(index_select2, index_alter2)) / length(index_select2)
    res_vec[c(5, 6)] <- c(all_pow2, all_fdp2)
  }
  ## return
  return(res_vec)
}

############################ adaptive ebh ############################
ebh_apa <- function(p1, p2, index_alter1, index_alter2,
                    alpha_bh = 0.05, alpha_ebh = 0.05) {
  ## initialization
  n1 <- length(p1)
  res_vec <- numeric(6)
  index_alter <- c(index_alter1, index_alter2 + n1)
  ## ada ebh results
  ada_res <- ebh_apa_cpp(p1, p2, alpha_bh, alpha_ebh)
  index_select <- ada_res$index_select + 1
  index_select1 <- index_select[which(index_select < n1 + 1)]
  index_select2 <- index_select[which(index_select > n1)] - n1
  ## save results
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
  }
  # group1
  if (length(index_select1) > 0) {
    all_pow1 <- length(intersect(index_select1, index_alter1)) / length(index_alter1)
    all_fdp1 <- length(setdiff(index_select1, index_alter1)) / length(index_select1)
    res_vec[c(3, 4)] <- c(all_pow1, all_fdp1)
  }
  # group2
  if (length(index_select2) > 0) {
    all_pow2 <- length(intersect(index_select2, index_alter2)) / length(index_alter2)
    all_fdp2 <- length(setdiff(index_select2, index_alter2)) / length(index_select2)
    res_vec[c(5, 6)] <- c(all_pow2, all_fdp2)
  }
  ## return
  return(res_vec)
}
