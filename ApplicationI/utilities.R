library(Rcpp)
library(RcppArmadillo)
sourceCpp("ApplicationI/utilities.cpp")

############################ joint BH ##################################
bc_joint <- function(pvalue_list, index_alter_list, alpha = 0.05) {
  ## initialization
  n_group <- length(pvalue_list)
  n_vec <- numeric(n_group)
  index_alter <- NULL
  pvalue_all <- NULL
  for (iter_group in seq_len(n_group)) {
    index_alter <- c(index_alter, index_alter_list[[iter_group]] + sum(n_vec))
    pvalue_tmp <- pvalue_list[[iter_group]]
    pvalue_all <- c(pvalue_all, pvalue_tmp)
    n_vec[iter_group] <- length(pvalue_tmp)
  }
  ## BC
  res <- bc_fun_cpp(pvalue_all, alpha)
  index_select <- res$index_select + 1
  ## save results
  res_vec <- numeric(2 + 2 * n_group)
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
    index_begin <- 0
    for (iter_group in seq_len(n_group)) {
      index_end <- index_begin + n_vec[iter_group]
      index_select_tmp <- index_select[intersect(
        which(index_select > index_begin),
        which(index_select < index_end + 1)
      )] -
        index_begin
      if (length(index_select_tmp) > 0) {
        index_alter_tmp <- index_alter_list[[iter_group]]
        pow <- length(intersect(index_select_tmp, index_alter_tmp)) / length(index_alter_tmp)
        fdp <- length(setdiff(index_select_tmp, index_alter_tmp)) / length(index_select_tmp)
        res_vec[c(2 * iter_group + 1, 2 * iter_group + 2)] <- c(pow, fdp)
      }
      index_begin <- index_end
    }
  }
  ## return
  return(res_vec)
}

############################ separate BH ############################
bc_separate <- function(pvalue_list, index_alter_list, alpha = 0.05) {
  ## initialization
  n_group <- length(pvalue_list)
  n_vec <- numeric(n_group)
  ## BC method
  index_select_list <- list()
  index_alter <- NULL
  index_select <- NULL
  res_vec <- numeric(2 + 2 * n_group)
  for (iter_group in seq_len(n_group)) {
    pvalue_tmp <- pvalue_list[[iter_group]]
    index_alter_tmp <- index_alter_list[[iter_group]]
    res_tmp <- bc_fun_cpp(pvalue_tmp, alpha)
    index_select_tmp <- res_tmp$index_select + 1
    if (length(index_select_tmp) > 0) {
      pow <- length(intersect(index_select_tmp, index_alter_tmp)) / length(index_alter_tmp)
      fdp <- length(setdiff(index_select_tmp, index_alter_tmp)) / length(index_select_tmp)
      res_vec[c(2 * iter_group + 1, 2 * iter_group + 2)] <- c(pow, fdp)
    }
    index_alter <- c(index_alter, index_alter_tmp + sum(n_vec))
    index_select <- c(index_select, index_select_tmp + sum(n_vec))
    n_vec[iter_group] <- length(pvalue_tmp)
  }
  ## save results
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
  }
  ## return
  return(res_vec)
}

############################ evalue method #############################
ebh_fun <- function(pvalue_list, index_alter_list, weight_list,
                    alpha_bc = 0.05, alpha_ebh = 0.05) {
  ## initialization
  n_group <- length(pvalue_list)
  n_vec <- numeric(n_group)
  ## BC method
  index_alter <- NULL
  evalue <- NULL
  for (iter_group in seq_len(n_group)) {
    pvalue_tmp <- pvalue_list[[iter_group]]
    index_alter_tmp <- index_alter_list[[iter_group]]
    index_alter <- c(index_alter, index_alter_tmp + sum(n_vec))
    res_tmp <- bc_fun_cpp(pvalue_tmp, alpha_bc)
    evalue_tmp <- res_tmp$evalue
    weight <- weight_list[[iter_group]]
    evalue <- c(evalue, weight * evalue_tmp)
    n_vec[iter_group] <- length(pvalue_tmp)
  }
  ## eBH
  index_select <- ebh_fun_cpp(evalue, alpha_ebh) + 1
  ## save results
  res_vec <- numeric(2 + 2 * n_group)
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
    index_begin <- 0
    for (iter_group in seq_len(n_group)) {
      index_end <- index_begin + n_vec[iter_group]
      index_select_tmp <- index_select[intersect(
        which(index_select > index_begin),
        which(index_select < index_end + 1)
      )] -
        index_begin
      if (length(index_select_tmp) > 0) {
        index_alter_tmp <- index_alter_list[[iter_group]]
        pow <- length(intersect(index_select_tmp, index_alter_tmp)) / length(index_alter_tmp)
        fdp <- length(setdiff(index_select_tmp, index_alter_tmp)) / length(index_select_tmp)
        res_vec[c(2 * iter_group + 1, 2 * iter_group + 2)] <- c(pow, fdp)
      }
      index_begin <- index_end
    }
  }
  ## return
  return(res_vec)
}

############################ adaptive ebh ############################
ebh_apa <- function(pvalue_list, index_alter_list,
                    alpha_bc = 0.05, alpha_ebh = 0.05) {
  ## initialization
  n_group <- length(pvalue_list)
  n_vec <- numeric(n_group)
  index_alter <- NULL
  for (iter_group in seq_len(n_group)) {
    index_alter <- c(index_alter, index_alter_list[[iter_group]] + sum(n_vec))
    pvalue_tmp <- pvalue_list[[iter_group]]
    n_vec[iter_group] <- length(pvalue_tmp)
  }
  ## BC
  res <- ebh_apa_cpp(pvalue_list, alpha_bc, alpha_ebh)
  index_select <- res$index_select + 1
  ## save results
  res_vec <- numeric(2 + 2 * n_group)
  # overall
  if (length(index_select) > 0) {
    all_pow <- length(intersect(index_select, index_alter)) / length(index_alter)
    all_fdp <- length(setdiff(index_select, index_alter)) / length(index_select)
    res_vec[c(1, 2)] <- c(all_pow, all_fdp)
    index_begin <- 0
    for (iter_group in seq_len(n_group)) {
      index_end <- index_begin + n_vec[iter_group]
      index_select_tmp <- index_select[intersect(
        which(index_select > index_begin),
        which(index_select < index_end + 1)
      )] -
        index_begin
      if (length(index_select_tmp) > 0) {
        index_alter_tmp <- index_alter_list[[iter_group]]
        pow <- length(intersect(index_select_tmp, index_alter_tmp)) / length(index_alter_tmp)
        fdp <- length(setdiff(index_select_tmp, index_alter_tmp)) / length(index_select_tmp)
        res_vec[c(2 * iter_group + 1, 2 * iter_group + 2)] <- c(pow, fdp)
      }
      index_begin <- index_end
    }
  }
  ## return
  return(res_vec)
}