## BH method
bh_fun <- function(p_all, index_alter, alpha = 0.05) {
  p_adj <- p.adjust(p_all, method = "BH")
  index_select <- which(p_adj <= alpha)
  ## return results
  res_vec <- numeric(2)
  if (length(index_select) > 0) {
    res_vec[1] <- length(intersect(index_select, index_alter)) / length(index_alter)
    res_vec[2] <- length(setdiff(index_select, index_alter)) / length(index_select)
  }
  return(res_vec)
}

## BC method
bc_fun <- function(p_all, index_alter, alpha = 0.05) {
  n <- length(p_all)
  ordered_p <- sort(p_all, decreasing = TRUE)
  for (iter in seq_len(n + 1)) {
    if (iter == n + 1) {
      index_select <- NULL
    } else {
      p_tmp <- ordered_p[iter]
      n_num <- 1 + sum(p_all >= 1 - p_tmp)
      n_dom <- max(1, sum(p_all <= p_tmp))
      hat_fdp <- n_num / n_dom
      if (hat_fdp <= alpha) {
        index_select <- which(p_all <= p_tmp)
        break
      }
    }
  }
  ## return results
  res_vec <- numeric(2)
  if (length(index_select) > 0) {
    res_vec[1] <- length(intersect(index_select, index_alter)) / length(index_alter)
    res_vec[2] <- length(setdiff(index_select, index_alter)) / length(index_select)
  }
  return(res_vec)
}

## BHBC method
bhbc_fun <- function(p_all, index_alter, alpha_bh = 0.05, alpha_ebh = 0.05) {
  n <- length(p_all)
  ## bh method
  e_tmp_bh <- 0
  evalue_bh <- numeric(n)
  p_adj <- p.adjust(p_all, method = "BH")
  index_select_bh <- which(p_adj <= alpha_bh)
  if (length(index_select_bh) > 0) {
    e_tmp_bh <- 1 / max(p_all[index_select_bh])
    evalue_bh[index_select_bh] <- e_tmp_bh
  }

  ## bc method
  e_tmp_bc <- 0
  evalue_bc <- numeric(n)
  ordered_p <- sort(p_all, decreasing = TRUE)
  for (iter in seq_len(n + 1)) {
    if (iter == n + 1) {
      index_select_bc <- NULL
    } else {
      p_tmp <- ordered_p[iter]
      n_num <- 1 + sum(p_all >= 1 - p_tmp)
      n_dom <- max(1, sum(p_all <= p_tmp))
      hat_fdp <- n_num / n_dom
      if (hat_fdp <= alpha_bh) {
        index_select_bc <- which(p_all <= p_tmp)
        e_tmp_bc <- n / (1 + sum(p_all >= 1 - p_tmp))
        evalue_bc[index_select_bc] <- e_tmp_bc
        break
      }
    }
  }

  ## combine e-value
  evalue <- (evalue_bc + evalue_bh) / 2

  ## e-BH
  orderd_evalue <- sort(evalue, decreasing = TRUE)
  for (iter in seq(n, 0, by = -1)) {
    if (iter == 0) {
      index_select <- NULL
      break
    } else {
      e_tmp <- orderd_evalue[iter]
      thres <- n / (alpha_ebh * iter)
      if (e_tmp >= thres) {
        index_select <- which(evalue >= thres)
        break
      }
    }
  }
  ## save results
  res_vec <- numeric(2)
  if (length(index_select) > 0) {
    res_vec[1] <- length(intersect(index_select, index_alter)) / length(index_alter)
    res_vec[2] <- length(setdiff(index_select, index_alter)) / length(index_select)
  }
  return(res_vec)
}

## BHBC cross method
bhbc_cross <- function(p_all, index_alter, n_try = 1, n_split = 2,
                       alpha_bh = 0.05, alpha_ebh = 0.05) {
  ## initialization
  n <- length(p_all)
  evalue_mat <- matrix(NA, nrow = n_try, ncol = n)
  for (iter_try in seq_len(n_try)) {
    shuffled_index <- sample(seq_len(n))
    split_groups <- split(shuffled_index, ceiling(seq_len(n) * n_split / n))
    ## calculate evalue
    ebh_group <- ebc_group <- numeric(n_split)
    evaluebh_group <- evaluebc_group <- list()
    for (iter_g in seq_len(n_split)) {
      index_use <- split_groups[[iter_g]]
      nuse <- length(index_use)
      puse <- p_all[index_use]
      ## bh method
      e_tmp_bh <- 0
      etmp_bh_vec <- numeric(nuse)
      p_adj <- p.adjust(puse, method = "BH")
      index_select_bh <- which(p_adj <= alpha_bh)
      if (length(index_select_bh) > 0) {
        e_tmp_bh <- 1 / max(puse[index_select_bh])
      }
      ebh_group[iter_g] <- e_tmp_bh
      etmp_bh_vec[index_select_bh] <- e_tmp_bh
      evaluebh_group[[iter_g]] <- etmp_bh_vec
      ## bc method
      e_tmp_bc <- 0
      etmp_bc_vec <- numeric(nuse)
      ordered_p <- sort(puse, decreasing = TRUE)
      for (iter_use in seq_len(nuse)) {
        p_tmp <- ordered_p[iter_use]
        n_num <- 1 + sum(puse >= 1 - p_tmp)
        n_dom <- max(1, sum(puse <= p_tmp))
        hat_fdp <- n_num / n_dom
        if (hat_fdp <= alpha_bh) {
          index_select_bc <- which(puse <= p_tmp)
          e_tmp_bc <- nuse / (1 + sum(puse >= 1 - p_tmp))
          etmp_bc_vec[index_select_bc] <- e_tmp_bc
          break
        }
      }
      ebc_group[iter_g] <- e_tmp_bc
      evaluebc_group[[iter_g]] <- etmp_bc_vec
    }
    ## calculate weight
    wbh_group <- wbc_group <- numeric(n_split)
    for (iter_g in seq_len(n_split)) {
      index_use <- unlist(split_groups[-iter_g])
      nuse <- length(index_use)
      puse <- p_all[index_use]
      ## bh method
      e_tmp_bh <- 0
      p_adj <- p.adjust(puse, method = "BH")
      index_select_bh <- which(p_adj <= alpha_bh)
      if (length(index_select_bh) > 0) {
        e_tmp_bh <- 1 / max(puse[index_select_bh])
      }
      ## bc method
      e_tmp_bc <- 0
      ordered_p <- sort(puse, decreasing = TRUE)
      for (iter_use in seq_len(nuse)) {
        p_tmp <- ordered_p[iter_use]
        n_num <- 1 + sum(puse >= 1 - p_tmp)
        n_dom <- max(1, sum(puse <= p_tmp))
        hat_fdp <- n_num / n_dom
        if (hat_fdp <= alpha_bh) {
          e_tmp_bc <- nuse / (1 + sum(p_all >= 1 - p_tmp))
          break
        }
      }
      ## calculate weight
      if (e_tmp_bh + e_tmp_bc == 0) {
        wbh_group[iter_g] <- wbc_group[iter_g] <- 0.5
      } else if (e_tmp_bh == 0) {
        wbh_group[iter_g] <- 0
        wbc_group[iter_g] <- 1
      } else if (e_tmp_bc == 0) {
        wbh_group[iter_g] <- 1
        wbc_group[iter_g] <- 0
      } else {
        wbh_group[iter_g] <- e_tmp_bh / (e_tmp_bh + e_tmp_bc)
        wbc_group[iter_g] <- e_tmp_bc / (e_tmp_bh + e_tmp_bc)
      }
    }

    ## combine e-value
    evalue <- numeric(n)
    for (iter_g in seq_len(n_split)) {
      index_use <- split_groups[[iter_g]]
      evalue[index_use] <- wbh_group[iter_g] * evaluebh_group[[iter_g]] +
        wbc_group[iter_g] * evaluebc_group[[iter_g]]
    }
    evalue_mat[iter_try, ] <- evalue
  }
  eavlue <- colMeans(evalue_mat)
  ## e-BH
  orderd_evalue <- sort(evalue, decreasing = TRUE)
  for (iter in seq(n, 0, by = -1)) {
    if (iter == 0) {
      index_select <- NULL
      break
    } else {
      e_tmp <- orderd_evalue[iter]
      thres <- n / (alpha_ebh * iter)
      if (e_tmp >= thres) {
        index_select <- which(evalue >= thres)
        break
      }
    }
  }
  ## save results
  res_vec <- numeric(2)
  if (length(index_select) > 0) {
    res_vec[1] <- length(intersect(index_select, index_alter)) / length(index_alter)
    res_vec[2] <- length(setdiff(index_select, index_alter)) / length(index_select)
  }
  return(res_vec)
}
