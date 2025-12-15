library(Rcpp)
library(RcppArmadillo)
sourceCpp("ApplicationIII/utilities.cpp")
#### get initial value
solvek_fdr <- function(pi0, pvals, pvals_cutoff = 1e-15) {
  pvals[pvals <= pvals_cutoff] <- pvals_cutoff
  pvals[pvals >= 1 - pvals_cutoff] <- 1 - pvals_cutoff

  pi1 <- 1 - pi0

  fun <- function(k) {
    pk <- pvals^(-k)
    sum(-pi1 * pk * (1 + (1 - k) * log(pvals)) / (pi0 + pi1 * (1 - k) * pk))
  }

  err <- try(
    {
      k <- uniroot(fun, c(0.01, 0.99))$root
    },
    silent = TRUE
  )

  if (class(err) == "try-error") {
    k <- 0.75
  }
  return(k)
}

## main function
ebh_gbc <- function(pvalue_vec, pi0_var = NULL, f1_var = NULL, n_try = 10, n_split = 10,
                    alpha_ebh = 0.05, alpha_bc = 0.025, pvals_cutoff = 1e-15,
                    EM_paras = list(iterlim = 50, tol = 1e-5, k_init = NULL, pi0_init = NULL, nlm_iter = 5),
                    weight_bool = TRUE, parallel = FALSE, n_cores = 1, test = FALSE) {
  #### initialization
  n <- length(pvalue_vec)
  ## EM parameters
  if (is.null(EM_paras$pi0_init)) {
    pi0_est <- qvalue::qvalue(pvalue_vec, pi0.method = "bootstrap")$pi0
    if (pi0_est >= 0.99) pi0_est <- 0.99
    EM_paras$pi0_init <- pi0_est
  } else {
    pi0_est <- EM_paras$pi0_init
  }
  if (is.null(EM_paras$k_init)) {
    k_init <- solvek_fdr(pi0_est, pvalue_vec, pvals_cutoff)
    EM_paras$k_init <- k_init
  } else {
    k_init <- EM_paras$k_init
  }
  if (is.null(EM_paras$iterlim)) {
    iterlim <- 50
    EM_paras$iterlim <- iterlim
  } else {
    iterlim <- EM_paras$iterlim
  }
  if (is.null(EM_paras$nlm_iter)) {
    nlm_iter <- 5
    EM_paras$nlm_iter <- nlm_iter
  } else {
    nlm_iter <- EM_paras$nlm_iter
  }
  if (is.null(EM_paras$tol)) {
    tol <- 1e-5
    EM_paras$tol <- tol
  } else {
    tol <- EM_paras$tol
  }
  ## covariate information
  if (is.null(pi0_var)) {
    pi0_var <- matrix(rep(1, n))
  } else {
    pi0_var <- cbind(rep(1, n), pi0_var)
  }
  if (is.null(f1_var)) {
    f1_var <- matrix(rep(1, n))
  } else {
    f1_var <- cbind(rep(1, n), f1_var)
  }
  ## test covariate information
  if (test) {
    test_res <- suppressWarnings(cov_test_n(pvalue_vec, covariate = pi0_var[, 2]))
  } else {
    test_res <- list(p.value = 0) 
  }
  ## gbh-ebh
  if (test_res$p.value < 5e-3) {
    ## gbh-ebh
    if (parallel) {
      library(parallel)
      out_res <- mclapply(seq_len(n_try), function(iter_try) {
        split_groups <- shuffle_fun(n, n_split)
        evalue_tmp <- crossfit_fun(
          pvalue_vec, pi0_var, f1_var, alpha_bc, pi0_est,
          k_init, tol, iterlim, nlm_iter, n_split,
          split_groups, weight_bool
        )
        return(evalue_tmp)
      }, mc.cores = min(n_cores, n_try))
      evalue_mat <- matrix(unlist(out_res), nrow = n_try, byrow = TRUE)
      evalue_ave <- colMeans(evalue_mat)
      index_select <- ebh_fun_cpp(evalue_ave, alpha_ebh)
      outres <- list(
        index_select = index_select,
        evalue_mat = evalue_mat,
        evalue_ave = evalue_ave
      )
      outres$index_select <- outres$index_select + 1
    } else {
      outres <- suppressWarnings(ebh_gbc_cpp(
        pvalue_vec, pi0_var, f1_var, n_try, n_split,
        alpha_ebh, alpha_bc, pi0_est, k_init, tol, iterlim, nlm_iter,
        weight_bool
      ))
      outres$index_select <- outres$index_select + 1
    }
  } else {
    cat("There is no information in the covariate, the BC method is implement.")
    outres <- bc_fun_cpp(pvalue_vec, alpha_ebh)
    outres$index_select <- outres$index_select + 1
  }
  return(outres)
}

############################################################################################################
# Testing the informativeness of the covariate
cochran_armitage_test <- function(x, alternative = c("two.sided", "increasing", "decreasing")) {
  DNAME <- deparse(substitute(x))
  if (!(any(dim(x) == 2))) {
    stop("Cochran-Armitage test for trend must be used with rx2-table",
      call. = FALSE
    )
  }
  if (dim(x)[2] != 2) {
    x <- t(x)
  }
  nidot <- apply(x, 1, sum)
  n <- sum(nidot)
  Ri <- 1:dim(x)[1]
  Rbar <- sum(nidot * Ri) / n
  s2 <- sum(nidot * (Ri - Rbar)^2)
  pdot1 <- sum(x[, 1]) / n
  z <- sum(x[, 1] * (Ri - Rbar)) / sqrt(pdot1 * (1 - pdot1) *
    s2)
  STATISTIC <- z
  alternative <- match.arg(alternative)
  PVAL <- switch(alternative,
    two.sided = 2 * pnorm(abs(z),
      lower.tail = FALSE
    ),
    increasing = pnorm(z),
    decreasing = pnorm(z,
      lower.tail = FALSE
    )
  )
  PARAMETER <- dim(x)[1]
  names(STATISTIC) <- "Z"
  names(PARAMETER) <- "dim"
  METHOD <- "Cochran-Armitage test for trend"
  structure(list(
    statistic = STATISTIC, parameter = PARAMETER,
    alternative = alternative, p.value = PVAL, method = METHOD,
    data.name = DNAME
  ), class = "htest")
}

cov_test_n <- function(pvals, covariate, cutoffs = quantile(pvals, c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)),
                       grps = c(2, 4, 8, 16, 32), perm.no = 999, n.max = 100000, silence = TRUE) {
  save.seed <- get(".Random.seed", .GlobalEnv)

  n <- length(pvals)
  if (n > n.max) {
    ind <- sample(1:n, n.max)
    pvals <- pvals[ind]
    covariate <- covariate[ind]
  }


  stat.p.a <- stat.p.a1 <- stat.p.a2 <- array(NA, c(length(cutoffs), length(grps), perm.no))
  stat.o.a <- stat.o.a1 <- stat.o.a2 <- array(NA, c(length(cutoffs), length(grps)))

  for (i in 1:length(cutoffs)) {
    cutoff <- cutoffs[i]
    x <- as.numeric(pvals <= cutoff)
    for (j in 1:length(grps)) {
      if (!silence) cat(".")
      grp <- grps[j]
      y <- cut(covariate, c(min(covariate) - 0.01, quantile(covariate, 1 / grp * (1:(grp - 1))), max(covariate) + 0.01))

      if (!silence) cat(".")
      mat <- table(x, y)
      test.obj1 <- chisq.test(mat)
      test.obj2 <- cochran_armitage_test(mat)
      stat.o.a1[i, j] <- -pchisq(test.obj1$statistic, df = test.obj1$parameter, lower.tail = FALSE, log.p = TRUE)
      stat.o.a2[i, j] <- -log(2) - pnorm(abs(test.obj2$statistic), lower.tail = FALSE, log.p = TRUE)
      stat.o.a[i, j] <- max(c(stat.o.a1[i, j], stat.o.a2[i, j]))

      assign(".Random.seed", save.seed, .GlobalEnv)

      lpvs <- sapply(1:perm.no, function(k) {
        xp <- sample(x)
        mat <- table(xp, y)
        test.obj1 <- chisq.test(mat)
        test.obj2 <- cochran_armitage_test(mat)
        c(
          -pchisq(test.obj1$statistic, df = test.obj1$parameter, lower.tail = FALSE, log.p = TRUE),
          -log(2) - pnorm(abs(test.obj2$statistic), lower.tail = FALSE, log.p = TRUE)
        )
      })

      stat.p.a1[i, j, ] <- lpvs[1, ]
      stat.p.a2[i, j, ] <- lpvs[2, ]
      stat.p.a[i, j, ] <- pmax(lpvs[1, ], lpvs[2, ])
    }
  }

  stat.o1 <- matrixStats::colMaxs(stat.o.a)
  stat.o2 <- matrixStats::rowMaxs(stat.o.a)
  stat.o <- max(stat.o1)
  x.cut.optim <- grps[which.max(stat.o1)]
  p.cut.optim <- cutoffs[which.max(stat.o2)]
  stat.p <- apply(stat.p.a, 3, max)
  p.value <- mean(c(stat.p, stat.o) >= stat.o)

  stat.o1 <- max(stat.o.a1)
  stat.p1 <- apply(stat.p.a1, 3, max)
  p.value1 <- mean(c(stat.p1, stat.o1) >= stat.o1)

  stat.o2 <- max(stat.o.a2)
  stat.p2 <- apply(stat.p.a2, 3, max)
  p.value2 <- mean(c(stat.p2, stat.o2) >= stat.o2)

  return(list(
    stat.o = stat.o, stat.p = stat.p, p.value = p.value,
    x.cut.optim = x.cut.optim, p.cut.optim = p.cut.optim
  ))
}
