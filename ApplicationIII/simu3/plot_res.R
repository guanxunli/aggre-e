set.seed(1)
library(ggplot2)
library(gridExtra)
n_simu <- 500
a0_vec <- c(3.5, 2.5, 1.5)
a1_vec <- c(1.5, 2, 2.5)
signa_sten_vec <- seq(2.5, 3.4, by = 0.15)
n_signa <- length(signa_sten_vec)
method_vec <- c(
  "BH", "IHW_storey", "IHW_betamix", "adaPT", "SABHA", "eBH_FBC"
)
n_method <- length(method_vec)

######################## plot one figures ########################
df_plot_all <- NULL
for (a0 in a0_vec) {
  for (a1 in a1_vec) {
    ## load results
    res_mat <- matrix(NA, nrow = n_signa, ncol = 2 * n_method)
    ci_mat_l <- matrix(NA, nrow = n_signa, ncol = 2 * n_method)
    ci_mat_r <- matrix(NA, nrow = n_signa, ncol = 2 * n_method)
    var_mat <- matrix(NA, nrow = n_signa, ncol = 2 * n_method)
    for (iter_sign in seq_len(n_signa)) {
      signa_stre <- signa_sten_vec[iter_sign]
      res_tmp <- readRDS(paste0(
        "ApplicationIII/simu3/results/a0", a0, "a1", a1,
        "sign", signa_stre, ".rds"
      ))
      res_mat[iter_sign, ] <- colMeans(res_tmp)[-c(3, 4)]
      ci_mat_l[iter_sign, ] <- apply(res_tmp, 2, quantile, 0.025)[-c(3, 4)]
      ci_mat_r[iter_sign, ] <- apply(res_tmp, 2, quantile, 0.975)[-c(3, 4)]
      var_mat[iter_sign, ] <- apply(res_tmp, 2, var)[-c(3, 4)]
    }
    df_plot <- data.frame(
      power = as.numeric(res_mat[, seq(1, 2 * n_method, by = 2)]),
      power_ci_l = as.numeric(ci_mat_l[, seq(1, 2 * n_method, by = 2)]),
      power_ci_r = as.numeric(ci_mat_r[, seq(1, 2 * n_method, by = 2)]),
      power_var = as.numeric(var_mat[, seq(1, 2 * n_method, by = 2)]),
      fdr = as.numeric(res_mat[, seq(2, 2 * n_method, by = 2)]),
      fdr_ci_l = as.numeric(ci_mat_l[, seq(2, 2 * n_method, by = 2)]),
      fdr_ci_r = as.numeric(ci_mat_r[, seq(2, 2 * n_method, by = 2)]),
      fdr_var = as.numeric(var_mat[, seq(2, 2 * n_method, by = 2)]),
      method = rep(method_vec, each = n_signa),
      signa = rep(signa_sten_vec, n_method),
      signal_den = rep(a0, n_signa * n_method),
      inform = rep(a1, n_signa * n_method)
    )
    df_plot$method <- factor(df_plot$method, levels = method_vec)
    df_plot_all <- rbind(df_plot_all, df_plot)
  }
}
df_plot_all$signal_den[which(df_plot_all$signal_den == 3.5)] <- "Sparse Signal"
df_plot_all$signal_den[which(df_plot_all$signal_den == 2.5)] <- "Medium Signal"
df_plot_all$signal_den[which(df_plot_all$signal_den == 1.5)] <- "Dense Signal"
df_plot_all$signal_den <- factor(df_plot_all$signal_den,
  levels = c("Sparse Signal", "Medium Signal", "Dense Signal")
)
df_plot_all$inform[which(df_plot_all$inform == 1.5)] <- "Less informative"
df_plot_all$inform[which(df_plot_all$inform == 2)] <- "Moderate informative"
df_plot_all$inform[which(df_plot_all$inform == 2.5)] <- "Strong informative"
df_plot_all$inform <- factor(df_plot_all$inform,
  levels = c("Less informative", "Moderate informative", "Strong informative")
)
## power plot
ppower <- ggplot(data = df_plot_all, aes(x = signa, y = power, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  # geom_errorbar(aes(ymin = power_ci_l, ymax = power_ci_r), width = 0.05) +
  geom_point(size = 0.5) +
  xlab("Signal strength") +
  ylab("Power") +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu3/results/figures/power.pdf"), width = 15, height = 10)
print(ppower)
dev.off()

## Power variance
ppower_var <- ggplot(data = df_plot_all, aes(x = signa, y = power_var, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 0.5) +
  xlab("Signal strength") +
  ylab("Variance of Power") +
  scale_y_continuous(limits = c(0, 0.012)) +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu3/results/figures/power_var.pdf"), width = 15, height = 10)
print(ppower_var)
dev.off()

## FDR plot
pfdr <- ggplot(data = df_plot_all, aes(x = signa, y = fdr, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 0.5) +
  # geom_errorbar(aes(ymin = fdr_ci_l, ymax = fdr_ci_r), width = 0.05) +
  xlab("Signal strength") +
  ylab("FDR") +
  # scale_y_continuous(limits = c(0, 0.2)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", linewidth = 0.5) +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu3/results/figures/fdr.pdf"), width = 15, height = 10)
print(pfdr)
dev.off()

## FDR variance
pfdr_var <- ggplot(data = df_plot_all, aes(x = signa, y = fdr_var, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 0.5) +
  xlab("Signal strength") +
  ylab("Variance of FDP") +
  scale_y_continuous(limits = c(0, 0.0025)) +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu3/results/figures/fdr_var.pdf"), width = 15, height = 10)
print(pfdr_var)
dev.off()