set.seed(1)
library(ggplot2)
library(gridExtra)
n_simu <- 100
a0_vec <- c(3.5, 2.5, 1.5)
a1_vec <- c(1.5, 2, 2.5)
signa_sten_vec <- seq(1.9, 3.4, by = 0.15)
n_signa <- length(signa_sten_vec)
method_vec <- c(
  "BH", "IHW", "adaPT", "SABHA", "eBH_FBC"
)
n_method <- length(method_vec)
dir.create("ApplicationIII/simu2/results/figures", recursive = TRUE, showWarnings = FALSE)

######################## plot one figures ########################
df_plot_all <- NULL
for (a0 in a0_vec) {
  for (a1 in a1_vec) {
    ## load results
    res_mat <- matrix(NA, nrow = n_signa, ncol = 2 * n_method)
    sd_mat <- matrix(NA, nrow = n_signa, ncol = 2 * n_method)
    for (iter_sign in seq_len(n_signa)) {
      signa_stre <- signa_sten_vec[iter_sign]
      res_tmp <- readRDS(paste0(
        "ApplicationIII/simu2/results/a0", a0, "a1", a1,
        "sign", signa_stre, ".rds"
      ))
      res_mat[iter_sign, ] <- colMeans(res_tmp)
      sd_mat[iter_sign, ] <- apply(res_tmp, 2, sd)
    }
    df_plot <- data.frame(
      power = as.numeric(res_mat[, seq(1, 2 * n_method, by = 2)]),
      sd_power = as.numeric(sd_mat[, seq(1, 2 * n_method, by = 2)]),
      fdr = as.numeric(res_mat[, seq(2, 2 * n_method, by = 2)]),
      sd_fdr = as.numeric(sd_mat[, seq(2, 2 * n_method, by = 2)]),
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
  geom_point(size = 1) +
  xlab("Signal strength") +
  ylab("Power") +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu2/results/figures/power.pdf"), width = 20, height = 10)
print(ppower)
dev.off()
## FDR plot
pfdr <- ggplot(data = df_plot_all, aes(x = signa, y = fdr, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 1.5) +
  xlab("Signal strength") +
  ylab("FDR") +
  geom_hline(yintercept = 0.1, linetype = "solid", linewidth = 1.5) +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu2/results/figures/fdr.pdf"), width = 20, height = 10)
print(pfdr)
dev.off()
## power variance
ppower_sd <- ggplot(data = df_plot_all, aes(x = signa, y = sd_power, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 1) +
  xlab("Signal strength") +
  ylab("Power standard error") +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu2/results/figures/power_sd.pdf"), width = 20, height = 10)
print(ppower_sd)
dev.off()
## FDR variance
pfdr_sd <- ggplot(data = df_plot_all, aes(x = signa, y = sd_fdr, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 1.5) +
  xlab("Signal strength") +
  ylab("FDR") +
  theme_bw(base_size = 22) +
  theme(legend.position = "right") +
  facet_grid(signal_den ~ inform, scales = "free_y")
## save plots
pdf(paste0("ApplicationIII/simu2/results/figures/fdr_sd.pdf"), width = 20, height = 10)
print(pfdr_sd)
dev.off()
