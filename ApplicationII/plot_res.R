set.seed(1)
library(ggplot2)
library(gridExtra)
dir.create("ApplicationII/results", recursive = TRUE, showWarnings = FALSE)

######################## first figure ########################
ratio <- 2.5e-1
signa_stre_vec <- c(2.2, 2.35, 2.5, 2.65, 2.8, 2.95, 3.1, 3.25, 3.5)
n_signa <- length(signa_stre_vec)
res_mat <- matrix(NA, nrow = n_signa, ncol = 12)
for (iter_sign in seq_len(n_signa)) {
  signa_stre <- signa_stre_vec[iter_sign]
  res_tmp <- readRDS(paste0(
    "ApplicationII/results/bhbc_ratio",
    ratio, "sign", signa_stre, ".rds"
  ))
  stopifnot(ncol(res_tmp) == 12) # Rerun simulations if these are old 10-column results.
  res_mat[iter_sign, ] <- colMeans(res_tmp)
}
df_plot1 <- data.frame(
  power = as.numeric(res_mat[, c(1, 3, 5, 7, 9)]),
  fdr = as.numeric(res_mat[, c(2, 4, 6, 8, 10)]),
  method = rep(c("BH", "BC", "eBH_Ave", "eBH_Ada", "fast_eBH_Ada"), each = n_signa),
  signa = rep(signa_stre_vec, 5)
)
df_plot1$method <- factor(df_plot1$method, levels = c("BH", "BC", "eBH_Ave", "eBH_Ada", "fast_eBH_Ada"))
ppower1 <- ggplot(data = df_plot1, aes(x = signa, y = power, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 2) +
  xlab("Mean signal strength") +
  ylab("Power") +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")
pfdr1 <- ggplot(data = df_plot1, aes(x = signa, y = fdr, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 2) +
  xlab("Mean signal strength") +
  ylab("FDR") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")
## generate legend
p_legend <- ggplot(data = df_plot1, aes(x = signa, y = power, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 3) +
  xlab("Signal strength") +
  ylab("Power") +
  theme_bw(base_size = 22) +
  theme(legend.position = "bottom")
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
shared_legend <- extract_legend(p_legend)

######################## second figure ########################
ratio <- 0.5
signa_stre_vec <- c(1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2)
n_signa <- length(signa_stre_vec)
res_mat <- matrix(NA, nrow = n_signa, ncol = 12)
for (iter_sign in seq_len(n_signa)) {
  signa_stre <- signa_stre_vec[iter_sign]
  res_tmp <- readRDS(paste0(
    "ApplicationII/results/bcbh_ratio",
    ratio, "sign", signa_stre, ".rds"
  ))
  stopifnot(ncol(res_tmp) == 12) # Rerun simulations if these are old 10-column results.
  res_mat[iter_sign, ] <- colMeans(res_tmp)
}
df_plot2 <- data.frame(
  power = as.numeric(res_mat[, c(1, 3, 5, 7, 9)]),
  fdr = as.numeric(res_mat[, c(2, 4, 6, 8, 10)]),
  method = rep(c("BH", "BC", "eBH_Ave", "eBH_Ada", "fast_eBH_Ada"), each = n_signa),
  signa = rep(signa_stre_vec, 5)
)
df_plot2$method <- factor(df_plot2$method, levels = c("BH", "BC", "eBH_Ave", "eBH_Ada", "fast_eBH_Ada"))
ppower2 <- ggplot(data = df_plot2, aes(x = signa, y = power, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 2) +
  xlab("Signal strength") +
  ylab("Power") +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")
pfdr2 <- ggplot(data = df_plot2, aes(x = signa, y = fdr, color = method)) +
  geom_line(aes(linetype = method), linewidth = 1) +
  geom_point(size = 2) +
  xlab("Signal strength") +
  ylab("FDR") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw(base_size = 22) +
  theme(legend.position = "none")

## save plots
pdf("ApplicationII/results/normal.pdf", width = 10, height = 6)
print(grid.arrange(arrangeGrob(pfdr1, pfdr2, ppower1, ppower2, nrow = 2),
  shared_legend,
  nrow = 2, heights = c(10, 2)
))
dev.off()
