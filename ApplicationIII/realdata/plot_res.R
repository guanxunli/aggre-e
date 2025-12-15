set.seed(1)
library(ggplot2)
library(gridExtra)
dataset_vec <- c("mwas", "pasilla", "hammer")
n_dta <- length(dataset_vec)
alpha_vec <- seq(0.01, 0.1, by = 0.01)
n_alpha <- length(alpha_vec)
method_vec <- c("BH", "IHW", "adaPT", "SABHA", "eBH_FBC")
n_method <- length(method_vec)

######################## plot figures ########################
plot_list <- list()
for (iter_dta in seq_len(n_dta)) {
  dtaname_use <- dataset_vec[iter_dta]
  res <- readRDS(paste0("ApplicationIII/real_data/results/", dtaname_use, "_res.rds"))
  df_plot <- data.frame(
    n_rej = as.numeric(res),
    method = rep(method_vec, n_alpha),
    alpha = rep(alpha_vec, each = n_method)
  )
  df_plot$method <- factor(df_plot$method, levels = method_vec)
  ## power plot
  prej <- ggplot(data = df_plot, aes(x = alpha, y = n_rej, color = method)) +
    geom_line(aes(linetype = method), linewidth = 1.5) +
    geom_point(size = 2) +
    xlab("Target FDR level") +
    ylab("Number of rejections") +
    ggtitle(dtaname_use) +
    theme_bw(base_size = 22) +
    theme(legend.position = "bottom")
  plot_list[[iter_dta]] <- prej
}


library(patchwork)
combined <- wrap_plots(plotlist = plot_list, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

# save the combined figure
dir.create("ApplicationIII/real_data/results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("ApplicationIII/real_data/results/figures/combined_three_datasets.pdf",
  combined,
  width = 16, height = 6
)
