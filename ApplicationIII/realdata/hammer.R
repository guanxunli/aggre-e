# ---- Setup ----
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "Biobase"))
library(DESeq2)
library(Biobase)
library(dplyr)

# ---- 1) Load hammer.eset from your local RData ----
# If your file path is different, change "hammer_eset.RData" accordingly.
e <- new.env(parent = emptyenv())
load("ApplicationIII/realdata/datasets/hammer_eset.RData", envir = e) # expects object named 'hammer.eset'
hammer.eset <- e$hammer.eset
stopifnot(inherits(hammer.eset, "ExpressionSet"))

# ---- 2) Extract counts and sample annotation ----
cts <- exprs(hammer.eset) # genes x samples (counts)
pd <- pData(hammer.eset) # sample metadata

# ---- 3) Standardize factors for DESeq2 design ~ time ----
# protocol: control vs SNL
pd$protocol <- factor(
  ifelse(grepl("SNL", pd$protocol, ignore.case = TRUE), "snl", "c"),
  levels = c("c", "snl")
)

# time: 2 weeks vs 2 months (handle label variants like "2months", "2 m")
time_std <- tolower(gsub("\\s+", "", as.character(pd$Time)))
pd$time <- factor(ifelse(grepl("^2m", time_std), "2m", "2w"),
  levels = c("2w", "2m")
)

# Ensure sample order matches
stopifnot(identical(colnames(cts), rownames(pd)))

# ---- 4) Build DESeq2 object & run DE ----
dds <- DESeqDataSetFromMatrix(
  countData = round(cts), # DESeq2 expects integer counts
  colData   = pd,
  design    = ~time
)
dds <- DESeq(dds)

# ---- 5) Choose contrasts and compute results ----

# Time effect: 2 months vs 2 weeks
res_early <- as.data.frame(results(dds))

# ---- 6) Match your airway-style extraction ----
pvals <- res_early$pvalue
ind <- which(!is.na(pvals))
pvals <- pvals[ind]

x <- log(res_early$baseMean)
x <- x[ind]

###################################################
source("ApplicationIII/realdata/utility.R")
pvalue_vec <- pvals
x_covariate <- x
alpha_vec <- seq(0.01, 0.1, by = 0.01)
n_alpha <- length(alpha_vec)
n_method <- 5
outres <- list()

set.seed(1991)
for (iter_alpha in seq_len(n_alpha)) {
  print(iter_alpha)
  alpha_use <- alpha_vec[iter_alpha]
  outres[[iter_alpha]] <- simu_fun(alpha_use)
}
res_mat <- matrix(unlist(outres), nrow = n_method, ncol = n_alpha)
print(res_mat)

saveRDS(res_mat, paste0("ApplicationIII/realdata/results/hammer_res.rds"))
