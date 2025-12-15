library(DESeq2)
library(dplyr)

# ---- 1) Load Pasilla counts + sample annotation ----
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
  package = "pasilla", mustWork = TRUE
)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv",
  package = "pasilla", mustWork = TRUE
)

cts <- as.matrix(read.csv(pasCts, sep = "\t", row.names = "gene_id"))
coldata <- read.csv(pasAnno, row.names = 1, check.names = FALSE)

# Harmonize labels: drop the trailing "fb" in sample names from annotation
coldata$file <- sub("fb$", "", rownames(coldata))

# Make factors and simplify 'type' from "single-read"/"paired-end" to "single"/"paired"
coldata$condition <- factor(coldata$condition, levels = c("untreated", "treated"))
coldata$type <- factor(sub("-.*", "", coldata$type), levels = c("single", "paired"))

# Align rows of coldata to columns of counts
rownames(coldata) <- coldata$file
coldata <- coldata[colnames(cts), c("condition", "type")]

stopifnot(all(rownames(coldata) == colnames(cts))) # sanity check

# ---- 2) Build DESeq2 object & run DE ----
# Include 'type' (single vs paired) as a nuisance covariate, analogous to '~ cell + dex' in airway
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~type
) %>%
  DESeq()

# Contrast: treated vs untreated (this name is created from the design)
deRes <- as.data.frame(results(dds))
# deRes <- na.omit(deRes)

# ---- 3) Extract p-values and matched covariate like your airway snippet ----
pvals <- deRes$pvalue
ind <- which(!is.na(pvals))
pvals <- pvals[ind]

x <- log(deRes$baseMean)
x <- x[ind]


###################################################
pvalue_vec <- pvals
x_covariate <- x
alpha_vec <- seq(0.01, 0.1, by = 0.01)
n_alpha <- length(alpha_vec)
n_method <- 5
source("ApplicationIII/real_data/utility.R")

set.seed(2011)
outres <- list()
for (iter_alpha in seq_len(n_alpha)) {
  print(iter_alpha)
  alpha_use <- alpha_vec[iter_alpha]
  outres[[iter_alpha]] <- simu_fun(alpha_use)
}
res_mat <- matrix(unlist(outres), nrow = n_method, ncol = n_alpha)
print(res_mat)
saveRDS(res_mat, paste0("ApplicationIII/real_data/results/pasilla_res.rds"))
