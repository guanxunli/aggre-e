set.seed(1)
# Load required libraries
library(Rcpp)
library(RcppArmadillo)
sourceCpp("ApplicationI/utilities.cpp")
library(MicrobiomeStat)
alpha <- 0.2
dir.create("ApplicationI/results", recursive = TRUE, showWarnings = FALSE)
otu_file <- "ApplicationI/datasets/cdi_schubert_results/RDP/cdi_schubert.otu_table.100.denovo.rdp_assigned"
meta_file <- "ApplicationI/datasets/cdi_schubert_results/cdi_schubert.metadata.txt"

# Extract count matrix (features x samples)
otu <- read.delim(otu_file, header=TRUE, sep="\t", row.names=1, comment.char="", stringsAsFactors=FALSE)
meta <- read.table(meta_file, sep="\t")
colnames(meta) <- meta[1, ]
meta <- meta[-1, ]
index_col <- which(colnames(otu) %in% meta$sample_id == TRUE)
otu <- otu[, index_col]
index_row <- which(meta$sample_id %in% colnames(otu) == TRUE)
meta <- meta[index_row, ]
idx <- match(colnames(otu), meta$sample_id)
meta <- meta[idx, ]
# Parse taxonomic information from OTU identifiers to get phylum for each OTU
tax_split <- strsplit(rownames(otu), ";")
phylum_vec <- sapply(tax_split, function(x) {                      # extract the phylum component
  phylum_field <- x[2]  # second field is p__<PhylumName>
  sub("^p__*", "", phylum_field)  # remove "p__" prefix to get phylum name
})
phylum_use <- c("Bacteroidetes", "Firmicutes", "Proteobacteria")
n_group <- length(phylum_use)
taxa_use <- which(phylum_vec %in% phylum_use)
phylum_vector <- phylum_vec[taxa_use]
otu_counts <- otu[taxa_use, ]
# filtering
index_use <- which(apply(otu_counts, 1, function(x) {
  sum(x != 0)
}) >= 10)
otu_counts <- otu_counts[index_use, ]
ntaxon <- nrow(otu_counts)
nsample <- ncol(otu_counts)
phylum_vector <- phylum_vector[index_use]
## run linda
formula <- "~u"
Y <- otu_counts
rownames(Y) <- paste0("taxa", 1:ntaxon)
colnames(Y) <- paste0("sample", 1:nsample)
Z <- as.factor(ifelse(meta$DiseaseState == "nonCDI", 0, 1))
Z <- as.data.frame(Z)
colnames(Z) <- "u"
options(warn = -1)
res <- linda(feature.dat = Y, meta.dat = Z, formula = formula)
pval_vec <- res$output$u1$pvalue
options(warn = 0)

#### run BC method
## joint BC
res_joint <- bc_fun_cpp(p_value = pval_vec, alpha = alpha)
index_joint <- res_joint$index_select + 1L
table(phylum_vector[index_joint])

## average eBH
eval_ave <- numeric(ntaxon)
for (iter_group in seq_len(n_group)) {
  index_tmp <- which(phylum_vector == phylum_use[iter_group])
  p_tmp <- pval_vec[index_tmp]
  res_tmp <- bc_fun_cpp(p_value = p_tmp, alpha = alpha)
  eval_ave[index_tmp] <- res_tmp$evalue
}
index_ave <- ebh_fun_cpp(eval_ave, alpha = alpha) + 1L
table(phylum_vector[index_ave])

## weighted eBH
eval_weight <- numeric(ntaxon)
for (iter_group in seq_len(n_group)) {
  index_tmp <- which(phylum_vector == phylum_use[iter_group])
  n_tmp <- length(index_tmp)
  p_tmp <- pval_vec[index_tmp]
  res_tmp <- bc_fun_cpp(p_value = p_tmp, alpha = alpha)
  eval_weight[index_tmp] <- ntaxon / (n_group * n_tmp) * res_tmp$evalue
}
index_weight <- ebh_fun_cpp(eval_weight, alpha = alpha) + 1L
table(phylum_vector[index_weight])

## adaptive eBH
pvalue_list <- split(pval_vec, phylum_vector)
index_list <- split(seq_len(ntaxon), phylum_vector)
res_weight <- ebh_apa_cpp(
  pvalue_list = pvalue_list, alpha_bc = alpha,
  alpha_ebh = alpha
)
index_tmp <- res_weight$index_select + 1L
index_weight <- unlist(index_list)[index_tmp]
table(phylum_vector[index_weight])
