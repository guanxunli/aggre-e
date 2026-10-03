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
out_list <- list()
for (iter in seq_len(2)) {
  out_list[[iter]] <- matrix(NA, nrow = 100, ncol = 3)
}

Y <- otu_counts
rownames(Y) <- paste0("taxa", 1:ntaxon)
colnames(Y) <- paste0("sample", 1:nsample)

cond_vec <- as.factor(ifelse(meta$DiseaseState == "nonCDI", 0, 1))
cond_vec1 <- which(cond_vec == 1)
n1 <- length(cond_vec1)
cond_vec0 <- which(cond_vec == 0)
n0 <- length(cond_vec0)

for (iter in seq_len(100)) {
  print(iter)
  
  idx1 <- sample(seq_len(n1), size = floor(n1 / 2))
  idx0 <- sample(seq_len(n0), size = floor(n0 / 2))
  idx <- c(cond_vec1[idx1], cond_vec0[idx0])
  Y1 <- Y[, idx]
  Y2 <- Y[, -idx]
  
  Z1 <- cond_vec[idx]
  Z1 <- as.data.frame(Z1)
  colnames(Z1) <- "u"
  
  Z2 <- cond_vec[-idx]
  Z2 <- as.data.frame(Z2)
  colnames(Z2) <- "u"
  
  options(warn = -1)
  res1 <- linda(feature.dat = Y1, meta.dat = Z1, formula = formula)
  pval_vec1 <- res1$output$u1$pvalue
  res2 <- linda(feature.dat = Y2, meta.dat = Z2, formula = formula)
  pval_vec2 <- res2$output$u1$pvalue
  options(warn = 0)
  
  #### run BC method
  ## joint BC
  res_joint1 <- bc_fun_cpp(p_value = pval_vec1, alpha = alpha)
  index_joint1 <- res_joint1$index_select + 1L
  res_joint2 <- bc_fun_cpp(p_value = pval_vec2, alpha = alpha)
  index_joint2 <- res_joint2$index_select + 1L
  
  for (iter_p in seq_len(3)) {
    idx_det1 <- index_joint1[phylum_vector[index_joint1] == phylum_use[iter_p]]
    idx_det2 <- index_joint2[phylum_vector[index_joint2] == phylum_use[iter_p]]
    idx_inter <- intersect(idx_det1, idx_det2)
    idx_union <- union(idx_det1, idx_det2)
    if (length(idx_union) > 0) {
      out_list[[1]][iter, iter_p] <- length(idx_inter) / length(idx_union)
    } else {
      out_list[[1]][iter, iter_p] <- 1
    }
  }

  ## adaptive eBH
  pvalue_list <- split(pval_vec1, phylum_vector)
  index_list <- split(seq_len(ntaxon), phylum_vector)
  res_weight <- ebh_apa_cpp(
    pvalue_list = pvalue_list, alpha_bc = alpha,
    alpha_ebh = alpha
  )
  index_tmp <- res_weight$index_select + 1L
  index_weight1 <- unlist(index_list)[index_tmp]
  
  pvalue_list <- split(pval_vec2, phylum_vector)
  index_list <- split(seq_len(ntaxon), phylum_vector)
  res_weight <- ebh_apa_cpp(
    pvalue_list = pvalue_list, alpha_bc = alpha,
    alpha_ebh = alpha
  )
  index_tmp <- res_weight$index_select + 1L
  index_weight2 <- unlist(index_list)[index_tmp]
  
  for (iter_p in seq_len(3)) {
    idx_det1 <- index_weight1[phylum_vector[index_weight1] == phylum_use[iter_p]]
    idx_det2 <- index_weight2[phylum_vector[index_weight2] == phylum_use[iter_p]]
    idx_inter <- intersect(idx_det1, idx_det2)
    idx_union <- union(idx_det1, idx_det2)
    if (length(idx_union) > 0) {
      out_list[[2]][iter, iter_p] <- length(idx_inter) / length(idx_union)
    } else {
      out_list[[2]][iter, iter_p] <- 1
    }
  }
}
saveRDS(out_list, "ApplicationI/results/cdi_schubert_reproducibility.rds")
lapply(out_list, function(x) round(colMeans(x), 4))
