source("ApplicationIII/realdata/utility.R")
#############################################################################
alpha_vec <- seq(0.01, 0.1, by = 0.01)
n_alpha <- length(alpha_vec)
n_method <- 5
# The downloaded data were already normalized
# Rare OTUs occured in less than or equal to 20 subjects were excluded (>0.2% prevalence) in the data
meta.dat <- readRDS(file = "ApplicationIII/realdata/datasets/MWAS/amgut.meta.dat.rds")
otu.tab <- readRDS(file = "ApplicationIII/realdata/datasets/MWAS/amgut.otu.dat.rds")
otu.name <- readRDS(file = "ApplicationIII/realdata/datasets/MWAS/amgut.otu.name.rds")

# # We select subjects from United States and adults
country_residence_use <- unique(meta.dat$country_residence)[-1]
ind <- meta.dat$country_residence %in% country_residence_use &
  meta.dat$age_cat %in% c("teen", "20s", "30s", "40s", "50s", "60s", "70+") & meta.dat$sex %in% c("female", "male")
# ind <- meta.dat$country_residence %in% country_residence_use &
#   meta.dat$age_cat %in% c('20s', '30s', '40s', '50s', '60s', '70+') & meta.dat$sex %in% c('female', 'male')
otu.tab <- otu.tab[, ind]
meta.dat <- meta.dat[ind, ]
sex <- meta.dat$sex
sex <- factor(sex)

# Further discard OTU occuring in less than 10 subjects
ind <- rowSums(otu.tab != 0) >= 10
otu.tab <- otu.tab[ind, ]
otu.name <- otu.name[ind, ]

# Wilcox rank sum test
pvals <- apply(otu.tab, 1, function(x) wilcox.test(x ~ sex)$p.value)
x <- log(rowSums(otu.tab))
# x <- log(rowSums(otu.tab != 0))
mwas <- data.frame(pvalue = pvals, covariate = x)
saveRDS(mwas, "ApplicationIII/realdata/datasets/mwas.p.value.rds")

pvalue_vec <- pvals
x_covariate <- x

set.seed(1991)
outres <- list()
for (iter_alpha in seq_len(n_alpha)) {
  alpha_use <- alpha_vec[iter_alpha]
  outres[[iter_alpha]] <- simu_fun(alpha_use)
}
res_mat <- matrix(unlist(outres), nrow = n_method, ncol = n_alpha)
print(res_mat)
saveRDS(res_mat, paste0("ApplicationIII/realdata/results/mwas_res.rds"))
