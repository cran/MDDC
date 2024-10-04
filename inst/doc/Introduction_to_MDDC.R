## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("MDDC")

## -----------------------------------------------------------------------------
library(MDDC)

## -----------------------------------------------------------------------------
data("statin49")
head(statin49)

data("statin49_AE_idx")
head(statin49_AE_idx)

## -----------------------------------------------------------------------------
set.seed(42)
test1 <- mddc_boxplot(
  contin_table = statin49,
  col_specific_cutoff = T,
  separate = T,
  if_col_cor = F,
  cor_lim = 0.8,
  coef = 1.5
)

## -----------------------------------------------------------------------------
head(test1$boxplot_signal)[, 1:5]

## -----------------------------------------------------------------------------
round(head(test1$corr_signal_pval)[, 1:5], digits = 3)

## -----------------------------------------------------------------------------
set.seed(42)
find_optimal_coef(
  contin_table = statin49,
  n_sim = 1000,
  target_fdr = 0.05,
  grid = 0.1,
  col_specific_cutoff = TRUE,
  exclude_small_count = TRUE
)

## -----------------------------------------------------------------------------
set.seed(42)
test2 <- mddc_mc(
  contin_table = statin49,
  quantile = 0.95,
  rep = 10000,
  exclude_same_drug_class = T,
  col_specific_cutoff = T,
  separate = T,
  if_col_cor = F,
  cor_lim = 0.8
)

## -----------------------------------------------------------------------------
test3 <- report_drug_AE_pairs(
  contin_table = statin49,
  contin_table_signal = test2$mc_signal
)

head(test3)

## -----------------------------------------------------------------------------
report_drug_AE_pairs(
  contin_table = statin49,
  contin_table_signal = test2$corr_signal_pval < 0.05
)

## -----------------------------------------------------------------------------
# create a matrix indicating signal strength
sig_mat <- matrix(1,
  nrow = nrow(statin49),
  ncol = ncol(statin49)
)

# assign (Rhabdomyolysis, Atorvastain) as a signal
# with a signal strength 4
sig_mat[1, 1] <- 4

## -----------------------------------------------------------------------------
head(statin49_AE_idx)

## -----------------------------------------------------------------------------
sim_dat <- generate_contin_table_with_clustered_AE(
  contin_table = statin49,
  n_rep = 3,
  AE_idx = statin49_AE_idx,
  rho = 0.5,
  signal_mat = sig_mat,
  seed = 42
)

## -----------------------------------------------------------------------------
test5 <- mddc_mc(sim_dat[[1]], seed = 42)
report_drug_AE_pairs(
  contin_table = sim_dat[[1]],
  contin_table_signal = test5$mc_signal
)

## ----fig.width=6, fig.height=10-----------------------------------------------
plot_heatmap(test2$mc_pval[-50, -7])

