devtools::load_all()
data_matrix = readRDS("large_test/yeast_data.rds")

global_na = c(NA, Inf, 0)
perspective = "global"
scale_max = TRUE
diag_good = TRUE
include_only = NULL
check_timing = FALSE

data_matrix = data_matrix[1:200, ]


# checking timing
library(ICIKendallTau)
library(furrr)
plan(multicore)
tmp_vals = matrix(rnorm(1000 * 100), 1000, 100)
colnames(tmp_vals) = paste0("s", seq_len(ncol(tmp_vals)))
rownames(tmp_vals) = paste0("f", seq_len(nrow(tmp_vals)))

out_cor = ici_kendalltau(tmp_vals)
out_cor$run_time / 3600

# reinstall
library(ICIKendallTau)
library(furrr)
plan(multicore)
tmp_vals = matrix(rnorm(1000 * 100), 1000, 100)
colnames(tmp_vals) = paste0("s", seq_len(ncol(tmp_vals)))
rownames(tmp_vals) = paste0("f", seq_len(nrow(tmp_vals)))

out_cor = ici_kendalltau(tmp_vals, return_matrix = FALSE)
out_cor$run_time / 3600
