devtools::load_all()
data_matrix = readRDS("large_test/yeast_data.rds")

global_na = c(NA, Inf, 0)
perspective = "global"
scale_max = TRUE
diag_good = TRUE
include_only = NULL
check_timing = FALSE

data_matrix = data_matrix[1:200, ]
