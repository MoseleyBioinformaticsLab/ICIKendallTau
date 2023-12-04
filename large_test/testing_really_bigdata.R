<<<<<<< Updated upstream
devtools::load_all()
=======
<<<<<<< Updated upstream
devtools::load_all()
=======

>>>>>>> Stashed changes
>>>>>>> Stashed changes
data_matrix = readRDS("large_test/yeast_data.rds")

global_na = c(NA, Inf, 0)
perspective = "global"
scale_max = TRUE
diag_good = TRUE
include_only = NULL
check_timing = FALSE

data_matrix = data_matrix[1:200, ]


# checking timing
set.seed(1234)
<<<<<<< Updated upstream
=======
<<<<<<< Updated upstream
>>>>>>> Stashed changes
library(ICIKendallTau)
library(furrr)
plan(multicore)
tmp_vals = matrix(rnorm(1000 * 100), 1000, 100)
colnames(tmp_vals) = paste0("s", seq_len(ncol(tmp_vals)))
rownames(tmp_vals) = paste0("f", seq_len(nrow(tmp_vals)))

library(tictoc)
tic()
out_cor = ici_kendalltau(tmp_vals)
<<<<<<< Updated upstream
=======
=======
devtools::load_all()
library(furrr)
library(tictoc)
plan(multicore)
tmp_vals = matrix(rnorm(2000 * 100), 2000, 100)
colnames(tmp_vals) = paste0("s", seq_len(ncol(tmp_vals)))
rownames(tmp_vals) = paste0("f", seq_len(nrow(tmp_vals)))


tic()
out_cor = ici_kendalltau(tmp_vals, return_matrix = FALSE)
>>>>>>> Stashed changes
>>>>>>> Stashed changes
toc()
out_cor$run_time

# reinstall
library(ICIKendallTau)
library(furrr)
plan(multicore)
tmp_vals = matrix(rnorm(1000 * 100), 1000, 100)
colnames(tmp_vals) = paste0("s", seq_len(ncol(tmp_vals)))
rownames(tmp_vals) = paste0("f", seq_len(nrow(tmp_vals)))

out_cor = ici_kendalltau(tmp_vals, return_matrix = FALSE)
out_cor$run_time / 3600

# testing another approache
getData <- function(i){
  x <- data.frame(x=rnorm(100000), y=rnorm(100000), z=rnorm(100000))
  x
}

large_list = purrr::map(seq_len(20), getData)

unlist_it = function(in_list){
  list2DF(
    lapply(setNames(seq_along(in_list[[1]]), names(in_list[[1]])), \(i){
      unlist(lapply(in_list, `[[`, i), FALSE, FALSE)
    })
  )               
}

out_df = unlist_it(large_list)

microbenchmark::microbenchmark(
  plist = purrr::list_rbind(large_list),
  ut = unlist_it(large_list),
  times = 100
)
