#' information-content-informed kendall tau
#' 
#' Given a data-matrix, computes the information-content-informed (ICI) Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix samples are rows, features are columns
#' @param exclude_na should NA values be treated as NA?
#' @param exclude_inf should Inf values be treated as NA?
#' @param exclude_0 should zero values be treated as NA?
#' @param zero_value what is the actual zero value?
#' @param perspective how to treat missing data in denominator and ties, see details
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_not_na should the diagonal entries reflect how many entries in the sample were "good"?
#' 
#' @details For more details, see the ICI-Kendall-tau vignette:
#'   \href{../doc/ici-kendalltau.html}{\code{vignette("ici-kendalltau", package = "visualizationQualityControl")}}
#' 
#' @return numeric
#' @keywords internal
#' 
ici_kendalltau_ref = function(data_matrix, 
                             exclude_na = TRUE, 
                             exclude_inf = TRUE, 
                             exclude_0 = TRUE, 
                             zero_value = 0, 
                             perspective = "global",
                             scale_max = TRUE,
                             diag_good = TRUE,
                             progress = FALSE){
  
  # assume row-wise (because that is what the description states), so need to transpose
  # because `cor` actually does things columnwise.
  data_matrix <- t(data_matrix)
  na_loc <- matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc
  
  if (exclude_na) {
    na_loc <- is.na(data_matrix)
  }
  
  if (exclude_inf) {
    inf_loc <- is.infinite(data_matrix)
  }
  
  if (exclude_0) {
    zero_loc <- data_matrix == zero_value
  }
  
  exclude_loc <- na_loc | zero_loc | inf_loc
  
  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  # set everything to NA and let R take care of it
  
  cor_matrix = matrix(NA, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
  ntotal = 0
  for (icol in seq(1, ncol(exclude_data))) {
    for (jcol in seq(icol, ncol(exclude_data))) {
      cor_matrix[icol, jcol] = cor_matrix[jcol, icol] = ici_kt_pairs(exclude_data[, icol], exclude_data[, jcol], perspective = perspective)
      # ntotal = ntotal + 1
      # message(ntotal)
    }
  }
  
  # calculate the max-cor value for use in scaling across multiple comparisons
  n_observations = nrow(exclude_data)
  n_na = sort(colSums(exclude_loc))
  m_value = floor(sum(n_na[1:2]) / 2)
  n_m = n_observations - m_value
  max_cor_denominator = choose(n_m, 2) + n_observations * m_value
  max_cor_numerator = choose(n_m, 2) + n_observations * m_value + choose(m_value, 2)
  max_cor = max_cor_denominator / max_cor_numerator
  
  if (scale_max) {
    out_matrix = cor_matrix / max_cor
  } else {
    out_matrix = cor_matrix
  }
  
  if (diag_good) {
    n_good = colSums(!exclude_loc)
    diag(out_matrix) = n_good / max(n_good)
  }
  
  return(list(cor = out_matrix, raw = cor_matrix, keep = t(!exclude_loc)))
}


#' information-content-informed kendall tau
#' 
#' Given a data-matrix, computes the information-theoretic Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix samples are rows, features are columns
#' @param exclude_na should NA values be treated as NA?
#' @param exclude_inf should Inf values be treated as NA?
#' @param exclude_0 should zero values be treated as NA?
#' @param zero_value what is the actual zero value?
#' @param perspective how to treat missing data in denominator and ties, see details
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_not_na should the diagonal entries reflect how many entries in the sample were "good"?
#' @param check_timings should we try to estimate run time for full dataset? (default is FALSE)
#' 
#' @details For more details, see the ICI-Kendall-tau vignette:
#'   \href{../doc/ici-kendalltau.html}{\code{vignette("ici-kendalltau", package = "visualizationQualityControl")}}
#'   
#'   When \code{check_timings = TRUE}, 5 random pairwise comparisons will be run to generate timings on a single core, and then estimates of how long the full set will take are calculated. The data is returned as a data.frame, and will be on the low side, but it should provide you with a good idea of how long your data will take.
#' 
#' @return numeric
#' @export
#' 
ici_kendalltau = function(data_matrix, 
                             exclude_na = TRUE, 
                             exclude_inf = TRUE, 
                             exclude_0 = TRUE, 
                             zero_value = 0, 
                             perspective = "global",
                             scale_max = TRUE,
                             diag_good = TRUE,
                             check_timing = FALSE){
  
  # assume row-wise (because that is what the description states), so need to transpose
  data_matrix <- t(data_matrix)
  na_loc <- matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc
  
  if (exclude_na) {
    na_loc <- is.na(data_matrix)
  }
  
  if (exclude_inf) {
    inf_loc <- is.infinite(data_matrix)
  }
  
  if (exclude_0) {
    zero_loc <- data_matrix == zero_value
  }
  
  exclude_loc <- na_loc | zero_loc | inf_loc
  
  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  n_sample = ncol(exclude_data)
  # set everything to NA and let R take care of it
  
  if (suppressWarnings(require("furrr", quietly = TRUE))) {
    ncore = future::nbrOfWorkers()
    names(ncore) = NULL
    split_fun = furrr::future_map
  } else {
    ncore = 1
    split_fun = purrr::map
  }
  
  
  pairwise_comparisons = combn(n_sample, 2)
  
  if (!diag_good) {
    extra_comparisons = matrix(rep(seq(1, n_sample), each = 2), nrow = 2, ncol = n_sample, byrow = FALSE)
    pairwise_comparisons = cbind(pairwise_comparisons, extra_comparisons)
  }
  
  n_todo = ncol(pairwise_comparisons)
  n_each = ceiling(n_todo / ncore)
  
  if (check_timing) {
    sample_compare = sample(ncol(pairwise_comparisons), 5)
    tmp_pairwise = pairwise_comparisons[, sample_compare]
    
    run_tmp = check_icikt_timing(exclude_data, tmp_pairwise, perspective, n_todo, ncore)
    return(run_tmp)
  }
  
  split_comparisons = vector("list", ncore)
  start_loc = 1
  
  for (isplit in seq_along(split_comparisons)) {
    stop_loc = min(start_loc + n_each, n_todo)
    
    split_comparisons[[isplit]] = pairwise_comparisons[, start_loc:stop_loc, drop = FALSE]
    start_loc = stop_loc + 1
    
    if (start_loc > n_todo) {
      break()
    }
  }
  
  null_comparisons = purrr::map_lgl(split_comparisons, is.null)
  split_comparisons = split_comparisons[!null_comparisons]
  
  do_split = function(do_comparisons, exclude_data, perspective) {
    #seq_range = seq(in_range[1], in_range[2])
    #print(seq_range)
    tmp_cor = matrix(0, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
    rownames(tmp_cor) = colnames(tmp_cor) = colnames(exclude_data)
    tmp_pval = tmp_cor
    
    for (icol in seq(1, ncol(do_comparisons))) {
      iloc = do_comparisons[1, icol]
      jloc = do_comparisons[2, icol]
      ici_res = ici_kt(exclude_data[, iloc], exclude_data[, jloc], perspective = perspective)
      tmp_cor[iloc, jloc] = tmp_cor[jloc, iloc] = ici_res["tau"]
      tmp_pval[iloc, jloc] = tmp_pval[jloc, iloc] = ici_res["pvalue"]
    }
    list(cor = tmp_cor, pval = tmp_pval)
  }
  # we record how much time is actually spent doing ICI-Kt
  # itself, as some of the other operations will add a bit of time
  # 
  t1 = Sys.time()
  split_cor = split_fun(split_comparisons, do_split, exclude_data, perspective)
  t2 = Sys.time()
  t_diff = as.numeric(difftime(t2, t1, units = "secs"))

  cor_matrix = matrix(0, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
  pvalue_matrix = cor_matrix
  for (isplit in split_cor) {
    cor_matrix = cor_matrix + isplit$cor
    pvalue_matrix = pvalue_matrix + isplit$pval
  }
  
  
  # calculate the max-cor value for use in scaling across multiple comparisons
  n_observations = nrow(exclude_data)
  n_na = sort(colSums(exclude_loc))
  m_value = floor(sum(n_na[1:2]) / 2)
  n_m = n_observations - m_value
  max_cor_denominator = choose(n_m, 2) + n_observations * m_value
  max_cor_numerator = choose(n_m, 2) + n_observations * m_value + choose(m_value, 2)
  max_cor = max_cor_denominator / max_cor_numerator
  
  if (scale_max) {
    out_matrix = cor_matrix / max_cor
  } else {
    out_matrix = cor_matrix
  }
  
  if (diag_good) {
    n_good = colSums(!exclude_loc)
    diag(out_matrix) = n_good / max(n_good)
  }
  
  return(list(cor = out_matrix, raw = cor_matrix, pval = pvalue_matrix, keep = t(!exclude_loc),
              run_time = t_diff))
}


check_icikt_timing = function(exclude_data, tmp_pairwise, perspective, n_todo, ncore){
  t_start = Sys.time()
  for (icol in seq_len(ncol(tmp_pairwise))) {
    iloc = tmp_pairwise[1, icol]
    jloc = tmp_pairwise[2, icol]
    tmp_val = ici_kt(exclude_data[, iloc], exclude_data[, jloc], perspective = perspective)
  }
  t_stop = Sys.time()
  t_total = as.numeric(difftime(t_stop, t_start, units = "secs"))
  n_comp = ncol(tmp_pairwise)
  
  t_each = t_total / n_comp
  
  t_theoretical = t_each * n_todo
  t_cores = t_theoretical / ncore
  
  data.frame(which =
               c("n_tested",
                 "n_todo",
                 "time_tested",
                 "time_single",
                 "time_all",
                 "time_across_cores",
                 "time_minutes",
                 "time_hours",
                 "time_days"),
             value = c(n_comp,
                n_todo,
                t_total,
                t_each,
                t_theoretical,
                t_cores,
                t_cores / 60,
                t_cores / (60 * 60),
                t_cores / (60 * 60 * 60)))
  
}
