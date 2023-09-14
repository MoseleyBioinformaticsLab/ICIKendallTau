#' information-content-informed kendall tau
#' 
#' Given a data-matrix, computes the information-content-informed (ICI) Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix samples are rows, features are columns
#' @param global_na what values should be treated as missing (NA)?
#' @param zero_value what is the actual zero value?
#' @param perspective how to treat missing data in denominator and ties, see details
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_good should the diagonal entries reflect how many entries in the sample were "good"?
#' @param progress should progress be displayed.
#' 
#' @details For more details, see the ICI-Kendall-tau vignette:
#'   \href{../doc/ici-kendalltau.html}{\code{vignette("ici-kendalltau", package = "ICIKendallTau")}}
#' 
#' @return numeric
#' @keywords internal
#' 
ici_kendalltau_ref = function(data_matrix, 
                             global_na = c(NA, Inf, 0),
                             zero_value = 0, 
                             perspective = "global",
                             scale_max = TRUE,
                             diag_good = TRUE,
                             progress = FALSE){
  
  # assume row-wise (because that is what the description states), so need to transpose
  # because `cor` actually does things columnwise.
  data_matrix <- t(data_matrix)
  exclude_loc = matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  
  # Actual NA and Inf values are special cases, so we do
  # this very specifically.
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
    
  }
  
  # now that we've done the NA and Inf values, we can go
  # ahead and take care of the rest.
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }
  
  
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

missing_either = function(in_x, in_y){
  not_in_one = sum(in_x | in_y)
  not_in_one
}

#' pairwise completeness
#' 
#' Calculates the completeness between any two samples using "or", is an
#' entry missing in either X "or" Y.
#' 
#' @param data_matrix samples are rows, features are columns
#' @param global_na globally, what should be treated as NA?
#' @param zero_value what is the actual zero value?
#' 
#' @export
#' 
#' @return matrix of degree of completeness
pairwise_completeness = function(data_matrix,
                                global_na = c(NA, Inf, 0),
                                zero_value = 0){
  
  data_matrix <- t(data_matrix)
  exclude_loc = matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  
  # Actual NA and Inf values are special cases, so we do
  # this very specifically.
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
    
  }
  
  # now that we've done the NA and Inf values, we can go
  # ahead and take care of the rest.
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }
  
  n_sample = ncol(exclude_loc)
  
  if ("furrr" %in% utils::installed.packages()) {
    ncore = future::nbrOfWorkers()
    names(ncore) = NULL
    split_fun = furrr::future_map
  } else {
    ncore = 1
    split_fun = purrr::map
  }
  
  
  pairwise_comparisons = utils::combn(n_sample, 2)
  
  extra_comparisons = matrix(rep(seq(1, n_sample), each = 2), nrow = 2, ncol = n_sample, byrow = FALSE)
  pairwise_comparisons = cbind(pairwise_comparisons, extra_comparisons)
  
  n_todo = ncol(pairwise_comparisons)
  n_each = ceiling(n_todo / ncore)
  
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
  
  do_split = function(do_comparisons, exclude_loc) {
    #seq_range = seq(in_range[1], in_range[2])
    #print(seq_range)
    tmp_missing = matrix(0, nrow = ncol(exclude_loc), ncol = ncol(exclude_loc))
    rownames(tmp_missing) = colnames(tmp_missing) = colnames(exclude_loc)
    
    for (icol in seq(1, ncol(do_comparisons))) {
      iloc = do_comparisons[1, icol]
      jloc = do_comparisons[2, icol]
      missing_res = missing_either(exclude_loc[, iloc], exclude_loc[, jloc])
      tmp_missing[iloc, jloc] = tmp_missing[jloc, iloc] = missing_res
    }
    tmp_missing
  }
  
  split_missing = split_fun(split_comparisons, do_split, exclude_loc)
  
  missing_matrix = matrix(0, nrow = ncol(exclude_loc), ncol = ncol(exclude_loc))
  rownames(missing_matrix) = colnames(missing_matrix) = colnames(exclude_loc)
  for (isplit in split_missing) {
    missing_matrix = missing_matrix + isplit
  }
  missing_matrix = missing_matrix / nrow(exclude_loc)
  1 - missing_matrix
}

#' information-content-informed kendall tau
#' 
#' Given a data-matrix, computes the information-theoretic Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix samples are rows, features are columns
#' @param global_na globally, what should be treated as NA?
#' @param perspective how to treat missing data in denominator and ties, see details
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_good should the diagonal entries reflect how many entries in the sample were "good"?
#' @param include_only only run the correlations that include the members (as a vector) or combinations (as a list or data.frame)
#' @param check_timing should we try to estimate run time for full dataset? (default is FALSE)
#' @param return_matrix should the data.frame or matrix result be returned?
#' 
#' @details For more details, see the ICI-Kendall-tau vignette 
#' 
#' \code{browseVignettes("ICIKendallTau")}
#' 
#'   The default for \code{global_na} includes what values in the data to replace with NA for the Kendall-tau calculation. By default these are \code{global_na = c(NA, Inf, 0)}. If you want to replace something other than 0, for example, you might use \code{global_na = c(NA, Inf, -2)}, and all values of -2 will be replaced instead of 0.
#'   
#'   When \code{check_timing = TRUE}, 5 random pairwise comparisons will be run to generate timings on a single core, and then estimates of how long the full set will take are calculated. The data is returned as a data.frame, and will be on the low side, but it should provide you with a good idea of how long your data will take.
#'   
#'   Returned is a list containing matrices with:
#'   
#'   * cor: scaled correlations
#'   * raw: raw kendall-tau correlations
#'   * pval: p-values
#'   * taumax: the theoretical maximum kendall-tau value possible
#'   
#'   Eventually, we plan to provide two more parameters for replacing values, \code{feature_na} for feature specific NA values and \code{sample_na} for sample specific NA values.
#' 
#' @return list with cor, raw, pval, taumax
#' 
#' @examples
#' \dontrun{
#' # not run
#' set.seed(1234)
#' s1 = sort(rnorm(1000, mean = 100, sd = 10))
#' s2 = s1 + 10 
#' 
#' matrix_1 = cbind(s1, s2)
#' 
#' r_1 = ici_kendalltau(t(matrix_1))
#' r_1$cor
#' 
#' #    s1 s2
#' # s1  1  1
#' # s2  1  1
#' names(r_1)
#' # "cor", "raw", "pval", "taumax", "keep", "run_time"
#' 
#' s3 = s1
#' s3[sample(100, 50)] = NA
#' 
#' s4 = s2
#' s4[sample(100, 50)] = NA
#' 
#' matrix_2 = cbind(s3, s4)
#' r_2 = ici_kendalltau(t(matrix_2))
#' r_2$cor
#' #           s3        s4
#' # s3 1.0000000 0.9944616
#' # s4 0.9944616 1.0000000
#' 
#' # using include_only
#' set.seed(1234)
#' x = matrix(rnorm(5000), nrow = 100, ncol = 50)
#' rownames(x) = paste0("s", seq(1, nrow(x)))
#' 
#' # only calculate correlations of other columns with "s1"
#' include_s1 = "s1"
#' s1_only = ici_kendalltau(x, include_only = include_s1)
#' 
#' # include s1 and s3 things both
#' include_s1s3 = c("s1", "s3")
#' s1s3_only = ici_kendalltau(x, include_only = include_s1s3)
#' 
#' # only specify certain pairs either as a list
#' include_pairs = list(g1 = "s1", g2 = c("s2", "s3"))
#' s1_other = ici_kendalltau(x, include_only = include_pairs)
#' 
#' # or a data.frame
#' include_df = as.data.frame(list(g1 = "s1", g2 = c("s2", "s3")))
#' s1_df = ici_kendalltau(x, include_only = include_df)
#' 
#' }
#' @export
#' 
ici_kendalltau = function(data_matrix, 
                             global_na = c(NA, Inf, 0),
                             perspective = "global",
                             scale_max = TRUE,
                             diag_good = TRUE,
                             include_only = NULL,
                             check_timing = FALSE,
                             return_matrix = TRUE){
  
  # assume row-wise (because that is what the description states), so need to transpose
  data_matrix = t(data_matrix)
  exclude_loc = matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  
  if (is.null(colnames(data_matrix))) {
    stop("rownames of data_matrix cannot be NULL!")
  }
  
  message("Processing missing values ...\n")
  
  # Actual NA and Inf values are special cases, so we do
  # this very specifically.
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
    
  }
  
  # now that we've done the NA and Inf values, we can go
  # ahead and take care of the rest.
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }

  # eventually we should be able to handle sample specific and 
  # feature specific values as well
  
  # after we've figured out what should be excluded, we actually set it to NA
  # and let the icikt code worry about it
  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  n_sample = ncol(exclude_data)
  
  # figure out if we are using furrr to process these
  if ("furrr" %in% utils::installed.packages()) {
    ncore = future::nbrOfWorkers()
    names(ncore) = NULL
    split_fun = furrr::future_map
  } else {
    ncore = 1
    split_fun = purrr::map
  }
  
  # generate the array of comparisons, 2 x ...,
  # where each column is a comparison between two columns of data
  message("Figuring out comparisons to do ...")
  pairwise_comparisons = utils::combn(n_sample, 2)
  
  if (!diag_good) {
    extra_comparisons = matrix(rep(seq(1, n_sample), each = 2), nrow = 2, ncol = n_sample, byrow = FALSE)
    pairwise_comparisons = cbind(pairwise_comparisons, extra_comparisons)
  }
  
  # create a data.frame of the comparisons by the names of the columns instead,
  # this enables indexed comparisons and named comparisons, because
  # we can have row / column names in R
  # This is now n_comparisons x 2
  named_comparisons = data.frame(s1 = colnames(data_matrix)[pairwise_comparisons[1, ]],
                                 s2 = colnames(data_matrix)[pairwise_comparisons[2, ]])
  
  if (!is.null(include_only)) {
    if (is.character(include_only) || is.numeric(include_only)) {
      #message("a vector!")
      # Check each of the comparison vectors against the include_only variable
      # This returns TRUE where they match
      # Use OR to make sure we return everything that should be returned
      s1_include = named_comparisons$s1 %in% include_only
      s2_include = named_comparisons$s2 %in% include_only
      named_comparisons = named_comparisons[(s1_include | s2_include), ]
    } else if (is.list(include_only)) {
      if (length(include_only) == 2) {
        #message("a list!")
        # In this case the include_only is a list of things, so we have to check both
        # of the sets against each of the lists. Again, this returns TRUE where
        # they match. Because we want the things where they
        # are both TRUE (assuming l1[1] goes with l2[1]), we use the AND at the end.
        l1_include = (named_comparisons$s1 %in% include_only[[1]]) | (named_comparisons$s2 %in% include_only[[1]])
        l2_include = (named_comparisons$s1 %in% include_only[[2]]) | (named_comparisons$s2 %in% include_only[[2]])
        named_comparisons = named_comparisons[(l1_include & l2_include), ]
      } else {
        stop("include_only must either be a single vector, or a list of 2 vectors!")
      }
    }
  }
  
  if (nrow(named_comparisons) == 0) {
    stop("nrow(named_comparisons) == 0, did you create include_only correctly?")
  }
  
  n_todo = nrow(named_comparisons)
  if (check_timing) {
    sample_compare = sample(nrow(named_comparisons), 5)
    tmp_pairwise = named_comparisons[, sample_compare]
    
    run_tmp = check_icikt_timing(exclude_data, tmp_pairwise, perspective, n_todo, ncore)
    return(run_tmp)
  }
  
  # in my experience, the list of actual comparisons to do does not get split
  # up well across compute cores (cores, machines, etc), and we end up waiting
  # on stuff to complete.
  # Therefore, this code actually does the splitting up of comparisons across
  # the number of cores (ncore) in a list. 
  
  message("Splitting up across compute ...")
  n_each = ceiling(n_todo / ncore)
  
  which_core = rep(seq(1, ncore), each = n_each)
  which_core = which_core[1:nrow(named_comparisons)]
  
  which_core = sample(which_core, length(which_core))
  
  named_comparisons$core = which_core
  named_comparisons$raw = Inf
  named_comparisons$pvalue = Inf
  named_comparisons$taumax = Inf
  
  
  split_comparisons = split(named_comparisons, named_comparisons$core)
  
  
  # finally, this is the function for doing each set of icikt comparisons,
  # possibly across cores.
  do_split = function(do_comparisons, exclude_data, perspective) {
    #seq_range = seq(in_range[1], in_range[2])
    #print(seq_range)
    
    raw = vector("numeric", nrow(do_comparisons))
    pvalue = raw
    taumax = raw
    
    for (irow in seq_len(nrow(do_comparisons))) {
      iloc = do_comparisons[irow, 1]
      jloc = do_comparisons[irow, 2]
      ici_res = ici_kt(exclude_data[, iloc], exclude_data[, jloc], perspective = perspective)
      raw[irow] = ici_res["tau"]
      pvalue[irow] = ici_res["pvalue"]
      taumax[irow] = ici_res["tau_max"]
    }
    do_comparisons$raw = raw
    do_comparisons$pvalue = pvalue
    do_comparisons$taumax = taumax
    #return(ls())
    do_comparisons
  }
  # we record how much time is actually spent doing ICI-Kt
  # itself, as some of the other operations will add a bit of time
  # 
  message("Running correlations ...")
  t1 = Sys.time()
  # note here, this takes our list of comparisons, and then calls the do_split
  # function above on each of them.
  split_cor = split_fun(split_comparisons, do_split, exclude_data, perspective)
  t2 = Sys.time()
  t_diff = as.numeric(difftime(t2, t1, units = "secs"))

  # put all the results back together again into one data.frame
  message("Recombining results ...")
  all_cor = purrr::list_rbind(split_cor)
  rownames(all_cor) = NULL
  
  # calculate the max-cor value for use in scaling across multiple comparisons
  # n_observations = nrow(exclude_data)
  # n_na = sort(colSums(exclude_loc))
  # m_value = floor(sum(n_na[1:2]) / 2)
  # n_m = n_observations - m_value
  # max_cor_denominator = choose(n_m, 2) + n_observations * m_value
  # max_cor_numerator = choose(n_m, 2) + n_observations * m_value + choose(m_value, 2)
  # max_cor = max_cor_denominator / max_cor_numerator
  
  if (scale_max) {
    max_cor = max(all_cor$taumax, na.rm = TRUE)
    all_cor$cor = all_cor$raw / max_cor
  } else {
    all_cor$cor = all_cor$raw
  }
  
  if (diag_good) {
    n_good = colSums(!exclude_loc)
    names(n_good) = colnames(data_matrix)
    extra_cor = data.frame(s1 = names(n_good),
                           s2 = names(n_good),
                           core = 0,
                           raw = n_good / max(n_good),
                           pvalue = 0,
                           taumax = 1,
                           cor = n_good / max(n_good))
    rownames(extra_cor) = NULL
    all_cor = rbind(all_cor, extra_cor)
  }
  
  # if the user asks for the matrix back, we give the matrices, otherwise we
  # leave them in the data.frame
  if (return_matrix) {
    message("Generating the output matrix ...")
    cor_matrix = matrix(0, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
    rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
    raw_matrix = cor_matrix
    pvalue_matrix = cor_matrix
    taumax_matrix = cor_matrix
    
    for (irow in seq_len(nrow(all_cor))) {
      rowloc = all_cor$s1[irow]
      colloc = all_cor$s2[irow]
      raw_matrix[rowloc, colloc] = raw_matrix[colloc, rowloc] = all_cor$raw[irow]
      cor_matrix[rowloc, colloc] = cor_matrix[colloc, rowloc] = all_cor$cor[irow]
      pvalue_matrix[rowloc, colloc] = pvalue_matrix[colloc, rowloc] = all_cor$pvalue[irow]
      taumax_matrix[rowloc, colloc] = taumax_matrix[colloc, rowloc] = all_cor$taumax[irow]
    }
    return(list(cor = cor_matrix, raw = raw_matrix, pval = pvalue_matrix, taumax = taumax_matrix, keep = t(!exclude_loc), run_time = t_diff))
  } else {
    return(list(cor = all_cor, run_time = t_diff))
  }
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
