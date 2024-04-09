#' Test for left censorship
#' 
#' Does a binomial test to check if the most likely cause of missing values
#' is due to values being below the limit of detection, or coming from a 
#' left-censored distribution.
#' 
#' @param in_data matrix or data.frame of numeric data
#' @param sample_classes which samples are in which class
#' @param global_na what represents zero or missing?
#' 
#' @details
#' For each feature that is missing in a group of samples, we save as a possibility
#' to test. For each sample, we calculate the median value with any missing values
#' removed. Each feature that had a missing value, we test whether the remaining 
#' non-missing values are below the sample median for those samples where the 
#' feature is non-missing. A binomial test considers the total number of features
#' instances (minus missing values) as the number of trials, and the number of
#' of features below the sample medians as the number of successes.
#' 
#' @seealso [vignette("testing-for-left-censorship", package = "ICIKendallTau")]
#' 
#' @examples
#' # this example has 80% missing due to left-censorship
#' data(missing_dataset)
#' missingness = test_left_censorship(missing_dataset)
#' missingness$values
#' missingness$binomial_test
#' 
#' @export
#' @return data.frame of trials / successes, and binom.test result
test_left_censorship = function(in_data, 
                                sample_classes = NULL, 
                                global_na = c(NA, Inf, 0))
{
  if (is.null(sample_classes)) {
    sample_classes = rep("A", ncol(in_data))
  }
  
  split_indices = split(seq_len(ncol(in_data)), sample_classes)
  missing_loc = setup_missing_matrix(in_data, global_na)
  in_data_missing = in_data
  in_data_missing[missing_loc] = NA
  
  # split the dataset by group
  split_counts = purrr::imap(split_indices, \(in_split, split_id){
    # in_split = split_indices[[1]]
    
    # grab the group we want to work with
    split_missing = in_data_missing[, in_split, drop = FALSE]
    
    # count the number of missing samples for each feature,
    # and keep those that have at least one
    n_miss = rowSums(is.na(split_missing))
    keep_miss = split_missing[n_miss > 0, ]
    
    # get sample medians
    sample_medians = calculate_matrix_medians(split_missing, use = "col", na.rm = TRUE)
    
    # turn the medians into a matrix to make life easier
    median_matrix = matrix(sample_medians, nrow = nrow(keep_miss),
                           ncol = ncol(keep_miss), byrow = TRUE)
    # do the comparison
    keep_miss_updown = keep_miss < median_matrix
    
    # count how many trials we ran, and how many successes we have
    all_trials = (nrow(keep_miss_updown) * ncol(keep_miss_updown)) - sum(is.na(keep_miss_updown))
    all_success = sum(keep_miss_updown, na.rm = TRUE)
    
    data.frame(trials = all_trials, success = all_success, class = split_id)
  }) |>
    purrr::list_rbind()
  
  total_trials = sum(split_counts$trials)
  total_success = sum(split_counts$success)
  
  binom_res = stats::binom.test(total_success, total_trials, p = 0.5, alternative = "greater")
  return(list(values = split_counts,
              binomial_test = binom_res))
}


#' Calculate matrix medians
#' 
#' Given a matrix of data, calculates the median value in each column or row.
#' 
#' @param in_matrix numeric matrix of values
#' @param use character of "col" or "row" defining columns or rows
#' @param ... extra parameters to the median function
#' 
#' @export
#' @return numeric
calculate_matrix_medians = function(in_matrix, use = "col", ...)
{
  if (use %in% "row") {
    in_matrix = t(in_matrix)
  } 
  out_medians = purrr::map_dbl(seq_len(ncol(in_matrix)), \(in_col){
    stats::median(in_matrix[, in_col], ...)
  })
  return(out_medians)
}

#' Add uniform noise
#' 
#' Adds uniform noise to values, generating replicates with noise added to
#' the original.
#' 
#' @param value a single or vector of numeric values
#' @param n_rep the number of replicates to make (numeric). Default is 1.
#' @param sd the standard deviation of the data
#' @param use_zero logical, should returned values be around zero or not?
#' 
#' @export
#' 
#' @return numeric matrix
add_uniform_noise = function(value, n_rep, sd, use_zero = FALSE){
  n_value = length(value)
  
  n_sd = n_rep * n_value
  
  out_sd = stats::rnorm(n_sd, 0, sd)
  out_sd = matrix(out_sd, nrow = n_value, ncol = n_rep)
  
  if (!use_zero) {
    tmp_value = matrix(value, nrow = n_value, ncol = n_rep, byrow = FALSE)
    out_value = tmp_value + out_sd
  } else {
    out_value = out_sd
  }
  
  return(out_value)
}

#' Example Dataset With Missingness
#'
#' An example dataset that has missingness from left-censorship
#'
#' @format ## `missing_dataset`
#' A matrix with 1000 rows and 20 columns, where rows are features and
#' columns are samples.
#' 
#' @source Robert M Flight
"missing_dataset"