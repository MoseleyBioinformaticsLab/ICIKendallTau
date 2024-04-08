test_that("test_left_censorship works", {
  n_feature = 1000
  n_sample = 20
  n_miss = 100
  n_low = 80
  test_dataset = withr::with_seed(1234, rlnorm(n_feature, 10, 1))
  test_dataset = sort(test_dataset)
  
  noisy_dataset = withr::with_seed(1234, add_uniform_noise(n_sample, log(test_dataset), 0.1))
  sample_medians = calculate_matrix_medians(noisy_dataset)
  
  # check that we get the counts we expect
  # we are going to pick n_miss features, and of those, n_low will be less than 
  # the medians, hopefully by ordering.
  low_indices = withr::with_seed(1234, sample(seq_len(300), n_low))
  hi_indices = withr::with_seed(1234, sample(seq(800, 1000), n_miss - n_low))
  
  all_indices = c(low_indices, hi_indices)
  
  sample_indices = withr::with_seed(1234, sample(n_sample, n_miss, replace = TRUE))
  
  zero_dataset = noisy_dataset
  for (i_loc in seq_along(all_indices)) {
    zero_dataset[all_indices[i_loc], sample_indices[i_loc]] = 0
  }
  
  zero_binom = test_left_censorship(zero_dataset, global_na = c(0, NA))
  
  expect_equal(zero_binom$values$success / zero_binom$values$trials, n_low / n_miss)
  
  na_dataset = zero_dataset
  na_dataset[zero_dataset == 0] = NA
  
  na_binom = test_left_censorship(na_dataset, global_na = c(0, NA))
  expect_equal(na_binom$values$success / na_binom$values$trials, n_low / n_miss)
  
  # check that zero and NA work the same
  expect_equal(zero_binom, na_binom)
  
  na_loc = which(is.na(na_dataset))
  back_zero = withr::with_seed(1234, sample(na_loc, 20))
  mix_dataset = na_dataset
  mix_dataset[back_zero] = 0
  mix_binom = test_left_censorship(mix_dataset, global_na = c(0, NA))
  expect_equal(na_binom, mix_binom)
})
