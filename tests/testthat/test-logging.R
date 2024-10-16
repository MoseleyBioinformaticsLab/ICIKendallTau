test_that("logger detection works", {
  new_lib = tempfile()
  dir.create(new_lib)
  withr::with_libpaths(new_lib, {
    expect_error(enable_logging(), "Logging requested")
  })
  unlink(new_lib, recursive = TRUE)
})

test_that("logging works", {
  skip_on_os(c("windows", "mac"))
  enable_logging("nomemory_logs.log")
  
  set.seed(1234)
  x = matrix(rnorm(5000), nrow = 50, ncol = 100)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  ici_res = ici_kendalltau(x)
  
  log_contents = readLines("nomemory_logs.log")
  
  expect_false(all(grepl("Memory", log_contents)))
  
  enable_logging("withmemory_logs.log", memory = TRUE)
  ici_res = ici_kendalltau(x)
  log_contents = readLines("withmemory_logs.log")  
  expect_true(any(grepl("Memory", log_contents)))
  unlink("nomemory_logs.log", recursive = TRUE)
  unlink("withmemory_logs.log", recursive = TRUE)
})

test_that("progress works", {
  show_progress()
  
  set.seed(1234)
  x = matrix(rnorm(5000), nrow = 50, ncol = 100)
  colnames(x) = paste0("s", seq(1, ncol(x)))
  expect_message(ici_kendalltau(x), "Processing missing values")
})
