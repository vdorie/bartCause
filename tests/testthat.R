if (require(testthat, quietly = TRUE)) {
  require(bartCause)
  test_check("bartCause")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}
