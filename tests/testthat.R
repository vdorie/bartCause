if (require(testthat, quietly = TRUE)) {
  require(BartCause)
  test_check("BartCause")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}
