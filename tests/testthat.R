if (require(testthat, quietly = TRUE)) {
  require(cibart)
  test_check("cibart")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}
