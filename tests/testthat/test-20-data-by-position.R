context("data given by position")

test_that("data given by position is refused with the name to pass it by", {
  expect_error(bartc(y, z, x1 + x2, formData), "the fourth argument of bartc is 'parametric', which was given a data frame; 'data' must be passed by name, data = formData", fixed = TRUE)
  expect_error(bartc(y, z, x1 + x2, formData, group.by = country),
               "the fourth argument of bartc is 'parametric'", fixed = TRUE)
  expect_error(bartc(y, z, x1 + x2, formData[formData$x1 > 0, ]),
               "data = formData[formData$x1 > 0, ]", fixed = TRUE)
  
  ## by name it is data
  expect_true(inherits(fit(), "bartcFit"))
  
  ## a data frame is a parametric component of its own where the call is made of vectors
  x <- as.matrix(formData[c("x1", "x2")])
  w <- data.frame(w = formData$x2)
  message <- tryCatch({
    bartc(formData$y, formData$z, x, w, n.samples = 10L, n.burn = 10L, n.chains = 1L, n.threads = 1L,
          n.trees = 10L, verbose = FALSE)
    ""
  }, error = function(e) conditionMessage(e))
  expect_false(grepl("fourth argument", message, fixed = TRUE))
})
