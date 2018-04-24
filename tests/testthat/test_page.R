context("Page test")
test_that("Exact p-value calculation", {
  expect_equal(1, page.compute.exact(6, 4, 224))
  expect_equal(0.03932265, page.compute.exact(6, 4, 322))
  expect_silent(page.test(rbind(1:10, 1:10), verbose=FALSE))
})
test_that("Approximate p-value calculation and border conditions", {
  expect_error(page.test(t(1:12)))
  expect_error(page.compute.exact(6, 4))
  expect_error(page.compute.exact(6, 4, 223))
  expect_error(page.compute.exact(6, 4, 365))
  expect_message(page.test(rbind(1:23, 1:23)))
})

