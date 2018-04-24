context("Mantel test")
test_that("Permutation test and plotting", {
  expect_silent({
  	result <- mantel.test(list(dist(1:8), dist(sample(8:1)), dist(runif(8))),
    hammingdists(enumerate.meaningcombinations(c(2, 2, 2))))
    plot(result)
    plot(result[1,])
  })
})
