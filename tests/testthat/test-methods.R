test_that("class plot", {

  expect_equal(
    class(
      create_histogram(
        tibble::tibble(counts = 1:100),
        counts
      ))[2],
     "ggplot"
  )


})
