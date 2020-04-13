test_that("class plot", {

  expect_equal(
    class(
      create_histogram(
        aupair::test_data,
        counts
      ))[2],
     "ggplot"
  )


})
