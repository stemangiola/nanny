test_that("class plot", {

  expect_equal(
    class(
      create_histogram(
        nanny::test_data,
        counts
      ))[2],
     "ggplot"
  )


})
