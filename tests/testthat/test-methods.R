test_that("cluster kmeans", {
  res = cluster_elements(nanny::test_data2 , cancer_ID, c(`ct 1`, `ct 2`), relation, centers = 2, method="kmeans")
  expect_equal(ncol(res) , 7)

  res = cluster_elements(nanny::test_data2 , c(`ct 1`, `ct 2`), cancer_ID, relation, centers = 2, method="kmeans")
  expect_equal(ncol(res) , 7)

})

test_that("reduce dimension PCA", {
  res = reduce_dimensions(nanny::test_data2 , cancer_ID, c(`ct 1`, `ct 2`), relation, method="PCA")
  expect_equal(ncol(res) , 8)
  
  res = reduce_dimensions(nanny::test_data2 , c(`ct 1`, `ct 2`), cancer_ID, relation, method="PCA")
  expect_equal(ncol(res) , 8)
  
})

test_that("reduce dimension MDS", {
  res = reduce_dimensions(nanny::test_data2 , cancer_ID, c(`ct 1`, `ct 2`), relation, method="MDS")
  expect_equal(ncol(res) , 8)
  
  res = reduce_dimensions(nanny::test_data2 , c(`ct 1`, `ct 2`), cancer_ID, relation, method="MDS")
  expect_equal(ncol(res) , 8)
  
})

test_that("as matrix", {
  tibble(a = 1:10, b = 1:10) %>%
    as_matrix()
  
})