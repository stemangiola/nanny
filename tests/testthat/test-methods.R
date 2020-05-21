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

test_that("rotate dimensions", {
  res = reduce_dimensions(nanny::test_data2 , cancer_ID, c(`ct 1`, `ct 2`), relation, method="MDS") %>%
    rotate_dimensions(dimension_1_column = Dim1, dimension_2_column = Dim2, rotation_degrees = 45, .element = cancer_ID)
  expect_equal(ncol(res) , 10)
  
  res = reduce_dimensions(nanny::test_data2 , c(`ct 1`, `ct 2`), cancer_ID, relation, method="MDS") %>%
    rotate_dimensions(dimension_1_column = Dim1, dimension_2_column = Dim2, rotation_degrees = 45, .element = c(`ct 1`, `ct 2`))
  expect_equal(ncol(res) , 10)
  
})

test_that("as matrix", {
  
  res = tibble(a = 1:10, b = 1:10) %>% as_matrix()
  expect_equal(ncol(res) , 2)
  
  res = 
    nanny::test_data2 %>% select(`ct 1` ,     `ct 2`  ,    relation, cancer_ID) %>% 
    spread(cancer_ID, relation) %>%
    as_matrix(rownames = c(`ct 1` , `ct 2` ))
  expect_equal(dim(res) , c(70, 32))
  

})
