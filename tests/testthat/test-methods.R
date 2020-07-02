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

# test_that("gate dimensions", {
# 
#   res =
#     reduce_dimensions(nanny::test_data2 , c(`ct 1`, `ct 2`), cancer_ID, relation, method="MDS") %>%
#     gate_dimensions(.element = c(`ct 1`, `ct 2`), Dim1, Dim2 )
#   expect_equal(ncol(res) , 10)
# 
# })

test_that("remove redundancy", {
  res = remove_redundancy(nanny::test_data2 , cancer_ID, c(`ct 1`, `ct 2`), relation) 
  expect_equal(ncol(res) , 6)
  
  res = remove_redundancy(nanny::test_data2 ,  c(`ct 1`, `ct 2`), cancer_ID, relation) 
  expect_equal(ncol(res) , 6)
  
})

test_that("subset", {
  res = subset(nanny::test_data2 , cancer_ID) 
  expect_equal(nrow(res) , 32)
  
  res = subset(nanny::test_data2 ,  c(`ct 1`, `ct 2`)) 
  expect_equal(nrow(res) , 70)
  
})

test_that("impute missing", {
  res = impute_missing(nanny::test_data2 , cancer_ID, c(`ct 1`, `ct 2`), relation, ~ 1) 
  expect_identical(res , nanny::test_data2 %>% mutate_if(is.factor, as.character))
  
  res = impute_missing(nanny::test_data2 ,  c(`ct 1`, `ct 2`),cancer_ID, relation, ~ 1)
  expect_identical(res , nanny::test_data2 )
  
  res = impute_missing(nanny::test_data2 %>% slice(-1),   c(`ct 1`, `ct 2`),cancer_ID, relation, ~ 1) 
  expect_equal(res %>% inner_join(nanny::test_data2 %>% slice(1) %>% select(-relation, -group)) %>% nrow , 1)
  
  # Test with formula and covariate
  res = impute_missing(
    nanny::test_data2 %>% 
      slice(-1) %>% 
      left_join( 
        (.) %>%
          distinct(`ct 1`, `ct 2`) %>%
          mutate(cov = sample(0:1, size = n(), replace = TRUE, prob = c(0.1, 0.9)) %>% as.factor)
      ),  c(`ct 1`, `ct 2`),cancer_ID, relation, ~ cov) 
  expect_equal(res %>% inner_join(nanny::test_data2 %>% slice(1) %>% select(-relation, -group)) %>% nrow , 1)
  
  
})

test_that("fill missing", {
  res = fill_missing(nanny::test_data2 ,cancer_ID, c(`ct 1`, `ct 2`), relation, fill_with = 0) 
  expect_identical(res , nanny::test_data2 %>% mutate_if(is.factor, as.character))
  
  res = fill_missing(nanny::test_data2 ,   c(`ct 1`, `ct 2`),cancer_ID, relation,  fill_with = 0) 
  expect_identical(res , nanny::test_data2 )
  
  res = fill_missing(nanny::test_data2 %>% slice(-1),  c(`ct 1`, `ct 2`),cancer_ID, relation,    fill_with = 0) 
  expect_equal(res %>% inner_join(nanny::test_data2 %>% slice(1) %>% select(-relation, -group)) %>% pull(relation) , 0)
  
})

test_that("permute nest", {
  res = permute_nest(nanny::test_data2 ,cancer_ID, relation) 
  expect_equal(nrow(res) ,992)
  
  res = permute_nest(nanny::test_data2 ,   `ct 1`, relation) 
  expect_equal(nrow(res) ,992)
  
  res = permute_nest(nanny::test_data2 ,   `ct 1`, c(group, cancer_ID)) 
  expect_equal(nrow(res) ,992)
  
})

test_that("combine nest", {
  res = combine_nest(nanny::test_data2 ,cancer_ID, relation) 
  expect_equal(nrow(res) ,496)
  
  res = combine_nest(nanny::test_data2 ,   `ct 1`, relation) 
  expect_equal(nrow(res) ,496)
  
  res = combine_nest(nanny::test_data2 ,   `ct 1`, c(group, cancer_ID)) 
  expect_equal(nrow(res) ,496)
  
})


test_that("keep variable", {
  res = keep_variable(nanny::test_data2 ,cancer_ID, c(`ct 1`, `ct 2`), relation, top = 10) 
  expect_equal(nrow(res) ,320)
  
  res = keep_variable(nanny::test_data2 ,   c(`ct 1`, `ct 2`),cancer_ID, relation, top=10) 
  expect_equal(nrow(res) ,700)
  
})

test_that("lower triangular", {
  res = lower_triangular(nanny::test_data2 %>% filter(cancer_ID == "ACC") ,`ct 2`, `ct 1`, relation) 
  expect_equal(nrow(res) ,35)
  
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

test_that("subset", {
  res = nest_subset(mtcars_tidy,data = -car_model)
  expect_equal(ncol(res) , 4)
  
  res = nest_subset(nanny::test_data2 ,  data = -c(`ct 1`, `ct 2`)) 
  expect_equal(ncol(res) , 4)
  
})

