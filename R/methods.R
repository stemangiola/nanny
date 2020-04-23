#' Get clusters of elements (e.g., elements or features)
#'
#' \lifecycle{maturing}
#'
#' @description cluster_elements() takes as imput a `tbl` formatted as | <element> | <feature> | <value> | <...> | and identify clusters in the data.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name cluster_elements
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column (normally elements).
#' @param .feature The name of the feature column (normally features)
#' @param .value The name of the column including the numerical value the clustering is based on (normally feature value)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_elements A boolean. In case the input is a nanny object, it indicates Whether the element column will be element or feature column
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans
#' 
#' @details identifies clusters in the data, normally of elements.
#' This function returns a tibble with additional columns for the cluster annotation.
#' At the moment only k-means clustering is supported, the plan is to introduce more clustering methods.
#'
#' @return A tbl object with additional columns with cluster labels
#'
#'
#' @examples
#'
#'
#'     cluster_elements(nanny::counts_mini, element, feature, count,	centers = 2, method="kmeans")
#'
#' @docType methods
#' @rdname cluster_elements-methods
#' @export
#'
setGeneric("cluster_elements", function(.data,
																				.element = NULL,
																				.feature = NULL,
																				.value = NULL,
																				method,
																				of_elements = TRUE,
																				transform = NULL,
																				action = "add",
																				...)
	standardGeneric("cluster_elements"))

# Set internal
.cluster_elements = 		function(.data,
															 .element = NULL,
															 .feature = NULL,
															 .value = NULL,
															 method ,
															 of_elements = TRUE,
															 transform = NULL,
															 action = "add",
															 ...)
{
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Validate data frame
	validation(.data, !!.element, !!.feature, !!.value)
	
	# # Check if data rectangular
	# ifelse_pipe(
	# 	(.) %>% check_if_data_rectangular(!!.element,!!.feature,!!.value, type = "soft"),
	# 	~ .x %>% eliminate_sparse_features(!!.feature)
	# ) %>%
	
	if (method == "kmeans") {
		if (action == "add"){
			
			.data %>%
				dplyr::left_join(
					(.) %>%
						get_clusters_kmeans_bulk(
							.value = !!.value,
							.element = !!.element,
							.feature = !!.feature,
							of_elements = of_elements,
							transform = transform,
							...
						)
				) 
		}
		else if (action == "get"){
			
			.data %>%
				
				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.value)$horizontal_cols
				) %>%
				distinct() %>%
				
				dplyr::left_join(
					.data %>%
						get_clusters_kmeans_bulk(
							.value = !!.value,
							.element = !!.element,
							.feature = !!.feature,
							of_elements = of_elements,
							transform = transform,
							...
						)
				) 
			
		}
		else if (action == "only")
			get_clusters_kmeans_bulk(
				.data,
				.value = !!.value,
				.element = !!.element,
				.feature = !!.feature,
				of_elements = of_elements,
				transform = transform,
				...
			)
		else
			stop(
				"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}
	else if (method == "SNN") {
		if (action == "add"){
			
			.data %>%
				dplyr::left_join(
					(.) %>%
						get_clusters_SNN_bulk(
							.value = !!.value,
							.element = !!.element,
							.feature = !!.feature,
							of_elements = of_elements,
							transform = transform,
							...
						)
				) 
			
		}
		else if (action == "get"){
			
			.data %>%
				
				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.value)$horizontal_cols
				) %>%
				distinct() %>%
				
				dplyr::left_join(
					.data %>%
						get_clusters_SNN_bulk(
							.value = !!.value,
							.element = !!.element,
							.feature = !!.feature,
							of_elements = of_elements,
							transform = transform,
							...
						)
				)
			
		}
		
		else if (action == "only")
			get_clusters_SNN_bulk(
				.data,
				.value = !!.value,
				.element = !!.element,
				.feature = !!.feature,
				of_elements = of_elements,
				transform = transform,
				...
			)
		else
			stop(
				"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}
	else
		stop("nanny says: the only supported methods are \"kmeans\" or \"SNN\" ")
	
}

#' cluster_elements
#' @inheritParams cluster_elements
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "spec_tbl_df", .cluster_elements)

#' cluster_elements
#' @inheritParams cluster_elements
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "tbl_df", .cluster_elements)

#' Dimension reduction of the feature value data
#'
#' \lifecycle{maturing}
#'
#' @description reduce_dimensions() takes as imput a `tbl` formatted as | <element> | <feature> | <value> | <...> | and calculates the reduced dimensional space of the feature value.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name reduce_dimensions
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column (normally elements).
#' @param .feature The name of the feature column (normally features)
#' @param .value The name of the column including the numerical value the clustering is based on (normally feature value)
#'
#' @param method A character string. The dimension reduction algorithm to use (PCA, MDS, tSNE).
#' @param top An integer. How many top genes to select for dimensionality reduction
#' @param of_elements A boolean. In case the input is a nanny object, it indicates Whether the element column will be element or feature column
#' @param .dims A list of integer vectors corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function prcomp if you choose method="PCA" or Rtsne if you choose method="tSNE"
#'
#' @details This function reduces the dimensions of the feature values.
#' It can use multi-dimensional scaling (MDS) of principal component analysis (PCA).
#'
#' @return A tbl object with additional columns for the reduced dimensions
#'
#'
#' @examples
#'
#'
#'
#' counts.MDS =  reduce_dimensions(nanny::counts_mini, element, feature, count, method="MDS", .dims = 3)
#'
#' counts.PCA =  reduce_dimensions(nanny::counts_mini, element, feature, count, method="PCA", .dims = 3)
#'
#'
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @export
#'
#'
setGeneric("reduce_dimensions", function(.data,
																				 .element = NULL,
																				 .feature = NULL,
																				 .value = NULL,
																				 method,
																				 .dims = 2,
																				 top = Inf,
																				 of_elements = TRUE,
																				 transform = NULL,
																				 scale = TRUE,
																				 action = "add",
																				 ...)
					 standardGeneric("reduce_dimensions"))

# Set internal
.reduce_dimensions = 		function(.data,
																.element = NULL,
																.feature = NULL,
																.value = NULL,
																method,
																.dims = 2,
																top = Inf,
																of_elements = TRUE,
																transform = NULL,
																scale = TRUE,
																action = "add",
																...)
{
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)

	# Validate data frame
	validation(.data, !!.element, !!.feature, !!.value)
	
	if (method == "MDS") {
		
		.data_processed =
			.data %>%
			get_reduced_dimensions_MDS_bulk(
				.value = !!.value,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_elements = of_elements,
				transform = transform,
				...
			)
		
		if (action == "add"){
			
			.data %>%	dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%
				reattach_internals(.data_processed)
			
		}
		else if (action == "get"){
			
			.data %>%
				
				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.value)$horizontal_cols
				) %>%
				distinct() %>%
				
				dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%
				reattach_internals(.data_processed)
			
		}
		
		else if (action == "only") .data_processed
		else
			stop(
				"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
	}
	else if (method == "PCA") {
		
		.data_processed =
			.data %>%
			get_reduced_dimensions_PCA_bulk(
				.value = !!.value,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_elements = of_elements,
				transform = transform,
				scale = scale,
				...
			)
		
		if (action == "add"){
			
			.data %>%
				dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%
				reattach_internals(.data_processed)
			
		}
		
		else if (action == "get"){
			
			.data %>%
				
				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.value)$horizontal_cols
				) %>%
				distinct() %>%
				
				dplyr::left_join(.data_processed,	by = quo_name(.element)) %>%
				reattach_internals(.data_processed)
			
		}
		
		else if (action == "only")	.data_processed
		else
			stop(
				"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
		
	}
	else if (method == "tSNE") {
		
		.data_processed =
			.data %>%
			get_reduced_dimensions_TSNE_bulk(
				.value = !!.value,
				.dims = .dims,
				.element = !!.element,
				.feature = !!.feature,
				top = top,
				of_elements = of_elements,
				transform = transform,
				...
			)
		
		if (action == "add"){
			
			.data %>%
				dplyr::left_join(.data_processed,	by = quo_name(.element)	) %>%
				reattach_internals(.data_processed)
			
		}
		else if (action == "get"){
			
			.data %>%
				
				# Selecting the right columns
				select(
					!!.element,
					get_x_y_annotation_columns(.data, !!.element,!!.feature, !!.value)$horizontal_cols
				) %>%
				distinct() %>%
				
				dplyr::left_join(.data_processed,	by = quo_name(.element)	) %>%
				reattach_internals(.data_processed)
			
		}
		else if (action == "only") .data_processed
		else
			stop(
				"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
			)
		
	}
	else
		stop("nanny says: method must be either \"MDS\" or \"PCA\"")
	
}

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "spec_tbl_df", .reduce_dimensions)

#' reduce_dimensions
#' @inheritParams reduce_dimensions
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "tbl_df", .reduce_dimensions)

#' Rotate two dimensions (e.g., principal components) of an arbitrary angle
#'
#' \lifecycle{maturing}
#'
#' @description rotate_dimensions() takes as imput a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the feature value.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name rotate_dimensions
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column (normally elements).
#'
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column  A character string. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param of_elements A boolean. In case the input is a nanny object, it indicates Whether the element column will be element or feature column
#' @param dimension_1_column_rotated A character string. The column of the rotated dimension 1 (optional)
#' @param dimension_2_column_rotated A character string. The column of the rotated dimension 2 (optional)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This function to rotate two dimensions such as the reduced dimensions.
#'
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
#'
#'
#' @examples
#'
#' counts.MDS =  reduce_dimensions(nanny::counts_mini, element, feature, count, method="MDS", .dims = 3)
#'
#' counts.MDS.rotated =  rotate_dimensions(counts.MDS, `Dim1`, `Dim2`, rotation_degrees = 45, .element = element)
#'
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#' @export
#'
setGeneric("rotate_dimensions", function(.data,
																				 dimension_1_column,
																				 dimension_2_column,
																				 rotation_degrees,
																				 .element = NULL,
																				 of_elements = TRUE,
																				 dimension_1_column_rotated = NULL,
																				 dimension_2_column_rotated = NULL,
																				 action = "add")
	standardGeneric("rotate_dimensions"))

# Set internal
.rotate_dimensions = 		function(.data,
																dimension_1_column,
																dimension_2_column,
																rotation_degrees,
																.element = NULL,
																of_elements = TRUE,
																dimension_1_column_rotated = NULL,
																dimension_2_column_rotated = NULL,
																action =	"add")
{
	# Get column names
	.element = enquo(.element)

	# Parse other colnames
	dimension_1_column = enquo(dimension_1_column)
	dimension_2_column = enquo(dimension_2_column)
	dimension_1_column_rotated = enquo(dimension_1_column_rotated)
	dimension_2_column_rotated = enquo(dimension_2_column_rotated)
	
	# Set default col names for rotated dimensions if not set
	if (quo_is_null(dimension_1_column_rotated))
		dimension_1_column_rotated = as.symbol(sprintf(
			"%s rotated %s",
			quo_name(dimension_1_column),
			rotation_degrees
		))
	if (quo_is_null(dimension_2_column_rotated))
		dimension_2_column_rotated = as.symbol(sprintf(
			"%s rotated %s",
			quo_name(dimension_2_column),
			rotation_degrees
		))
	
	.data_processed =
		get_rotated_dimensions(
			.data,
			dimension_1_column = !!dimension_1_column,
			dimension_2_column = !!dimension_2_column,
			rotation_degrees = rotation_degrees,
			.element = !!.element,
			of_elements = of_elements,
			dimension_1_column_rotated = !!dimension_1_column_rotated,
			dimension_2_column_rotated = !!dimension_2_column_rotated
		)
	
	if (action == "add"){
		
		.data %>%
			dplyr::left_join(	.data_processed,	by = quo_name(.element)	) 
		
	}
	else if (action == "get"){
		
		.data %>%
			
			# Selecting the right columns
			select(
				!!.element,
				get_specific_annotation_columns(.data, !!.element)
			) %>%
			distinct() %>%
			
			dplyr::left_join(	.data_processed,	by = quo_name(.element)	) 
		
	}
	else if (action == "only") .data_processed
	else
		stop(
			"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "spec_tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @inheritParams rotate_dimensions
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "tbl_df", .rotate_dimensions)

#' Drop redundant elements (e.g., elements) for which feature (e.g., feature/gene) aboundances are correlated
#'
#' \lifecycle{maturing}
#'
#' @description remove_redundancy() takes as imput a `tbl` formatted as | <element> | <feature> | <value> | <...> | for correlation method or | <DIMENSION 1> | <DIMENSION 2> | <...> | for reduced_dimensions method, and returns a `tbl` with dropped elements (e.g., elements).
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name remove_redundancy
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column (normally elements).
#' @param .feature The name of the feature column (normally features)
#' @param .value The name of the column including the numerical value the clustering is based on (normally feature value)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_elements A boolean. In case the input is a nanny object, it indicates Whether the element column will be element or feature column
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param top An integer. How many top genes to select for correlation based method
#' @param Dim_a_column A character string. For reduced_dimension based calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based calculation. The column of another principal component
#'
#'
#' @details This function removes redundant elements from the original data set (e.g., elements or features). For example, if we want to define cell-type specific signatures with low element redundancy. This function returns a tibble with dropped recundant elements (e.g., elements). Two redundancy estimation approaches are supported: (i) removal of highly correlated clusters of elements (keeping a representative) with method="correlation"; (ii) removal of most proximal element pairs in a reduced dimensional space.
#'
#' @return A tbl object with with dropped recundant elements (e.g., elements).
#'
#' @examples
#'
#'
#'
#'    remove_redundancy(
#'     nanny::counts_mini,
#' 	   .element = element,
#' 	   .feature = feature,
#' 	   	.value =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' counts.MDS =  reduce_dimensions(nanny::counts_mini, element, feature, count, method="MDS", .dims = 3)
#'
#' remove_redundancy(
#' 	counts.MDS,
#' 	Dim_a_column = `Dim1`,
#' 	Dim_b_column = `Dim2`,
#' 	.element = element,
#'   method = "reduced_dimensions"
#' )
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#' @export
#'
#'
setGeneric("remove_redundancy", function(.data,
																				 .element = NULL,
																				 .feature = NULL,
																				 .value = NULL,
																				 method,
																				 of_elements = TRUE,
																				 correlation_threshold = 0.9,
																				 top = Inf,
																				 transform = NULL,
																				 
																				 Dim_a_column,
																				 Dim_b_column)
					 standardGeneric("remove_redundancy"))

# Set internal
.remove_redundancy = 	 function(.data,
																.element = NULL,
																.feature = NULL,
																.value = NULL,
																method,
																of_elements = TRUE,
																correlation_threshold = 0.9,
																top = Inf,
																transform = NULL,
																Dim_a_column = NULL,
																Dim_b_column = NULL)
{
	# Make col names
	.value = enquo(.value)
	.element = enquo(.element)
	.feature = enquo(.feature)
	
	Dim_a_column = enquo(Dim_a_column)
	Dim_b_column = enquo(Dim_b_column)
	
	if (method == "correlation") {
		# Validate data frame
		validation(.data, !!.element, !!.feature, !!.value)
		
		remove_redundancy_elements_through_correlation(
			.data,
			.value = !!.value,
			.element = !!.element,
			.feature = !!.feature,
			correlation_threshold = correlation_threshold,
			top = top,
			of_elements = of_elements,
			transform = transform
		)
	}
	else if (method == "reduced_dimensions") {
		# Validate data frame
		# MISSING because feature not needed. I should build a custom funtion.
		
		remove_redundancy_elements_though_reduced_dimensions(
			.data,
			Dim_a_column = !!Dim_a_column,
			Dim_b_column = !!Dim_b_column,
			.element = !!.element,
			of_elements = of_elements
		)
	}
	else
		stop(
			"nanny says: method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
		)
	
}

#' remove_redundancy
#' @inheritParams remove_redundancy
#' @return A tbl object with with dropped recundant elements (e.g., elements).
setMethod("remove_redundancy", "spec_tbl_df", .remove_redundancy)

#' remove_redundancy
#' @inheritParams remove_redundancy
#' @return A tbl object with with dropped recundant elements (e.g., elements).
setMethod("remove_redundancy", "tbl_df", .remove_redundancy)

#' Extract sampe-wise information
#'
#' \lifecycle{maturing}
#'
#' @description subset() takes as imput a `tbl` formatted as | <element> | <ENSEMBL_ID> | <value> | <...> | and returns a `tbl` with only sampe-related columns
#'
#' @importFrom magrittr "%>%"
#'
#' @name subset
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column
#'
#'
#' @details This functon extracts only element-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to nanny function.
#'
#' @return A `tbl` object
#'
#'
#'
#'
#' @examples
#'
#'
#' 	subset(
#'			nanny::counts_mini,
#'			.element = element
#'		)
#'
#'
#' @docType methods
#' @rdname subset-methods
#' @export
#'
#'
setGeneric("subset", function(.data,
																		.column)
	standardGeneric("subset"))

# Set internal
.subset = 		function(.data,
										 .column)	{
	# Make col names
	.column = enquo(.column)
	
	.data %>%
		
		# Selecting the right columns
		select(	!!.column,	get_specific_annotation_columns(.data, !!.column)	) %>%
		distinct()
	
}

#' subset
#' @inheritParams subset
#' @return A `tbl` object
setMethod("subset",
					"spec_tbl_df",
					.subset)

#' subset
#' @inheritParams subset
#' @return A `tbl` object
setMethod("subset",
					"tbl_df",
					.subset)


#' Impute feature value if missing from element-feature pairs
#'
#' \lifecycle{maturing}
#'
#' @description impute_missing() takes as imput a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a `tbl` with an edditional adjusted value column. This method uses scaled counts if present.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name impute_missing
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_intrest + batch)
#' @param .element The name of the element column
#' @param .feature The name of the feature/gene column
#' @param .value The name of the feature/gene value column
#'
#' @details This function imputes the value of missing element-feature pair using the median of the element group defined by the formula
#'
#' @return A `tbl` non-sparse value
#'
#'
#'
#'
#' @examples
#'
#'
#' res =
#' 	impute_missing(
#' 		nanny::counts_mini,
#' 	~ condition,
#' 	.element = element,
#' 	.feature = feature,
#' 	.value = count
#' )
#'
#'
#' @docType methods
#' @rdname impute_missing-methods
#'
#' @export
#'
#'
setGeneric("impute_missing", function(.data,
																				.formula,
																				.element = NULL,
																				.feature = NULL,
																				.value = NULL)
	standardGeneric("impute_missing"))

# Set internal
.impute_missing = 	function(.data,
															.formula,
															.element = NULL,
															.feature = NULL,
															.value = NULL)
{
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Validate data frame
	validation(.data, !!.element, !!.feature, !!.value)
	
	fill_NA_using_formula(
		.data,
		.formula,
		.element = !!.element,
		.feature = !!.feature,
		.value = !!.value)
}

#' impute_missing
#' @inheritParams impute_missing
#' @return A `tbl` with imputed abundnce
setMethod("impute_missing", "spec_tbl_df", .impute_missing)

#' impute_missing
#' @inheritParams impute_missing
#' @return A `tbl` with imputed abundnce
setMethod("impute_missing", "tbl_df", .impute_missing)


#' Fill feature value if missing from element-feature pairs
#'
#' \lifecycle{maturing}
#'
#' @description fill_missing() takes as imput a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a `tbl` with an edditional adjusted value column. This method uses scaled counts if present.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name fill_missing
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column
#' @param .feature The name of the feature/gene column
#' @param .value The name of the feature/gene value column
#' @param fill_with A numerical value with which fill the mssing data points
#'
#' @details This function fills the value of missing element-feature pair using the median of the element group defined by the formula
#'
#' @return A `tbl` non-sparse value
#'
#'
#'
#'
#' @examples
#'
#'
#' res =
#' 	fill_missing(
#' 		nanny::counts_mini,
#' 	~ condition,
#' 	.element = element,
#' 	.feature = feature,
#' 	.value = count
#' )
#'
#'
#' @docType methods
#' @rdname fill_missing-methods
#'
#' @export
#'
#'
setGeneric("fill_missing", function(.data,
																			.element = NULL,
																			.feature = NULL,
																			.value = NULL,
																			fill_with)
	standardGeneric("fill_missing"))

# Set internal
.fill_missing = 	function(.data,
														.element = NULL,
														.feature = NULL,
														.value = NULL,
													fill_with)
{
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Validate data frame
	validation(.data, !!.element, !!.feature, !!.value)
	
	fill_NA_using_value(
		.data,
		.element = !!.element,
		.feature = !!.feature,
		.value = !!.value,
		fill_with)
}

#' fill_missing
#' @inheritParams fill_missing
#' @return A `tbl` with filld abundnce
setMethod("fill_missing", "spec_tbl_df", .fill_missing)

#' fill_missing
#' @inheritParams fill_missing
#' @return A `tbl` with filld abundnce
setMethod("fill_missing", "tbl_df", .fill_missing)
