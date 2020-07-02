#' Get clusters of elements (e.g., elements or features)
#'
#' \lifecycle{maturing}
#'
#' @description cluster_elements() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and identify clusters in the data.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom tidygate gate
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
#' @param transform A function to use to transform the data internally (e.g., log1p)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans
#' 
#' @details identifies clusters in the data, normally of elements.
#' This function returns a tibble with additional columns for the cluster annotation.
#' At the moment only k-means clustering is supported, the plan is to introduce more clustering methods.
#'
#' @return A tbl object with additional columns with cluster labels
#'
#' @examples
#'
#'
#'    cluster_elements(mtcars_tidy, car_model, feature, value, method="kmeans",	centers = 2)
#'
#' @docType methods
#' @rdname cluster_elements-methods
#' @export
#'
setGeneric("cluster_elements", function(.data,
																				.element,
																				.feature,
																				.value,
																				method,
																				of_elements = TRUE,
																				transform = NULL,
																				action = "add",
																				...)
	standardGeneric("cluster_elements"))

# Set internal
.cluster_elements = 		function(.data,
															 .element,
															 .feature,
															 .value,
															 method ,
															 of_elements = TRUE,
															 transform = NULL,
															 action = "add",
															 ...)
{
	
	# Comply with CRAN NOTES
	. = NULL
	
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# # Check if data rectangular
	# ifelse_pipe(
	# 	(.) %>% check_if_data_rectangular(!!.element,!!.feature,!!.value, type = "soft"),
	# 	~ .x %>% eliminate_sparse_features(!!.feature)
	# ) %>%
	
	if (method == "kmeans") {
		
		# Validate data frame
		validation(.data, !!.element, !!.feature, !!.value)
		
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
						),
					by = quo_names(.element)
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
						),
					by = quo_names(.element)
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
		
		# Validate data frame
		validation(.data, !!.element, !!.feature, !!.value)
		
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
						),
					by = quo_names(.element)
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
						),
					by = quo_names(.element)
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
	else if (method == "gate") {
		
		
		if (!action %in% c("add", "get", "only")) 	stop(
			"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
		
		.feature_names = quo_names(.feature)
		if(length(.feature_names) != 2) stop("nanny says: for gate clustering .feature must include exactly two columns. For example the first two PCAs.")
		
		.data %>%
			gate(
			.element = !!.element,
			.dim1 = !!as.symbol(.feature_names[1]),
			.dim2 = !!as.symbol(.feature_names[2]),
			action = action,
			...
		)
			

}
	else
		stop("nanny says: the only supported methods are \"kmeans\", \"SNN\" and \"gate\" ")
	
}

#' cluster_elements
#' @docType methods
#' @rdname cluster_elements-methods
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "spec_tbl_df", .cluster_elements)

#' cluster_elements
#' @docType methods
#' @rdname cluster_elements-methods
#' @return A tbl object with additional columns with cluster labels
#'
setMethod("cluster_elements", "tbl_df", .cluster_elements)

#' Dimension reduction of the feature value data
#'
#' \lifecycle{maturing}
#'
#' @description reduce_dimensions() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and calculates the reduced dimensional space of the feature value. The functions available are PCA, MDS (Robinson et al., 2010, <doi:10.1093/bioinformatics/btp616>), tSNE (Laurens van der Maaten, 2009)
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
#'   reduce_dimensions(mtcars_tidy, car_model, feature, value, method="PCA")
#'   
#'   reduce_dimensions(mtcars_tidy, car_model, feature, value, method="MDS")
#'   
#'   reduce_dimensions(mtcars_tidy, car_model, feature, value, method="tSNE")
#'
#'
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @export
#'
#'
setGeneric("reduce_dimensions", function(.data,
																				 .element,
																				 .feature,
																				 .value,
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
																.element,
																.feature,
																.value,
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
			
			.data %>%	dplyr::left_join(.data_processed,	by = quo_names(.element)) %>%
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
				
				dplyr::left_join(.data_processed,	by = quo_names(.element)) %>%
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
				dplyr::left_join(.data_processed,	by = quo_names(.element)) %>%
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
				
				dplyr::left_join(.data_processed,	by = quo_names(.element)) %>%
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
				dplyr::left_join(.data_processed,	by = quo_names(.element)	) %>%
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
				
				dplyr::left_join(.data_processed,	by = quo_names(.element)	) %>%
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
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "spec_tbl_df", .reduce_dimensions)

#' reduce_dimensions
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @return A tbl object with additional columns for the reduced dimensions
setMethod("reduce_dimensions", "tbl_df", .reduce_dimensions)

#' Rotate two dimensions (e.g., principal components) of an arbitrary angle
#'
#' \lifecycle{maturing}
#'
#' @description rotate_dimensions() takes as input a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the feature value.
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
#'  mtcars_tidy_MDS = reduce_dimensions(mtcars_tidy, car_model, feature, value, method="MDS")
#'  
#'  rotate_dimensions(mtcars_tidy_MDS, `Dim1`, `Dim2`, .element = car_model, rotation_degrees = 45)
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
																				 .element,
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
																.element,
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
			quo_names(dimension_1_column),
			rotation_degrees
		))
	if (quo_is_null(dimension_2_column_rotated))
		dimension_2_column_rotated = as.symbol(sprintf(
			"%s rotated %s",
			quo_names(dimension_2_column),
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
			dplyr::left_join(	.data_processed,	by = quo_names(.element)	) 
		
	}
	else if (action == "get"){
		
		.data %>%
			
			# Selecting the right columns
			select(
				!!.element,
				get_specific_annotation_columns(.data, !!.element)
			) %>%
			distinct() %>%
			
			dplyr::left_join(	.data_processed,	by = quo_names(.element)	) 
		
	}
	else if (action == "only") .data_processed
	else
		stop(
			"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' rotate_dimensions
#' @docType methods
#' @rdname rotate_dimensions-methods
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "spec_tbl_df", .rotate_dimensions)

#' rotate_dimensions
#' @docType methods
#' @rdname rotate_dimensions-methods
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("rotate_dimensions", "tbl_df", .rotate_dimensions)

#' Label points within a scatter plot drawing a gate
#'
#' \lifecycle{maturing}
#'
#' @description gate_dimensions() takes as input a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the feature value.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name gate_dimensions
#'
#'
#' @param .data A tibble
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param .dim1 A column symbol. The x dimension
#' @param .dim2 A column symbol. The y dimension
#' @param name A character string. The name of the new column
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans
#'
#' @details This function allow the user to label data points in inside a 2D gate.
#'
#' @return A tbl object with additional columns for the inside gate information. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
#'
#'
#' @examples
#'
#' \dontrun{
#' 
#'  mtcars_tidy_MDS = reduce_dimensions(mtcars_tidy, car_model, feature, value, method="MDS")
#'  
#'  gate_dimensions(mtcars_tidy_MDS, car_model, `Dim1`, `Dim2`)
#'  
#' }
#'
#' @docType methods
#' @rdname gate_dimensions-methods
#' @export
#'
setGeneric("gate_dimensions", function(.data,
																			 .element,
																			 .dim1,
																			 .dim2, 
																			 .color = NULL,
																			 .shape = NULL,
																			 .size = NULL,
																			 name = "inside_gate",
																			 action =	"add", ...)
	standardGeneric("gate_dimensions"))

# Set internal
.gate_dimensions = 		function(.data,
															.element,
															.dim1,
															.dim2, 
															name = "inside_gate",
																action =	"add", ...)
{
	
	# Get column names
	.element = enquo(.element)
	.dim1 = enquo(.dim1)
	.dim2 = enquo(.dim2)
	.color = enquo(.color)
	.shape = enquo(.shape)
	.size = enquo(.size)
	
	.data_processed =
		
		.data %>% 

		# Run calculation
		gate_dimensions_(
			.element = !!.element,
			.dim1 = !!.dim1,
			.dim2 = !!.dim2,
			.color = !!.color,
			.shape = !!.shape,
			.size = !!.size,
			name = name,
			...
		)
	
	if (action == "add"){
		
		.data %>%
			dplyr::left_join(	.data_processed,	by = quo_names(.element)	) 
		
	}
	else if (action == "get"){
		
		.data %>%
			
			# Selecting the right columns
			select(
				!!.element,
				get_specific_annotation_columns(.data, !!.element)
			) %>%
			distinct() %>%
			
			dplyr::left_join(	.data_processed,	by = quo_names(.element)	) 
		
	}
	else if (action == "only") .data_processed
	else
		stop(
			"nanny says: action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
		)
}

#' gate_dimensions
#' @docType methods
#' @rdname gate_dimensions-methods
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("gate_dimensions", "spec_tbl_df", .gate_dimensions)

#' gate_dimensions
#' @docType methods
#' @rdname gate_dimensions-methods
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
setMethod("gate_dimensions", "tbl_df", .gate_dimensions)


#' Drop redundant elements (e.g., elements) for which feature (e.g., feature/gene) aboundances are correlated
#'
#' \lifecycle{maturing}
#'
#' @description remove_redundancy() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | for correlation method, and returns a `tbl` with dropped elements (e.g., elements). The backend function used is widyr::pairwise_cor (David Robinson, 2020)
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
#' @param of_elements A boolean. In case the input is a nanny object, it indicates Whether the element column will be element or feature column
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param top An integer. How many top genes to select for correlation based method
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
#' remove_redundancy(mtcars_tidy, car_model, feature, value)
#'
#'
#' @docType methods
#' @rdname remove_redundancy-methods
#' @export
#'
#'
setGeneric("remove_redundancy", function(.data,
																				 .element,
																				 .feature,
																				 .value,
																				 of_elements = TRUE,
																				 correlation_threshold = 0.9,
																				 top = Inf,
																				 transform = NULL)
					 standardGeneric("remove_redundancy"))

# Set internal
.remove_redundancy = 	 function(.data,
																.element,
																.feature,
																.value,
																of_elements = TRUE,
																correlation_threshold = 0.9,
																top = Inf,
																transform = NULL)
{
	# Make col names
	.value = enquo(.value)
	.element = enquo(.element)
	.feature = enquo(.feature)
	
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

#' remove_redundancy
#' @docType methods
#' @rdname remove_redundancy-methods
#' @return A tbl object with with dropped recundant elements (e.g., elements).
setMethod("remove_redundancy", "spec_tbl_df", .remove_redundancy)

#' remove_redundancy
#' @docType methods
#' @rdname remove_redundancy-methods
#' @return A tbl object with with dropped recundant elements (e.g., elements).
setMethod("remove_redundancy", "tbl_df", .remove_redundancy)

#' Extract selected-column-wise information
#'
#' \lifecycle{maturing}
#'
#' @description subset() takes as input a `tbl` and returns a `tbl` with only selected-column-related columns
#'
#' @importFrom magrittr "%>%"
#'
#' @name subset
#'
#' @param .data A `tbl` 
#' @param .column The name of the column of interest
#'
#'
#' @details This functon extracts only selected-column-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to nanny function.
#'
#' @return A `tbl` object
#'
#'
#'
#'
#' @examples
#'
#' subset(mtcars_tidy,car_model)
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
	
	# Check if column present
	if(quo_names(.column) %in% colnames(.data) %>% all %>% `!`)
		stop("nanny says: some of the .column specified do not exist in the input data frame.")
		
	.data %>%
		
		# Selecting the right columns
		select(	!!.column,	get_specific_annotation_columns(.data, !!.column)	) %>%
		distinct()
	
}

#' subset
#' @docType methods
#' @rdname subset-methods
#' @return A `tbl` object
setMethod("subset",		"spec_tbl_df",		.subset)

#' subset
#' @docType methods
#' @rdname subset-methods
#' @return A `tbl` object
setMethod("subset",		"tbl_df",				.subset)

#' subset
#' @docType methods
#' @rdname subset-methods
#' @return A `tbl` object
setMethod("subset",		"tbl",			.subset)

#' Nest according to selected-column-wise information
#'
#' \lifecycle{maturing}
#'
#' @description nest_subset() takes as input a `tbl` and returns a nested `tbl` according to only selected-column-related columns
#'
#' @importFrom magrittr "%>%"
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom purrr imap
#' @importFrom rlang set_names
#' @importFrom tidyselect eval_select
#'
#' @name nest_subset
#'
#' @param .data A `tbl` 
#' @param ... The name of the columns of interest
#' @param .names_sep Deprecated by tidyr
#'
#'
#' @details This function extracts only selected-column-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to nanny function.
#'
#' @return A `tbl` object
#'
#'
#'
#'
#' @examples
#'
#' nest_subset(mtcars_tidy,data = -car_model)
#'
#'
#' @docType methods
#' @rdname nest_subset-methods
#' @export
#'
#'
setGeneric("nest_subset", function(.data, ..., .names_sep = NULL)
	standardGeneric("nest_subset"))

# Set internal
.nest_subset = 		function(.data, ..., .names_sep = NULL)	{
	
	# Make col names - from tidyr
	cols = enquos(...)
	cols <- map(cols, ~ names(eval_select(.x, .data)))
	cols <- map(cols, set_names)
	if (!is.null(.names_sep)) cols <- imap(cols, strip_names, .names_sep)
	asis <- setdiff(names(.data), unlist(cols))
	
	# Check if column present
	if(asis %in% colnames(.data) %>% all %>% `!`)
		stop("nanny says: some of the .column specified do not exist in the input data frame.")
	
	# Get my subset columns
	asis_subset = asis %>% c(get_specific_annotation_columns(.data, asis))
	
	# Apply nest on those
	tidyr::nest(.data, data = -c(asis_subset))
	
}

#' nest_subset
#' @docType methods
#' @rdname nest_subset-methods
#' @return A `tbl` object
setMethod("nest_subset",		"spec_tbl_df",		.nest_subset)

#' nest_subset
#' @docType methods
#' @rdname nest_subset-methods
#' @return A `tbl` object
setMethod("nest_subset",		"tbl_df",				.nest_subset)

#' nest_subset
#' @inheritParams nest_subset
#' @return A `tbl` object
setMethod("nest_subset",		"tbl",			.nest_subset)


#' Impute feature value if missing from element-feature pairs
#'
#' \lifecycle{maturing}
#'
#' @description impute_missing() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a `tbl` with an edditional adjusted value column. This method uses scaled counts if present.
#'
#' @importFrom rlang enquo
#' @importFrom rlang is_formula
#' @importFrom magrittr "%>%"
#'
#' @name impute_missing
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column
#' @param .feature The name of the feature/gene column
#' @param .value The name of the feature/gene value column
#' @param .formula A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_intrest + batch)
#' 
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
#'  impute_missing(mtcars_tidy, car_model, feature, value, ~1)
#'
#' @docType methods
#' @rdname impute_missing-methods
#'
#' @export
#'
#'
setGeneric("impute_missing", function(.data,
																				.element,
																				.feature,
																				.value,
																				.formula)
	standardGeneric("impute_missing"))

# Set internal
.impute_missing = 	function(.data,
														.element,
														.feature,
														.value,
														.formula)
{
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Sanity check formula
	formula_error_message = "nanny says: your formula does not look like one. Check it with rlang::is_formula"
	if(
		tryCatch(!is_formula(.formula), error=function(x) stop(formula_error_message))
	) stop(formula_error_message)
	
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
#' @docType methods
#' @rdname impute_missing-methods
#' @return A `tbl` with imputed abundnce
setMethod("impute_missing", "spec_tbl_df", .impute_missing)

#' impute_missing
#' @docType methods
#' @rdname impute_missing-methods
#' @return A `tbl` with imputed abundnce
setMethod("impute_missing", "tbl_df", .impute_missing)


#' Fill feature value if missing from element-feature pairs
#'
#' \lifecycle{maturing}
#'
#' @description fill_missing() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a `tbl` with an edditional adjusted value column. This method uses scaled counts if present.
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
#' fill_missing(mtcars_tidy, car_model, feature, value, fill_with = 0)
#'
#'
#' @docType methods
#' @rdname fill_missing-methods
#'
#' @export
#'
#'
setGeneric("fill_missing", function(.data,
																			.element,
																			.feature,
																			.value,
																			fill_with)
	standardGeneric("fill_missing"))

# Set internal
.fill_missing = 	function(.data,
														.element,
														.feature,
														.value,
													fill_with)
{
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Check the value is set
	if(length(fill_with)==0) stop("nanny says: the argument fill_with must not be empty.")
	
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
#' @docType methods
#' @rdname fill_missing-methods
#' @return A `tbl` with filled abundance
setMethod("fill_missing", "spec_tbl_df", .fill_missing)

#' fill_missing
#' @docType methods
#' @rdname fill_missing-methods
#' @return A `tbl` with filled abundance
setMethod("fill_missing", "tbl_df", .fill_missing)




#' Permute columns and nest data for each permutation
#'
#' \lifecycle{maturing}
#'
#' @description permute_nest() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a `tbl` with data nested for each permutation. The package used in the backend is gtools (Gregory R. Warnes, Ben Bolker, and Thomas Lumley, 2020)
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom gtools permutations
#'
#' @name permute_nest
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .names_from The columns to build the permutations on (e.g., c(col1, col2))
#' @param .values_from The columns to be nested for each permutation (e.g., c(col3, col4, col5))
#'
#' @details ...
#'
#' @return A nested `tbl` 
#'
#'
#'
#'
#' @examples
#'
#' permute_nest(mtcars_tidy, car_model, c(feature,value))
#'
#' @docType methods
#' @rdname permute_nest-methods
#'
#' @export
#'
#'
setGeneric("permute_nest", function(.data, .names_from, .values_from)
	standardGeneric("permute_nest"))

# Set internal
.permute_nest = 	function(.data, .names_from, .values_from){
	
	# Comply with CRAN NOTES
	. = NULL
	run = NULL
	# V1 = NULL
	# V2 - NULL
	
	# Column names
	.names_from = enquo(.names_from)
	.values_from = enquo(.values_from)
	
	# Check if multiple column inputted
	if(length(quo_names(.names_from))>1) 
		stop("nanny says: At the moment only one names column can be used to permute")
	
	factor_levels = .data %>% pull(!!.names_from) %>% unique
	
	.data %>% 
		pull(!!.names_from) %>%
		unique() %>%
		as.character() %>%
		permutations(n = length(.), r = 2, v = .) %>%
		as_tibble() %>%
		unite("run", c("V1", "V2"), remove = FALSE, sep="___") %>%
		gather(which, !!.names_from, -run) %>%
		select(-which) %>%
		left_join(.data %>% select(!!.names_from, !!.values_from), by = quo_names(.names_from)) %>%
		nest(data = -run) %>%
		separate(run, sprintf("%s_%s", quo_names(.names_from), 1:2 ), sep="___") %>%
		
		# Introduce levels
		mutate_at(vars(1:2),function(x) factor(x, levels = factor_levels))
	
}

#' permute_nest
#' @docType methods
#' @rdname permute_nest-methods
#' @return A `tbl` with filled abundance
setMethod("permute_nest", "spec_tbl_df", .permute_nest)

#' permute_nest
#' @docType methods
#' @rdname permute_nest-methods
#' @return A `tbl` with filled abundance
setMethod("permute_nest", "tbl_df", .permute_nest)




#' Combine columns and nest data for each permutation
#'
#' \lifecycle{maturing}
#'
#' @description combine_nest() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a `tbl` with data nested for each combination  The package used in the backend is gtools (Gregory R. Warnes, Ben Bolker, and Thomas Lumley, 2020)
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom gtools combinations
#'
#' @name combine_nest
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .names_from The columns to build the permutations on (e.g., c(col1, col2))
#' @param .values_from The columns to be nested for each permutation (e.g., c(col3, col4, col5))
#'
#' @details ...
#'
#' @return A nested `tbl` 
#'
#'
#'
#'
#' @examples
#'
#' combine_nest(mtcars_tidy, car_model, c(feature,value))
#'
#'
#' @docType methods
#' @rdname combine_nest-methods
#'
#' @export
#'
#'
setGeneric("combine_nest", function(.data, .names_from, .values_from)
	standardGeneric("combine_nest"))

# Set internal
.combine_nest = function(.data, .names_from, .values_from){
	
	# Comply with CRAN NOTES
	. = NULL
	run = NULL
	# V1 = NULL
	# V2 - NULL
	
	# Column names
	.names_from = enquo(.names_from)
	.values_from = enquo(.values_from)
	
	factor_levels = .data %>% pull(!!.names_from) %>% unique
	
	# Check if multiple column inputted
	if(length(quo_names(.names_from))>1) 
		stop("nanny says: At the moment only one names column can be used to permute")
	
	factor_levels = .data %>% pull(!!.names_from) %>% unique
	
	.data %>% 
		pull(!!.names_from) %>%
		unique() %>%
		as.character() %>%
		gtools::combinations(n = length(.), r = 2, v = .) %>%
		as_tibble() %>%
		unite("run", c("V1", "V2"), remove = FALSE, sep="___") %>%
		gather(which, !!.names_from, -run) %>%
		select(-which) %>%
		left_join(.data %>% select(!!.names_from, !!.values_from), by = quo_names(.names_from)) %>%
		nest(data = -run) %>%
		separate(run, sprintf("%s_%s", quo_names(.names_from), 1:2), sep="___") %>%
		
		# Introduce levels
		mutate_at(vars(1:2),function(x) factor(x, levels = factor_levels))
	
}

#' combine_nest
#' @docType methods
#' @rdname combine_nest-methods
#' @return A `tbl` with filled abundance
setMethod("combine_nest", "spec_tbl_df", .combine_nest)

#' combine_nest
#' @docType methods
#' @rdname combine_nest-methods
#' @return A `tbl` with filled abundance
setMethod("combine_nest", "tbl_df", .combine_nest)



#' Keep top variable features across elements
#'
#' \lifecycle{maturing}
#'
#' @description keep_variable() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a `tbl` with the filtered most variable features. The formula used is from limma::plotMDS (Robinson et al., 2010, <doi:10.1093/bioinformatics/btp616>)
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name keep_variable
#'
#' @param .data A `tbl`
#' @param .element A character name of the element column
#' @param .feature A character name of the transcript/gene column
#' @param .value A character name of the read count column
#' @param top An integer. How many top genes to select
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#'
#' @details ...
#'
#' @return A `tbl` with filtered features
#'
#'
#'
#'
#' @examples
#'
#' keep_variable(mtcars_tidy, car_model, feature, value, top=10)
#'
#'
#' @docType methods
#' @rdname keep_variable-methods
#'
#' @export
#'
#'
setGeneric("keep_variable", function(.data,
																		 .element,
																		 .feature,
																		 .value,
																		 top = Inf,
																		 transform = NULL)
	standardGeneric("keep_variable"))

# Set internal
.keep_variable = function(.data,
													.element,
													.feature,
													.value,
													top = Inf,
													transform = NULL) {
	
	# Comply with CRAN NOTES
	. = NULL
	value = NULL
	variable = NULL
	
	
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Check that column names do not have the reserved pattern "___"
	if(.data %>% colnames %>% grep("___", .) %>% length %>% `>` (0))
		stop("nanny says: your column names cannot include the pattern \"___\" that is reserved for internal manipulation")
	
	
	# Manage Inf
	top = min(top, .data %>% select(!!.feature) %>% distinct %>% nrow)
	
	x =
		.data %>%
		select(!!.element, !!.feature, !!.value) %>%
		distinct %>%
		
		# Check if tranfrom is needed
		ifelse_pipe(
			is_function(transform),
			~ .x %>%
				mutate(!!.value := !!.value %>%  transform()) %>%
				
				# Check is log introduced -Inf
				ifelse_pipe(pull(.,!!.value) %>% min %>% equals(-Inf),
										~ stop(
											"nanny says: you applied a transformation that introduced negative infinite .value, was it log? If so please use log1p."
										))
		) %>%
		
		pivot_wider(names_from = !!.element, values_from = !!.value, names_sep = "___") %>%
		as_matrix(rownames = !!.feature)
	
	s <- rowMeans((x - rowMeans(x)) ^ 2)
	o <- order(s, decreasing = TRUE)
	x <- x[o[1L:top], , drop = FALSE]
	

	.data %>% inner_join(
		rownames(x) %>% as_tibble() %>% separate(col = value, into = quo_names(.feature), sep = "___")
	) 
	
}

#' keep_variable
#' @docType methods
#' @rdname keep_variable-methods
#' @return A `tbl` with filled abundance
setMethod("keep_variable", "spec_tbl_df", .keep_variable)

#' keep_variable
#' @docType methods
#' @rdname keep_variable-methods
#' @return A `tbl` with filled abundance
setMethod("keep_variable", "tbl_df", .keep_variable)



#' Keep rows corresponding of a lower triangular matrix built from two columns
#'
#' \lifecycle{maturing}
#'
#' @description lower_triangular() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and returns a filtered `tbl` 
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name lower_triangular
#'
#' @param .data A `tbl`
#' @param .col1 A column name
#' @param .col2 A column name
#' @param .value A column names of the value column
#'
#' @details ...
#'
#' @return A `tbl` with filtered rows
#'
#'
#'
#'
#' @examples
#'
#' library(dplyr)
#' library(purrr)
#' library(tidyr)
#' 
#' mtcars_tidy_permuted = 
#'   mtcars_tidy %>% 
#'   filter(feature == "mpg") %>% 
#'   head(5) %>% 
#'   permute_nest(car_model, c(feature,value))
#' 
#' mtcars_tidy_permuted %>%
#'  # Summarise mpg
#'  mutate(data = map(data, ~ .x %>% summarise(mean(value)))) %>%
#'	unnest(data) %>%
#'	
#'	# Lower triangular
#'	lower_triangular(car_model_1, car_model_2,  `mean(value)`)
#'
#' @docType methods
#' @rdname lower_triangular-methods
#'
#' @export
#'
#'
setGeneric("lower_triangular", function(.data, .col1, .col2, .value)
	standardGeneric("lower_triangular"))

# Set internal
.lower_triangular = function(.data, .col1, .col2, .value){
	
	# Comply with CRAN NOTES
	. = NULL
	
	# Column names
	.col1 = enquo(.col1)
	.col2 = enquo(.col2)
	.value = enquo(.value)

	
	#levs = .data %>% pull(!!.col1) %>% levels

	.data %>%
		
		# Check if duplicated elements
		error_if_duplicated_genes(!!.col1,!!.col2,!!.value) %>%
	
		select(!!.col1, !!.col2,    !!.value) %>%
		spread(!!.col2 ,   !!.value) %>%
		as_matrix(rownames = quo_names(.col1)) %>%

		# Drop upper triangular
		{ ma = (.); ma[lower.tri(ma)] <- NA; ma} %>%

		as_tibble(rownames = quo_names(.col1)) %>%
		gather(!!.col2, !!.value, -!!.col1) %>%
		mutate(
			# !!.col1 := factor(!!.col1, levels = levs),
			# !!.col2 := factor(!!.col2, levels = levs),
			
			!!.col1 := factor(!!.col1),
			!!.col2 := factor(!!.col2),
		) %>%
		drop_na %>%
		
	# Reattach col1 col2 wise annotation
	left_join(.data %>% select(-!!.value) %>% subset(c(!!.col1, !!.col2)), by=c(quo_name(.col1), quo_name(.col2)))
	
}

#' lower_triangular
#' @docType methods
#' @rdname lower_triangular-methods
#' @return A `tbl` with filled abundance
setMethod("lower_triangular", "spec_tbl_df", .lower_triangular)

#' lower_triangular
#' @docType methods
#' @rdname lower_triangular-methods
#' @return A `tbl` with filled abundance
setMethod("lower_triangular", "tbl_df", .lower_triangular)


#' Get matrix from tibble
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#' @importFrom rlang quo_is_symbolic
#' @importFrom purrr when
#' 
#'
#' @param .data A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#' @param sep_rownames A character with which multiple columns are united if rownames is a column array (e.g., rownames = c(col1, col2))
#'
#' @return A matrix
#'
#' @examples
#'
#'  library(dplyr)
#'  library(tidyr)
#'  select(mtcars_tidy, car_model, feature, value) %>%
#' 	spread(feature, value) %>%
#' 	as_matrix(rownames = car_model) 
#' 	
#' @docType methods
#' @rdname as_matrix-methods
#'
#' @export
setGeneric("as_matrix", function(.data,
																 rownames = NULL,
																 do_check = TRUE,
																 sep_rownames = "___")
	standardGeneric("as_matrix"))

# Set internal
.as_matrix = function(.data,
											rownames = NULL,
											do_check = TRUE,
											sep_rownames = "___") {
	
	# Comply with CRAN NOTES
	variable = NULL
	
	rownames = enquo(rownames)

	
	.data %>%
		
		# Through warning if data frame is not numerical beside the rownames column (if present)
		ifelse_pipe(
			do_check &&
				.data %>%
				# If rownames defined eliminate it from the data frame
				ifelse_pipe(!quo_is_null(rownames), ~ .x %>% select(-!!rownames), ~ .x) %>%
				dplyr::summarise_all(class) %>%
				tidyr::gather(variable, class) %>%
				pull(class) %>%
				unique() %>%
				`%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
			~ {
				warning("nanny says: there are NON-numerical columns, the matrix will NOT be numerical")
				.x
			}
		) %>%
		
		# If rownames multiple enquo (e.g., c(col1, col2)) merge them
		when(!quo_is_null(rownames) ~ (.) %>% unite(col = "rn", !!rownames, sep = sep_rownames), ~ (.)) %>%

		as.data.frame() %>%
		
		# Deal with rownames column if present
		ifelse_pipe(
			!quo_is_null(rownames),
			~ .x %>%
				set_rownames(.x %>% pull(rn)) %>%
				select(-rn)
		) %>%
		
		# Convert to matrix
		as.matrix()
}

#' as_matrix
#' @docType methods
#' @rdname as_matrix-methods
#' @return A `tbl` with filled abundance
setMethod("as_matrix", "spec_tbl_df", .as_matrix)

#' as_matrix
#' @docType methods
#' @rdname as_matrix-methods
#' @return A `tbl` with filled abundance
setMethod("as_matrix", "tbl_df", .as_matrix)