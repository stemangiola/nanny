#' Get K-mean clusters to a tibble
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats kmeans
#' @importFrom rlang :=
#' @importFrom rlang is_function
#' @importFrom magrittr `%$%`
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally elements)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_elements A boolean
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
#'
get_clusters_kmeans_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .value = NULL,
					 of_elements = TRUE,
					 transform = NULL,
					 ...) {
		
		# Comply with CRAN NOTES
		. = NULL
		seurat_clusters = NULL
		cluster = NULL
		cluster_kmeans = NULL
		
		# Check that column names do not have the reserved pattern "___"
		if(.data %>% colnames %>% grep("___", .) %>% length %>% `>` (0))
			stop("nanny says: your column names cannot include the pattern \"___\" that is reserved for internal manipulation")
		
		# Check if centers is in dots
		dots_args = rlang::dots_list(...)
		if ("centers" %in% names(dots_args) %>% `!`)
			stop("nanny says: for kmeans you need to provide the \"centers\" integer argument")
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		.data %>%
			
			# Prepare data frame
			select(!!.feature,!!.element,!!.value) %>%
			distinct() %>%
			
			# Check if tranfrom is needed
			ifelse_pipe(
				is_function(transform),
				~ .x %>% 
					mutate(!!.value := !!.value %>%  transform()) %>%
					
					# Check is log introduced -Inf
					ifelse_pipe(
						pull(., !!.value) %>% min %>% equals(-Inf), 
						~ stop("nanny says: you applied a transformation that introduced negative infinite .value, was it log? If so please use log1p.")
					)
			) %>%
			
			# Prepare data frame for return
			pivot_wider(names_from = !!.feature, values_from = !!.value, names_sep = "___") %>%
			as_matrix(rownames = !!.element) %>%
			
			# Wrap the do.call because of the centers check
			{
				do.call(kmeans, list(x = (.), iter.max = 1000) %>% c(dots_args))
			}	 %$%
			cluster %>%
			as.list() %>%
			as_tibble() %>%
			pivot_longer(names_to = quo_names(.element), cols=everything(), names_sep = when(length(quo_names(.element)), (.) > 1 ~ "___", ~ NULL), values_to = "cluster_kmeans") %>%
			mutate(cluster_kmeans = cluster_kmeans %>% as.factor()) %>%
			
			# Attach attributes
			reattach_internals(.data)
	}

#' Get SNN shared nearest neighbour clusters to a tibble
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom rlang is_function
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally elements)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_elements A boolean
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
get_clusters_SNN_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .value,
					 of_elements = TRUE,
					 transform = NULL,
					 ...) {
		
		# Comply with CRAN NOTES
		rn = NULL
		seurat_clusters = NULL
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		# Check if package is installed, otherwise install
		if (find.package("dbscan", quiet = T) %>% length %>% equals(0)) {
	 stop("nanny says: dbscan is necessary for this operation. Please install it with 	install.packages(\"dbscan\", repos = \"https://cloud.r-project.org\")")
		}
		
		# Check if centers is in dots
		dots_args = rlang::dots_list(...)
		#if ("k" %in% names(dots_args) %>% `!`) dots_args = dots_args %>% c(list(k = 20))
		if ("eps" %in% names(dots_args) %>% `!`) dots_args = dots_args %>% c(list(eps = 0.15))
		if ("minPts" %in% names(dots_args) %>% `!`) dots_args = dots_args %>% c(list(minPts = 5))
		
		my_df =
			.data %>%
			
			# Prepare data frame
			select(!!.element,!!.feature,!!.value) %>%
			distinct() %>%
			
			# Check if transform is needed
			ifelse_pipe(
				is_function(transform),
				~ .x %>% 
					mutate(!!.value := !!.value %>%  transform()) %>%
					
					# Check is log introduced -Inf
					ifelse_pipe(
						pull(., !!.value) %>% min %>% equals(-Inf), 
						~ stop("nanny says: you applied a transformation that introduced negative infinite .value, was it log? If so please use log1p.")
					)
			) %>%
			
			# Prepare data frame for return
			pivot_wider(names_from = !!.feature, values_from = !!.value, names_sep = "___") 


		.data %>%
			select(!!.element) %>%
			distinct() %>% 
			#arrange(!!.element) %>%
			bind_cols(
					my_df %>%
						as_matrix(rownames = !!.element) %>%
						
						# Scale
						scale() %>%
					
						# Wrap the do.call because of the centers check
						{
							do.call(dbscan::dbscan, list(x = (.)) %>% c(dots_args))
						}	 %$%
							
						cluster %>%
						as_tibble() %>%
						rename(cluster_SNN = value) %>%
						mutate(cluster_SNN = as.factor(cluster_SNN))
			)

	}




#' Get dimensionality information to a tibble using MDS
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom rlang is_function
#' @importFrom magrittr `%$%`
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param top An integer. How many top genes to select
#' @param of_elements A boolean
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#'
#' @return A tibble with additional columns
#'
#'
get_reduced_dimensions_MDS_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .value = NULL,
					 .dims = 2,
					 top = Inf,
					 of_elements = TRUE,
					 transform = NULL) {
		
		# Comply with CRAN NOTES
		. = NULL
		cmdscale.out = NULL
		rn = NULL
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		# Get components from dims
		components = 1:.dims
		
		mds_object =
			.data %>%
			
			# Filter lowly transcribed (I have to avoid the use of scaling function)
			# keep_abundant(!!.element, !!.feature,!!.value) %>%
			select(!!.feature,!!.element,!!.value) %>%
			distinct %>%
			
			# Check if tranfrom is needed
			ifelse_pipe(
				is_function(transform),
				~ .x %>% 
					mutate(!!.value := !!.value %>%  transform()) %>%
					
					# Check is log introduced -Inf
					ifelse_pipe(
						pull(., !!.value) %>% min %>% equals(-Inf), 
						~ stop("nanny says: you applied a transformation that introduced negative infinite .value, was it log? If so please use log1p.")
					)
			) %>%
			
			# Stop any column is not if not numeric or integer
			ifelse_pipe(
				(.) %>% select(!!.value) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
				~ stop("nanny says: .value must be numerical or integer")
			) %>%
			pivot_wider(names_from = !!.element, values_from = !!.value, names_sep = "___") %>%
			as_matrix(rownames = !!.feature, do_check = FALSE) %>%
			limma::plotMDS(ndim = .dims, plot = FALSE, top = top)
		
		# Pase results
		mds_object %$%	cmdscale.out %>%
			as.data.frame %>%
			as_tibble(rownames = "rn") %>%
			separate(col = rn, into = quo_names(.element), sep = "___") %>%

			setNames(c(quo_names(.element), sprintf("Dim%s", 1:.dims))) %>%
			
			
			# Attach attributes
			reattach_internals(.data) %>%
			
			# Add raw object
			attach_to_internals(mds_object, "MDS") %>%
			# Communicate the attribute added
			{
				message("nanny says: to access the raw results do `attr(..., \"internals\")$MDS`")
				(.)
			}
		
	}

#' Get principal component information to a tibble using PCA
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats prcomp
#' @importFrom rlang is_function
#' @importFrom magrittr `%$%`
#' @importFrom utils capture.output
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param top An integer. How many top genes to select
#' @param of_elements A boolean
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param scale A boolean
#' @param ... Further parameters passed to the function prcomp
#'
#' @return A tibble with additional columns
#'
#'
get_reduced_dimensions_PCA_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .value  = NULL,
					 .dims = 2,
					 top = Inf,
					 of_elements = TRUE,
					 transform = NULL,
					 scale = FALSE,
					 ...) {
		
		# Comply with CRAN NOTES
		. = NULL
		sdev = NULL
		name = NULL
		value = NULL
		rn = NULL
		Y = NULL
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		# Get components from dims
		components = 1:.dims
		
		prcomp_obj =
			.data %>%
			
			# Prepare data frame
			select(!!.feature,!!.element,!!.value) %>%
			distinct %>%
			
			# Check if tranfrom is needed
			ifelse_pipe(
				is_function(transform),
				~ .x %>% 
					mutate(!!.value := !!.value %>%  transform()) %>%
					
					# Check is log introduced -Inf
					ifelse_pipe(
						pull(., !!.value) %>% min %>% equals(-Inf), 
						~ stop("nanny says: you applied a transformation that introduced negative infinite .value, was it log? If so please use log1p.")
					)
			) %>%
			
			# Stop any column is not if not numeric or integer
			ifelse_pipe(
				(.) %>% select(!!.value) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
				~ stop("nanny says: .value must be numerical or integer")
			) %>%
			
			# Filter most variable genes
			keep_variable(!!.element,!!.feature,!!.value, top) %>%
			
			pivot_wider(names_from = !!.element, values_from = !!.value, names_sep = "___") %>%
			
			drop_na %>% # Is this necessary?
			
			# check that there are non-NA genes for enough elements
			ifelse2_pipe(# First condition
				(.) %>% nrow == 0,
				
				# Second condition
				(.) %>% nrow < 100,
				
				# First function
				~ stop(
					"nanny says: In calculating PCA there is no gene that have non NA values is all elements"
				),
				
				# Second function
				~ {
					warning(
						"
						nanny says: In PCA correlation there is < 100 genes that have non NA values is all elements.
						The correlation calculation would not be reliable,
						we suggest to partition the dataset for element clusters.
						"
					)
					.x
				}) %>%
			
			# Transform to matrix
			as_matrix(rownames = !!.feature, do_check = FALSE) %>%
			
			# Calculate principal components
			prcomp(scale = scale, ...)
		
		prcomp_obj %>%
			
			# Anonymous function - Prints fraction of variance
			# input: PCA object
			# output: PCA object
			{
				message("Fraction of variance explained by the selected principal components")
				
				(.) %$% sdev %>% `^` (2) %>% # Eigen value
					`/` (sum(.)) %>%
					`[` (components) %>%
					enframe() %>%
					select(-name) %>%
					rename(`Fraction of variance` = value) %>%
					mutate(PC = components) %>%
					as.data.frame() %>%
					
					# Print as message
					capture.output() %>% paste0(collapse = "\n") %>% message()
				
				(.)
				
			} %$%
			
			# Parse the PCA results to a tibble
			rotation %>%
			as_tibble(rownames = "rn") %>%
			separate(col = rn, into = quo_names(.element), sep = "___") %>%
			select(!!.element, sprintf("PC%s", components)) %>%
			
			# Attach attributes
			reattach_internals(.data) %>%
			
			# Add raw object
			attach_to_internals(prcomp_obj, "PCA") %>%
			# Communicate the attribute added
			{
				message("nanny says: to access the raw results do `attr(..., \"internals\")$PCA`")
				(.)
			}
		
	}

#' Get principal component information to a tibble using tSNE
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom rlang is_function
#' @importFrom magrittr `%$%`
#' @importFrom Rtsne Rtsne
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param top An integer. How many top genes to select
#' @param of_elements A boolean
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#' @param ... Further parameters passed to the function Rtsne
#'
#' @return A tibble with additional columns
#'
get_reduced_dimensions_TSNE_bulk <-
	function(.data,
					 .element = NULL,
					 .feature = NULL,
					 .value = NULL,
					 .dims = 2,
					 top = Inf,
					 of_elements = TRUE,
					 transform = NULL,
					 ...) {
		
		# Comply with CRAN NOTES
		. = NULL
		Y = NULL
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		# Evaluate ...
		arguments <- list(...)
		if (!"check_duplicates" %in% names(arguments))
			arguments = arguments %>% c(check_duplicates = TRUE)
		if (!"verbose" %in% names(arguments))
			arguments = arguments %>% c(verbose = TRUE)
		if (!"dims" %in% names(arguments))
			arguments = arguments %>% c(dims = .dims)
		
		# # Check if package is installed, otherwise install
		# if (find.package("Rtsne", quiet = T) %>% length %>% equals(0)) {
		# 	 stop("nanny says: Rtsne is necessary for this operation. Please install it with 	install.packages(\"Rtsne\")")
		# }
		
		# Set perprexity to not be too high
		if (!"perplexity" %in% names(arguments))
			arguments = arguments %>% c(perplexity = ((
				.data %>% select(!!.element) %>% distinct %>% nrow %>% sum(-1)
			) / 3 / 2) %>% floor() %>% min(30))
		
		# If not enough elements stop
		if (arguments$perplexity <= 2)
			stop("nanny says: You don't have enough elements to run tSNE")
		
		# Calculate the most variable genes, from plotMDS Limma
		
		
		df_tsne =
			.data %>%
			
			# Check if duplicates
			error_if_duplicated_genes(!!.element,!!.feature,!!.value)  %>%
			
			# Filter NA symbol
			filter(!!.feature %>% is.na %>% `!`) %>%
			
			# Prepare data frame
			select(!!.feature,!!.element,!!.value) %>%
			distinct %>%
			
			# Check if tranfrom is needed
			ifelse_pipe(
				is_function(transform),
				~ .x %>% 
					mutate(!!.value := !!.value %>%  transform()) %>%
					
					# Check is log introduced -Inf
					ifelse_pipe(
						pull(., !!.value) %>% min %>% equals(-Inf), 
						~ stop("nanny says: you applied a transformation that introduced negative infinite .value, was it log? If so please use log1p.")
					)
			) %>%
			
			# Filter most variable genes
			keep_variable(!!.element,!!.feature,!!.value, top) %>%
			
			pivot_wider(names_from = !!.feature, values_from = !!.value, names_sep = "___") %>%

			# select(-element) %>%
			# distinct %>%
			as_matrix(rownames = quo_names(.element))
		
		do.call(Rtsne, c(list(df_tsne), arguments)) %$%
			Y %>%
			as_tibble(.name_repair = "minimal") %>%
			setNames(c("tSNE1", "tSNE2")) %>%
			
			# add element name
			dplyr::mutate(!!.element := df_tsne %>% rownames) %>%
			select(!!.element, everything()) %>%
			
			# Attach attributes
			reattach_internals(.data)
		
	}

#' Get rotated dimensions of two principal components or MDS dimension of choice, of an angle
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang quo_is_null
#'
#'
#' @param .data A tibble
#' @param dimension_1_column A column symbol. The column of the dimension 1
#' @param dimension_2_column   A column symbol. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param of_elements A boolean
#' @param dimension_1_column_rotated A column symbol. The column of the dimension 1 rotated
#' @param dimension_2_column_rotated   A column symbol. The column of the dimension 2 rotated
#'
#' @return A tibble with additional rotated columns
#'
#'
get_rotated_dimensions =
	function(.data,
					 dimension_1_column,
					 dimension_2_column,
					 rotation_degrees,
					 .element = NULL,
					 of_elements = TRUE,
					 dimension_1_column_rotated = NULL,
					 dimension_2_column_rotated = NULL) {
		
		# Comply with CRAN NOTES
		. = NULL
		Y = NULL
		rotated_dimensions = NULL
		value = NULL
		
		# Get column names
		.element = enquo(.element)
		dimension_1_column = enquo(dimension_1_column)
		dimension_2_column = enquo(dimension_2_column)
		dimension_1_column_rotated = enquo(dimension_1_column_rotated)
		dimension_2_column_rotated = enquo(dimension_2_column_rotated)
		
		if (.data %>%
				select(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
				distinct %>%
				
				# Count
				group_by_at(vars(!!.element, !!dimension_1_column, !!dimension_2_column)) %>%
				tally() %>%
				ungroup() %>%
				
				# Check
				pull(n) %>%
				max %>%
				`>` (1))
			stop(sprintf(
				"nanny says: %s must be unique for each row for the calculation of rotation",
				quo_names(.element)
			))
		
		# Sanity check of the angle selected
		if (rotation_degrees %>% between(-360, 360) %>% `!`)
			stop("nanny says: rotation_degrees must be between -360 and 360")
		
		# Return
		.data %>%
			select(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
			distinct %>%
			as_matrix(rownames = !!.element) %>% t %>%
			rotation(rotation_degrees) %>%
			as_tibble() %>%
			mutate(`rotated_dimensions` =
						 	c(
						 		quo_names(dimension_1_column_rotated),
						 		quo_names(dimension_2_column_rotated)
						 	)) %>%
			pivot_longer(names_to = quo_names(.element),values_to = "value", cols = -`rotated_dimensions`, names_sep = when(length(quo_names(.element)), (.) > 1 ~ "___", ~ NULL)) %>%
			pivot_wider(names_from = `rotated_dimensions`, values_from = value) %>%
			

			
			# Attach attributes
			reattach_internals(.data)
		
	}

#' Get points within a user drawn gate
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom graphics plot
#'
#' @param .data A tibble
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param .dim1 A column symbol. The x dimension
#' @param .dim2 A column symbol. The y dimension
#' @param ... Further parameters passed to the function kmeans
#'
#' @return A tibble with additional columns
#'
gate_dimensions_ <-
	function(.data,
					 .element,
					 .dim1,
					 .dim2, name = "inside_gate", ...) {
		
		# Comply with CRAN NOTES
		. = NULL
		value = NULL
		
		# Get column names
		.element = enquo(.element)
		.dim1 = enquo(.dim1)
		.dim2 = enquo(.dim2)
		
		# Check if package is installed, otherwise install
		if (find.package("gatepoints", quiet = T) %>% length %>% equals(0)) {
			stop("nanny says: gatepoints is necessary for this operation. Please install it with 	install.packages(\"gatepoints\", repos = \"https://cloud.r-project.org\")")
		}
		
		if (.data %>%
				select(!!.element, !!.dim1, !!.dim2) %>%
				distinct %>%
				
				# Count
				group_by_at(vars(!!.element, !!.dim1, !!.dim2)) %>%
				tally() %>%
				ungroup() %>%
				
				# Check
				pull(n) %>%
				max %>%
				`>` (1))
		stop(sprintf(
			"nanny says: %s must be unique for each row for the calculation",
			quo_names(.element)
		))
		
		# Return
		my_df = 
			.data %>%
			select(!!.element, !!.dim1, !!.dim2) %>%
			distinct 
		
		my_matrix	=
			my_df %>%
			as_matrix(rownames = !!.element) 
		
		my_matrix %>% plot()
		
		my_df %>%
			select(-c(!!.dim1, !!.dim2)) %>%
			left_join(
				
				# Return clustering
				my_matrix %>% 
					gatepoints::fhs(mark = TRUE, ...) %>%
					as.character %>%
					as_tibble() %>%
					
					# Reconstitute columns
					separate(value, quo_names(.element), sep="___") %>%
					
					mutate(!!as.symbol(name) := T),
				by = quo_names(.element)
			) %>%
			
			mutate(!!as.symbol(name) := if_else(!!as.symbol(name) %>% is.na, F, !!as.symbol(name)))
		
	}

#' Drop redundant elements (e.g., elements) for which feature (e.g., genes) aboundances are correlated
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom rlang is_function
#' @importFrom widyr pairwise_cor
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param correlation_threshold A real number between 0 and 1
#' @param top An integer. How many top genes to select
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param of_elements A boolean
#' @param transform A function to use to tranforma the data internalli (e.g., log1p)
#'
#' @return A tibble with redundant elemens removed
#'
#'
remove_redundancy_elements_through_correlation <- function(.data,
																													 .element = NULL,
																													 .feature = NULL,
																													 .value = NULL,
																													 correlation_threshold = 0.9,
																													 top = Inf,
																													 of_elements = TRUE,
																													 transform = NULL) {
	# Comply with CRAN NOTES
	. = NULL
	value = NULL
	element = NULL
	feature = NULL
	correlation = NULL
	item1 = NULL
	
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Get the redundant data frame
	.data.correlated =
		.data %>%
		
		# Stop if there are duplicated features
		error_if_duplicated_genes(!!.element,!!.feature,!!.value) %>%
		
		# Prepare the data frame
		select(!!.feature,!!.element,!!.value) %>%
		
		# Filter variable genes
		keep_variable(!!.element,!!.feature,!!.value, top = top) %>%
		
		# Check if tranfrom is needed
		ifelse_pipe(
			is_function(transform),
			~ .x %>% 
				mutate(!!.value := !!.value %>%  transform()) %>%
				
				# Check is log introduced -Inf
				ifelse_pipe(
					pull(., !!.value) %>% min %>% equals(-Inf), 
					~ stop("nanny says: you applied a transformation that introduced negative infinite .value, was it log? If so please use log1p.")
				)
		) %>%
		
		distinct() %>%
		pivot_wider(names_from = !!.element, values_from = !!.value, names_sep = "___") %>%
		
		drop_na() %>%
		
		# check that there are non-NA genes for enough elements
		ifelse2_pipe(# First condition
			(.) %>% nrow == 0,
			
			# Second condition
			(.) %>% nrow < 100,
			
			# First function
			~ stop(
				"nanny says: In calculating correlation there is no gene that have non NA values is all elements"
			),
			
			# Second function
			~ {
				warning(
					"
					nanny says: In calculating correlation there is < 100 genes that have non NA values is all elements.
					The correlation calculation would not be reliable,
					we suggest to partition the dataset for element clusters.
					"
				)
				.x
			}) %>%
		
		# Prepare the data frame
		pivot_longer(names_to = quo_names(.element),values_to = quo_names(.value), cols = -!!.feature, names_sep = when(length(quo_names(.element)), (.) > 1 ~ "___", ~ NULL)) %>%
		mutate_if(is.factor, as.character) %>%
		
		# Unite columns and rename, needed for widyr
		unite("element", !!.element, sep="___") %>%
		unite("feature", !!.feature, sep="___") %>%
		rename(value := !!.value) %>%
		
		# Run pairwise correlation and return a tibble
		pairwise_cor(
			element,
			feature,
			value,
			sort = TRUE,
			diag = FALSE,
			upper = FALSE
		) %>%
		filter(correlation > correlation_threshold) %>%
		select(item1) %>%
		distinct %>%
		
		# Reconstitute columns
		separate(item1, quo_names(.element), sep="___") 
	
	# Return non redudant data frame
	.data %>% anti_join(.data.correlated) %>%
		
		# Attach attributes
		reattach_internals(.data)
}

#' Identifies the closest pairs in a MDS contaxt and return one of them
#' 
#' @keywords internal
#'
#' @importFrom stats setNames
#' @importFrom stats dist
#'
#' @param .data A tibble
#' @param Dim_a_column A column symbol. The column of one principal component
#' @param Dim_b_column A column symbol. The column of another principal component
#' @param .element A column symbol. The column that is represents entities to cluster (i.e., normally elements)
#' @param of_elements A boolean
#'
#' @return A tibble with pairs dropped
#'
#'
remove_redundancy_elements_though_reduced_dimensions <-
	function(.data,
					 Dim_a_column,
					 Dim_b_column,
					 .element = NULL,
					 of_elements = TRUE) {
		
		# Comply with CRAN NOTES
		. = NULL
		element_a = NULL
		element_b = NULL
		element_1 = NULL
		element_2 = NULL
		
		# This function identifies the closest pairs and return one of them
		
		# Get column names
		.element = enquo(.element)
		
		Dim_a_column = enquo(Dim_a_column)
		Dim_b_column = enquo(Dim_b_column)
		
		# Find redundant elements
		.data.redundant =
			
			# Calculate distances
			.data %>%
			select(!!.element,!!Dim_a_column,!!Dim_b_column) %>%
			distinct() %>%
			as_matrix(rownames = !!.element) %>%
			dist() %>%
			
			# Prepare matrix
			as.matrix() %>% as_tibble(rownames = "element_a") %>%
			gather(`element_b`, dist,-`element_a`) %>%
			filter(`element_a` != `element_b`) %>%
			
			# Sort the elements of the two columns to avoid eliminating all elements
			rowwise() %>%
			mutate(
				`element_1` = c(`element_a`, `element_b`) %>% sort() %>% `[`(1),
				`element_2` = c(`element_a`, `element_b`) %>% sort() %>% `[`(2)
			) %>%
			ungroup() %>%
			select(`element_1`, `element_2`, dist) %>%
			distinct() %>%
			
			# Select closestpairs
			select_closest_pairs %>%
			
			# Select pair to keep
			select(1) %>%
			
			# Set the column names
			setNames(quo_names(.element))
		
		# Drop elements that are correlated with others and return
		.data %>% anti_join(.data.redundant) %>%
			
			# Attach attributes
			reattach_internals(.data)
	}


#' This function is needed for DE in case the matrix is not rectangular, but includes NA
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, of the kind ~ factor_of_intrest + batch
#' @param .element The name of the element column
#' @param .feature The name of the feature/gene column
#' @param .value The name of the feature/gene value column
#' @param .value_scaled The name of the feature/gene scaled value column
#'
#'
#' @return A tibble with adjusted counts
#'
#'
fill_NA_using_formula = function(.data,
																 .formula,
																 .element = NULL,
																 .feature = NULL,
																 .value = NULL,
																 .value_scaled = NULL){
	
	# Comply with CRAN NOTES
	data = NULL
	. = NULL

	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	.value_scaled = enquo(.value_scaled)
	
	# Check that the covariate are unique to elements
	if(
		.data %>%
		select(!!.element, parse_formula(.formula)) %>%
		distinct() %>%
		nrow %>% `>` (.data %>% select(!!.element) %>% distinct() %>% nrow )
	) stop("nanny says: your covariate are not unique to elements. There are elements with multiple covariate values")
	
	# Parse formula
	df_formula =
		.data %>%
		select(parse_formula(.formula)) %>%
		distinct() 
	
	# Check that the at least one covariate is is.character(x) | is.logical(x) | is.factor(x)
	if(
		df_formula %>% lapply(class) %>% unlist %>% intersect(c("character", "logical", "factor")) %>% length %>% equals(0) & 
		length(parse_formula(.formula)) > 0
	) stop("nanny says: none of your covariate are type character, logical, factor, which is needed for element grouping for imputing missing values from formula")
	
	col_formula =
		df_formula %>%
		select_if(function(x) is.character(x) | is.logical(x) | is.factor(x)) %>%
		colnames
	
	# Create NAs for missing element/feature pair
	df_to_impute =
		.data %>%
		select(!!.element, !!.feature, !!.value, col_formula) %>%
		distinct %>%
		pivot_wider(names_from = !!.feature, values_from = !!.value, names_sep = "___") %>%
		pivot_longer(names_to = quo_names(.feature),values_to = quo_names(.value), cols = -c( !!.element, col_formula),	names_sep = when(quo_names(.feature), length(.) > 1 ~ "___", ~ NULL), 
) 
	
	# Select just features/covariates that have missing
	combo_to_impute = df_to_impute %>% anti_join(.data, by=c(quo_names(.element), quo_names(.feature))) %>% select(!!.feature, col_formula) %>% distinct()
	
	# Impute using median
	df_to_impute %>%
		inner_join(combo_to_impute, by=c(quo_names(.feature), col_formula)) %>%
		
		# Calculate median for NAs
		nest(data = -c(col_formula, !!.feature)) %>%
		mutate(data = map(data, ~
												.x %>%
												mutate(
													!!.value := ifelse(
														!!.value %>% is.na,
														median(!!.value, na.rm = TRUE),!!.value
													)
												) %>%
												
												# Impute scaled if exist
												ifelse_pipe(
													quo_is_symbol(.value_scaled),
													~ .x %>% mutate(
														!!.value_scaled := ifelse(
															!!.value_scaled %>% is.na,
															median(!!.value_scaled, na.rm = TRUE),!!.value_scaled
														)
													)
												) %>%
												
												# Throu warning if group of size 1
												ifelse_pipe((.) %>% nrow %>% `<` (2), warning("nanny says: According to your design matrix, u have element groups of size < 2, so you your dataset could still be sparse."))
		)) %>%
		unnest(data) %>%
		
		# Select only imputed data
		select(-col_formula) %>%
		
		# In next command avoid error if no data to impute
		ifelse_pipe(
			nrow(.) > 0,
			~ .x %>% left_join(.data %>% subset(!!.element), by=quo_names(.element))
		) %>%
		
		# Add oiginal dataset
		bind_rows(.data %>% anti_join(combo_to_impute, by=c(quo_names(.feature), col_formula))) %>%
		select(.data %>% colnames)
	
}

#' This function is needed for DE in case the matrix is not rectangular, but includes NA
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .element The name of the element column
#' @param .feature The name of the feature/gene column
#' @param .value The name of the feature/gene value column
#' @param fill_with A numerical value with which fill the mssing data points
#'
#'
#' @return A tibble with adjusted counts
#'
#'
fill_NA_using_value = function(.data,
																 .element = NULL,
																 .feature = NULL,
																 .value = NULL,
																 fill_with){
	
	# Comply with CRAN NOTES
	. = NULL
	 
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	
	# Create NAs for missing element/feature pair
	df_to_impute =
		.data %>%
		select(!!.element, !!.feature, !!.value) %>%
		distinct %>%
		pivot_wider(
			names_from = !!.feature,
			values_from = !!.value,
			names_sep = "___", 
			names_prefix = "fill_miss_"
		) %>%
		pivot_longer(
			names_to = .data %>% select(!!.feature) %>% names, 
			values_to = quo_names(.value), 
			names_sep = purrr::when(quo_names(.feature), length(.) > 1 ~ "___", ~ NULL), 
			names_prefix = "fill_miss_", 
			cols = contains("fill_miss_")
		)
		
	# Select just features/covariates that have missing
	combo_to_impute = df_to_impute %>% anti_join(.data, by=c(quo_names(.element), quo_names(.feature))) %>% select(!!.feature, !!.element) %>% distinct()
	
	# Impute using median
	df_to_impute %>%
		inner_join(combo_to_impute) %>%
		
		# Fill
		mutate(!!.value := ifelse(!!.value %>% is.na, fill_with, !!.value)) %>%
	
		# In next command avoid error if no data to impute
		ifelse_pipe(
			nrow(.) > 0,
			~ .x %>% left_join(.data %>% subset(!!.element), by=quo_names(.element))
		) %>%
		
		# Add oiginal dataset
		bind_rows(.data %>% anti_join(combo_to_impute, by=c(quo_names(.feature), quo_names(.element)))) %>%
		select(.data %>% colnames)
	
}




