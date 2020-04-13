#' Get K-mean clusters to a tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom stats kmeans
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally elements)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_elements A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
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
					 log_transform = TRUE,
					 ...) {
		# Check if centers is in dots
		dots_args = rlang::dots_list(...)
		if ("centers" %in% names(dots_args) %>% `!`)
			stop("nanny says: for kmeans you need to provide the \"centers\" integer argument")
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		.data %>%
			
			# Through error if some counts are NA
			error_if_counts_is_na(!!.value) %>%
			
			# Prepare data frame
			distinct(!!.feature,!!.element,!!.value) %>%
			
			# Check if log tranfrom is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.value := !!.value %>%  `+`(1) %>%  log())) %>%
			
			# Prepare data frame for return
			spread(!!.feature,!!.value) %>%
			as_matrix(rownames = !!.element) %>%
			
			# Wrap the do.call because of the centers check
			{
				do.call(kmeans, list(x = (.), iter.max = 1000) %>% c(dots_args))
			}	 %$%
			cluster %>%
			as.list() %>%
			as_tibble() %>%
			gather(!!.element, `cluster kmeans`) %>%
			mutate(`cluster kmeans` = `cluster kmeans` %>% as.factor()) %>%
			
			# Attach attributes
			reattach_internals(.data)
	}

#' Get SNN shared nearest neighbour clusters to a tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally elements)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally genes)
#' @param of_elements A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
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
					 log_transform = TRUE,
					 ...) {
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		# Check if package is installed, otherwise install
		if ("Seurat" %in% rownames(installed.packages()) == FALSE) {
			writeLines("Installing Seurat")
			install.packages("Seurat", repos = "https://cloud.r-project.org")
		}
		if ("KernSmooth" %in% rownames(installed.packages()) == FALSE) {
			writeLines("Installing KernSmooth")
			install.packages("KernSmooth", repos = "https://cloud.r-project.org")
		}
		
		my_df =
			.data %>%
			
			# Through error if some counts are NA
			error_if_counts_is_na(!!.value) %>%
			
			# Prepare data frame
			distinct(!!.element,!!.feature,!!.value) %>%
			
			# Check if log tranfrom is needed
			#ifelse_pipe(log_transform, ~ .x %>% dplyr::mutate(!!.value := !!.value %>%  `+`(1) %>%  log())) %>%
			
			# Prepare data frame for return
			spread(!!.element,!!.value)
		
		my_df %>%
			data.frame(row.names = quo_name(.feature)) %>%
			Seurat::CreateSeuratObject() %>%
			Seurat::ScaleData(display.progress = TRUE,
												num.cores = 4,
												do.par = TRUE) %>%
			Seurat::FindVariableFeatures(selection.method = "vst") %>%
			Seurat::RunPCA(npcs = 30) %>%
			Seurat::FindNeighbors() %>%
			Seurat::FindClusters(method = "igraph", ...) %>%
			`[[` ("seurat_clusters") %>%
			as_tibble(rownames = quo_name(.element)) %>%
			rename(`cluster SNN` = seurat_clusters) %>%
			dplyr::mutate(!!.element := gsub("\\.", "-",!!.element)) %>%
			
			# Attach attributes
			reattach_internals(.data)
	}

#' Get dimensionality information to a tibble using MDS
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang :=
#' @importFrom stats setNames
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param top An integer. How many top genes to select
#' @param of_elements A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
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
					 top = 500,
					 of_elements = TRUE,
					 log_transform = TRUE) {
		# Comply with CRAN NOTES
		. = NULL
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		# Get components from dims
		components = 1:.dims
		
		mds_object =
			.data %>%
			
			# Through error if some counts are NA
			error_if_counts_is_na(!!.value) %>%
			
			# Filter lowly transcribed (I have to avoid the use of scaling function)
			keep_abundant(!!.element, !!.feature,!!.value) %>%
			distinct(!!.feature,!!.element,!!.value) %>%
			
			# Check if logtansform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.value := !!.value %>% `+`(1) %>%  log())) %>%
			
			# Stop any column is not if not numeric or integer
			ifelse_pipe(
				(.) %>% select(!!.value) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
				~ stop("nanny says: .value must be numerical or integer")
			) %>%
			spread(!!.element,!!.value) %>%
			as_matrix(rownames = !!.feature, do_check = FALSE) %>%
			limma::plotMDS(ndim = .dims, plot = FALSE, top = top)
		
		# Pase results
		mds_object %$%	cmdscale.out %>%
			as.data.frame %>%
			as_tibble(rownames = quo_name(.element)) %>%
			setNames(c(quo_name(.element), sprintf("Dim%s", 1:.dims))) %>%
			
			
			# Attach attributes
			reattach_internals(.data) %>%
			
			# Add raw object
			attach_to_internals(mds_object, "MDS") %>%
			# Communicate the attribute added
			{
				message("nanny says: to access the raw results do `attr(..., \"tt_internals\")$MDS`")
				(.)
			}
		
	}

#' Get principal component information to a tibble using PCA
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats prcomp
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param top An integer. How many top genes to select
#' @param of_elements A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
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
					 top = 500,
					 of_elements = TRUE,
					 log_transform = TRUE,
					 scale = FALSE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL
		
		# Get column names
		.element = enquo(.element)
		.feature = enquo(.feature)
		.value = enquo(.value)
		
		# Get components from dims
		components = 1:.dims
		
		prcomp_obj =
			.data %>%
			
			# Through error if some counts are NA
			error_if_counts_is_na(!!.value) %>%
			
			# Prepare data frame
			distinct(!!.feature,!!.element,!!.value) %>%
			
			# Check if logtansform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.value := !!.value %>% `+`(1) %>%  log())) %>%
			
			# Stop any column is not if not numeric or integer
			ifelse_pipe(
				(.) %>% select(!!.value) %>% summarise_all(class) %>% `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
				~ stop("nanny says: .value must be numerical or integer")
			) %>%
			
			# Filter most variable genes
			keep_variable_features(!!.element,!!.feature,!!.value, top) %>%
			
			spread(!!.element,!!.value) %>%
			
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
				writeLines("Fraction of variance explained by the selected principal components")
				
				(.) %$% sdev %>% `^` (2) %>% # Eigen value
					`/` (sum(.)) %>%
					`[` (components) %>%
					enframe() %>%
					select(-name) %>%
					rename(`Fraction of variance` = value) %>%
					mutate(PC = components) %>%
					print(n = 9999999)
				
				(.)
				
			} %$%
			
			# Parse the PCA results to a tibble
			rotation %>%
			as_tibble(rownames = quo_name(.element)) %>%
			select(!!.element, sprintf("PC%s", components)) %>%
			
			# Attach attributes
			reattach_internals(.data) %>%
			
			# Add raw object
			attach_to_internals(prcomp_obj, "PCA") %>%
			# Communicate the attribute added
			{
				message("nanny says: to access the raw results do `attr(..., \"tt_internals\")$PCA`")
				(.)
			}
		
	}

#' Get principal component information to a tibble using tSNE
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param top An integer. How many top genes to select
#' @param of_elements A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
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
					 top = 500,
					 of_elements = TRUE,
					 log_transform = TRUE,
					 ...) {
		# Comply with CRAN NOTES
		. = NULL
		
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
		
		
		# Check if package is installed, otherwise install
		if ("Rtsne" %in% rownames(installed.packages()) == FALSE) {
			writeLines("Installing Rtsne")
			install.packages("Rtsne", repos = "https://cloud.r-project.org")
		}
		
		# Set perprexity to not be too high
		if (!"perplexity" %in% names(arguments))
			arguments = arguments %>% c(perplexity = ((
				.data %>% distinct(!!.element) %>% nrow %>% sum(-1)
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
			distinct(!!.feature,!!.element,!!.value) %>%
			
			# Check if data rectangular
			ifelse_pipe(
				(.) %>% check_if_data_rectangular(!!.element,!!.feature,!!.value, type = "soft"),
				~ .x %>% eliminate_sparse_features(!!.feature)
			) %>%
			
			# Check if logtansform is needed
			ifelse_pipe(log_transform,
									~ .x %>% dplyr::mutate(!!.value := !!.value %>% `+`(1) %>%  log())) %>%
			
			# Filter most variable genes
			keep_variable_features(!!.element,!!.feature,!!.value, top) %>%
			
			spread(!!.feature,!!.value) %>%
			# select(-element) %>%
			# distinct %>%
			as_matrix(rownames = quo_name(.element))
		
		do.call(Rtsne::Rtsne, c(list(df_tsne), arguments)) %$%
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
		# Get column names
		.element = enquo(.element)
		dimension_1_column = enquo(dimension_1_column)
		dimension_2_column = enquo(dimension_2_column)
		dimension_1_column_rotated = enquo(dimension_1_column_rotated)
		dimension_2_column_rotated = enquo(dimension_2_column_rotated)
		
		if (.data %>%
				distinct(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
				count(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
				pull(n) %>%
				max %>%
				`>` (1))
			stop(sprintf(
				"nanny says: %s must be unique for each row for the calculation of rotation",
				quo_name(.element)
			))
		
		# Function that rotates a 2D space of a arbitrary angle
		rotation = function(m, d) {
			r = d * pi / 180
			((dplyr::bind_rows(
				c(`1` = cos(r), `2` = -sin(r)),
				c(`1` = sin(r), `2` = cos(r))
			) %>% as_matrix) %*% m)
		}
		
		# Sanity check of the angle selected
		if (rotation_degrees %>% between(-360, 360) %>% `!`)
			stop("nanny says: rotation_degrees must be between -360 and 360")
		
		# Return
		.data %>%
			distinct(!!.element, !!dimension_1_column, !!dimension_2_column) %>%
			as_matrix(rownames = !!.element) %>% t %>%
			rotation(rotation_degrees) %>%
			as_tibble() %>%
			mutate(`rotated dimensions` =
						 	c(
						 		quo_name(dimension_1_column_rotated),
						 		quo_name(dimension_2_column_rotated)
						 	)) %>%
			gather(!!.element, value,-`rotated dimensions`) %>%
			spread(`rotated dimensions`, value) %>%
			
			# Attach attributes
			reattach_internals(.data)
		
	}

#' Drop redundant elements (e.g., elements) for which feature (e.g., genes) aboundances are correlated
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#'
#' @param .data A tibble
#' @param .value A column symbol with the value the clustering is based on (e.g., `count`)
#' @param correlation_threshold A real number between 0 and 1
#' @param top An integer. How many top genes to select
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally elements)
#' @param of_elements A boolean
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
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
																													 log_transform = FALSE) {
	# Comply with CRAN NOTES
	. = NULL
	
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	col_names = get_elements_features_value(.data, .element, .feature, .value, of_elements)
	.element = col_names$.element
	.feature = col_names$.feature
	.value = col_names$.value
	
	# Check if package is installed, otherwise install
	if ("widyr" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing widyr needed for correlation analyses")
		install.packages("widyr")
	}
	
	# Get the redundant data frame
	.data.correlated =
		.data %>%
		
		# Stop if any counts is NA
		error_if_counts_is_na(!!.value) %>%
		
		# Stop if there are duplicated features
		error_if_duplicated_genes(!!.element,!!.feature,!!.value) %>%
		
		# Prepare the data frame
		select(!!.feature,!!.element,!!.value) %>%
		
		# Filter variable genes
		keep_variable_features(!!.element,!!.feature,!!.value, top = top) %>%
		
		# Check if logtansform is needed
		ifelse_pipe(log_transform,
								~ .x %>% dplyr::mutate(!!.value := !!.value %>% `+`(1) %>%  log())) %>%
		distinct() %>%
		spread(!!.element,!!.value) %>%
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
		gather(!!.element,!!.value,-!!.feature) %>%
		dplyr::rename(rc := !!.value,
									element := !!.element,
									feature := !!.feature) %>% # Is rename necessary?
		mutate_if(is.factor, as.character) %>%
		
		# Run pairwise correlation and return a tibble
		widyr::pairwise_cor(
			element,
			feature,
			rc,
			sort = TRUE,
			diag = FALSE,
			upper = FALSE
		) %>%
		filter(correlation > correlation_threshold) %>%
		distinct(item1) %>%
		dplyr::rename(!!.element := item1)
	
	# Return non redudant data frame
	.data %>% anti_join(.data.correlated) %>%
		
		# Attach attributes
		reattach_internals(.data)
}

#' Identifies the closest pairs in a MDS contaxt and return one of them
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
		# This function identifies the closest pairs and return one of them
		
		# Get column names
		.element = enquo(.element)
		col_names = get_elements(.data, .element)
		.element = col_names$.element
		
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
			as.matrix() %>% as_tibble(rownames = "element a") %>%
			gather(`element b`, dist,-`element a`) %>%
			filter(`element a` != `element b`) %>%
			
			# Sort the elements of the two columns to avoid eliminating all elements
			rowwise() %>%
			mutate(
				`element 1` = c(`element a`, `element b`) %>% sort() %>% `[`(1),
				`element 2` = c(`element a`, `element b`) %>% sort() %>% `[`(2)
			) %>%
			ungroup() %>%
			select(`element 1`, `element 2`, dist) %>%
			distinct() %>%
			
			# Select closestpairs
			select_closest_pairs %>%
			
			# Select pair to keep
			select(1) %>%
			
			# Set the column names
			setNames(quo_name(.element))
		
		# Drop elements that are correlated with others and return
		.data %>% anti_join(.data.redundant) %>%
			
			# Attach attributes
			reattach_internals(.data)
	}

#' Get matrix from tibble
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' as_matrix(head(dplyr::select(nanny::counts_mini, feature, count)), rownames=feature)
#'
#' @export
as_matrix <- function(tbl,
											rownames = NULL,
											do_check = TRUE) {
	rownames = enquo(rownames)
	tbl %>%
		
		# Through warning if data frame is not numerical beside the rownames column (if present)
		ifelse_pipe(
			do_check &&
				tbl %>%
				# If rownames defined eliminate it from the data frame
				ifelse_pipe(!quo_is_null(rownames), ~ .x[,-1], ~ .x) %>%
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
		as.data.frame() %>%
		
		# Deal with rownames column if present
		ifelse_pipe(
			!quo_is_null(rownames),
			~ .x %>%
				magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
				select(-1)
		) %>%
		
		# Convert to matrix
		as.matrix()
}

#' This function is needed for DE in case the matrix is not rectangular, but includes NA
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @importFrom utils installed.packages
#' @importFrom utils install.packages
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
	
	# Get column names
	.element = enquo(.element)
	.feature = enquo(.feature)
	.value = enquo(.value)
	.value_scaled = enquo(.value_scaled)
	
	col_formula =
		.data %>%
		select(parse_formula(.formula)) %>%
		distinct() %>%
		select_if(function(x) is.character(x) | is.logical(x) | is.factor(x)) %>%
		colnames
	
	# Create NAs for missing element/feature pair
	df_to_impute =
		.data %>%
		select(!!.element, !!.feature, !!.value, col_formula) %>%
		distinct %>%
		spread(!!.feature, !!.value) %>%
		gather(!!.feature, !!.value, -!!.element, -col_formula)
	
	# Select just features/covariates that have missing
	combo_to_impute = df_to_impute %>% anti_join(.data, by=c(quo_name(.element), quo_name(.feature))) %>% select(!!.feature, col_formula) %>% distinct()
	
	# Impute using median
	df_to_impute %>%
		inner_join(combo_to_impute, by=c(quo_name(.feature), col_formula)) %>%
		
		# Calculate median for NAs
		nest(data = -c(col_formula, !!.feature)) %>%
		mutate(data = map(data, ~
												.x %>%
												mutate(
													!!.value := ifelse(
														!!.value %>% is.na,
														median(!!.value, na.rm = T),!!.value
													)
												) %>%
												
												# Impute scaled if exist
												ifelse_pipe(
													quo_is_symbol(.value_scaled),
													~ .x %>% mutate(
														!!.value_scaled := ifelse(
															!!.value_scaled %>% is.na,
															median(!!.value_scaled, na.rm = T),!!.value_scaled
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
			~ .x %>% left_join(.data %>% pivot_element(!!.element), by=quo_name(.element))
		) %>%
		
		# Add oiginal dataset
		bind_rows(.data %>% anti_join(combo_to_impute, by=c(quo_name(.feature), col_formula))) %>%
		select(.data %>% colnames)
	
}

permute_nest = function(.data, .names_from, .values_from){
	.names_from = enquo(.names_from)
	.values_from = enquo(.values_from)
	
	factor_levels = .data %>% pull(!!.names_from) %>% unique
	
	.data %>% 
		pull(!!.names_from) %>%
		unique() %>%
		gtools::permutations(n = length(.), r = 2, v = .) %>%
		as_tibble() %>%
		unite(run, c(V1, V2), remove = F, sep="___") %>%
		gather(which, !!.names_from, -run) %>%
		select(-which) %>%
		left_join(.data %>% select(!!.names_from, !!.values_from), by = quo_name(.names_from)) %>%
		nest(data = -run) %>%
		separate(run, sprintf("%s_%s", quo_name(.names_from), 1:2 ), sep="___") %>%
		
		# Introduce levels
		mutate_at(vars(1:2),function(x) factor(x, levels = factor_levels))
	
}

combine_nest = function(.data, .names_from, .values_from){
	.names_from = enquo(.names_from)
	.values_from = enquo(.values_from)
	
	factor_levels = .data %>% pull(!!.names_from) %>% unique
	
	.data %>% 
		pull(!!.names_from) %>%
		unique() %>%
		gtools::combinations(n = length(.), r = 2, v = .) %>%
		as_tibble() %>%
		unite(run, c(V1, V2), remove = F, sep="___") %>%
		gather(which, !!.names_from, -run) %>%
		select(-which) %>%
		left_join(.data %>% select(!!.names_from, !!.values_from), by = quo_name(.names_from)) %>%
		nest(data = -run) %>%
		separate(run, sprintf("%s_%s", quo_name(.names_from), 1:2), sep="___") %>%
		
		# Introduce levels
		mutate_at(vars(1:2),function(x) factor(x, levels = factor_levels))
	
}
