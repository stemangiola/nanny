add_get_only_error_message = 	"nanny says: action must be either \"add\" for adding this information to your data frame, \"get\" to get the new information along with the information element-wise of feature-wide (depending on what operation you are doing) or only to just get the new information"
add_get_error_message = 	"nanny says: action must be either \"add\" for adding this information to your data frame, \"get\" to get the new information along with the information element-wise of feature-wide (depending on what operation you are doing)"

my_stop = function() {
	stop("
        You should call nanny library *after* tidyverse libraries.
        nanny says: The function does not know what your element, feature and counts columns are.
        You have to either enter those as arguments, or use the funtion nanny() to pass your column names that will be remembered.
      ")
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @importFrom purrr as_mapper
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)
	
}

#' This is a generalisation of ifelse that accepts an object and return an objects
#' 
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#'
#' @param .x A tibble
#' @param .p1 A boolean
#' @param .p2 ELSE IF condition
#' @param .f1 A function
#' @param .f2 A function
#' @param .f3 A function
#'
#' @return A tibble
ifelse2_pipe = function(.x, .p1, .p2, .f1, .f2, .f3 = NULL) {
	# Nested switch
	switch(# First condition
		.p1 %>% `!` %>% sum(1),
		
		# First outcome
		as_mapper(.f1)(.x),
		switch(
			# Second condition
			.p2 %>% `!` %>% sum(1),
			
			# Second outcome
			as_mapper(.f2)(.x),
			
			# Third outcome - if there is not .f3 just return the original data frame
			if (.f3 %>% is.null %>% `!`)
				as_mapper(.f3)(.x)
			else
				.x
		))
}

#' Sub function of remove_redundancy_elements_though_reduced_dimensions
#' 
#' @keywords internal
#'
#' @importFrom stats dist
#' @importFrom utils head
#'
#' @param df A tibble
#'
#'
#' @return A tibble with pairs to drop
select_closest_pairs = function(df) {
	
	# Comply with CRAN NOTES
	. = NULL
	element_1 = NULL
	element_2 = NULL
	
	couples <- df %>% head(n = 0)
	
	while (df %>% nrow() > 0) {
		pair <- df %>%
			arrange(dist) %>%
			head(n = 1)
		couples <- couples %>% bind_rows(pair)
		df <- df %>%
			filter(
				!`element_1` %in% (pair %>% select(1:2) %>% as.character()) &
					!`element_2` %in% (pair %>% select(1:2) %>% as.character())
			)
	}
	
	couples
	
}

#' get_x_y_annotation_columns
#' 
#' @keywords internal
#'
#' @importFrom magrittr equals
#'
#' @param .data A `tbl` formatted as | <element> | <feature> | <value> | <...> |
#' @param .horizontal The name of the column horizontally presented in the heatmap
#' @param .vertical The name of the column vertically presented in the heatmap
#' @param .value The name of the feature/gene value column
#'
#' @description This function recognise what are the element-wise columns and transcrip-wise columns
#'
#' @return A list
#'
get_x_y_annotation_columns = function(.data, .horizontal, .vertical, .value){
	
	
	# Comply with CRAN NOTES
	. = NULL
	
	# Make col names
	.horizontal = enquo(.horizontal)
	.vertical = enquo(.vertical)
	.value = enquo(.value)

	# x-annotation df
	n_x = .data %>% select(!!.horizontal) %>% distinct() %>% nrow
	n_y = .data %>% select(!!.vertical) %>% distinct() %>% nrow
	
	# element wise columns
	horizontal_cols=
		.data %>%
		select(-!!.horizontal, -!!.vertical, -!!.value) %>%
		colnames %>%
		map(
			~
				.x %>%
				ifelse_pipe(
					.data %>%
						distinct(!!.horizontal, !!as.symbol(.x)) %>%
						nrow %>%
						equals(n_x),
					~ .x,
					~ NULL
				)
		) %>%
		
		# Drop NULL
		{	(.)[lengths((.)) != 0]	} %>%
		unlist
	
	# feature wise columns
	vertical_cols=
		.data %>%
		select(-!!.horizontal, -!!.vertical, -!!.value, -horizontal_cols) %>%
		colnames %>%
		map(
			~
				.x %>%
				ifelse_pipe(
					.data %>%
						distinct(!!.vertical, !!as.symbol(.x)) %>%
						nrow %>%
						equals(n_y),
					~ .x,
					~ NULL
				)
		) %>%
		
		# Drop NULL
		{	(.)[lengths((.)) != 0]	} %>%
		unlist
	
	# Counts wise columns, at the moment scaled counts is treated as special and not accounted for here
	counts_cols =
		.data %>%
		select(-!!.horizontal, -!!.vertical, -!!.value) %>%
		
		# Exclude horizontal
		ifelse_pipe(!is.null(horizontal_cols),  ~ .x %>% select(-horizontal_cols)) %>%
		
		# Exclude vertical
		ifelse_pipe(!is.null(vertical_cols),  ~ .x %>% select(-vertical_cols)) %>%
		
		# Select colnames
		colnames %>%
		
		# select columns
		map(
			~
				.x %>%
				ifelse_pipe(
					.data %>%
						distinct(!!.vertical, !!.horizontal, !!as.symbol(.x)) %>%
						nrow %>%
						equals(n_x * n_y),
					~ .x,
					~ NULL
				)
		) %>%
		
		# Drop NULL
		{	(.)[lengths((.)) != 0]	} %>%
		unlist
	
	list(  horizontal_cols = horizontal_cols,  vertical_cols = vertical_cols, counts_cols = counts_cols )
}

get_specific_annotation_columns = function(.data, .col){
	
	 
	# Comply with CRAN NOTES
	. = NULL
	
	# Make col names
	.col = enquo(.col)
	
	# x-annotation df
	n_x = .data %>% distinct_at(vars(!!.col)) %>% nrow
	
	# element wise columns
	.data %>%
		select(-!!.col) %>%
		colnames %>%
		map(
			~
				.x %>%
				ifelse_pipe(
					.data %>%
						distinct_at(vars(!!.col, .x)) %>%
						nrow %>%
						equals(n_x),
					~ .x,
					~ NULL
				)
		) %>%
		
		# Drop NULL
		{	(.)[lengths((.)) != 0]	} %>%
		unlist
	
}

initialise_internals = function(.data){
	
	# Comply with CRAN NOTES
	. = NULL
	
	.data %>%
		ifelse_pipe(
			"internals" %in% ((.) %>% attributes %>% names) %>% `!`,
			~ .x %>% add_attr(list(), "internals")
		)
}

reattach_internals = function(.data, .data_internals_from = NULL){
	if(.data_internals_from %>% is.null)
		.data_internals_from = .data
	
	.data %>% add_attr(.data_internals_from %>% attr("internals"), "internals")
}

attach_to_internals = function(.data, .object, .name){
	
	internals =
		.data %>%
		initialise_internals() %>%
		attr("internals")
	
	# Add tt_bolumns
	internals[[.name]] = .object
	
	.data %>% add_attr(internals, "internals")
}

drop_internals = function(.data){
	
	.data %>% drop_attr("internals")
}

#' Add attribute to abject
#' 
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
	attr(var, name) <- attribute
	var
}

#' Drop attribute to abject
#' 
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
drop_attr = function(var, name) {
	attr(var, name) <- NULL
	var
}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#' 
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {

	v = quo_name(quo_squash(v))
	gsub('^c\\(|`|\\)$', '', v) %>% 
		strsplit(', ') %>% 
		unlist 
}

# Function that rotates a 2D space of a arbitrary angle
rotation = function(m, d) {
	r = d * pi / 180
	((dplyr::bind_rows(
		c(`1` = cos(r), `2` = -sin(r)),
		c(`1` = sin(r), `2` = cos(r))
	) %>% as_matrix) %*% m)
}

#' .formula parser
#' 
#' @keywords internal
#'
#' @importFrom stats terms
#'
#' @param fm a formula
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
	if (attr(terms(fm), "response") == 1)
		stop("nanny says: The .formula must be of the kind \"~ covariates\" ")
	else
		as.character(attr(terms(fm), "variables"))[-1]
}

#' Remove class to abject
#' 
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
drop_class = function(var, name) {
	class(var) <- class(var)[!class(var)%in%name]
	var
}


nanny_to_tbl = function(.data) {
	.data %>%	drop_class(c("nanny", "tt"))
}

# From tidyr
strip_names <- function(df, base, names_sep) {
	base <- paste0(base, names_sep)
	names <- names(df)
	
	has_prefix <- regexpr(base, names, fixed = TRUE) == 1L
	names[has_prefix] <- substr(names[has_prefix], nchar(base) + 1, nchar(names[has_prefix]))
	
	set_names(df, names)
}

# Raise to the power
pow = function(a,b){	a^b }