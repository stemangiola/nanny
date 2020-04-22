#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param .data A tibble of read counts
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
#'
error_if_counts_is_na = function(.data, .abundance) {
	.abundance = enquo(.abundance)
	
	# Do the check
	if (.data %>% filter(!!.abundance %>% is.na) %>% nrow %>% `>` (0))
		stop("nanny says: You have NA values in your counts")
	
	# If all good return original data frame
	.data
}



#' Check whether a numeric vector has been log transformed
#'
#' @param x A numeric vector
#' @param .abundance A character name of the transcript/gene abundance column
#'
#' @return NA
error_if_log_transformed <- function(x, .abundance) {
	.abundance = enquo(.abundance)
	
	if (x %>% nrow %>% `>` (0))
		if (x %>% summarise(m = !!.abundance %>% max) %>% pull(m) < 50)
			stop(
				"nanny says: The input was log transformed, this algorithm requires raw (un-scaled) read counts"
			)
}

#' Check whether there are duplicated genes/transcripts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param .data A tibble of read counts
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
error_if_duplicated_genes <- function(.data,
																			.sample = `sample`,
																			.transcript = `transcript`,
																			.abundance = `read count`) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	duplicates <-
		distinct( .data, !!.sample,!!.transcript,!!.abundance) %>%
		count(!!.sample,!!.transcript) %>%
		filter(n > 1) %>%
		arrange(n %>% desc())
	
	if (duplicates %>% nrow() > 0) {
		writeLines("Those are the duplicated genes")
		duplicates %>% print()
		stop(
			"nanny says: Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding."
		)
	}
	
	.data
	
}



#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map
#'
#'
#' @param .data A tibble of read counts
#' @param list_input A list
#' @param expected_type A character string
#'
#' @return A tbl
#'
error_if_wrong_input = function(.data, list_input, expected_type) {
	
	
	
	
	# Do the check
	if (
		list_input %>%
		map(~ .x %>% class() %>% `[` (1)) %>%
		unlist %>%
		equals(expected_type) %>%
		`!`
	)
		stop("nanny says: You have passed the wrong argument to the function. Please check again.")
	
	# If all good return original data frame
	.data
}

#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map
#'
#'
#' @param .data A tibble of read counts
#' @param list_input A list
#' @param expected_type A character string
#'
#' @return A tbl
#'
check_if_wrong_input = function(.data, list_input, expected_type) {
	# Do the check
	if (list_input %>%
			map( ~ .x %>% class() %>% `[` (1)) %>%
			unlist %>%
			equals(expected_type) %>%
			`!`)
		stop("nanny says: You have passed the wrong argument to the function. Please check again.")
	
	# If all good return original data frame
	.data
}

#' Check whether there are duplicated genes/transcripts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param .data A tibble of read counts
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
check_if_duplicated_genes <- function(.data,
																			.sample = `sample`,
																			.transcript = `transcript`,
																			.abundance = `read count`) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	duplicates <-
		select(.data,!!.sample, !!.transcript, !!.abundance) %>%
		distinct() %>%
		group_by_at(vars(!!.sample, !!.transcript)) %>%
		tally() %>%
		filter(n > 1) %>%
		arrange(n %>% desc())
	
	is_unique = duplicates %>% nrow() %>% equals(0)
	
	if (!is_unique) {
		writeLines("nanny says: Those are the duplicated genes")
		duplicates %>% print()
	}
	
	is_unique
}

#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param .data A tibble of read counts
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
#'
check_if_counts_is_na = function(.data, .abundance) {
	.abundance = enquo(.abundance)
	
	.data %>% filter(!!.abundance %>% is.na) %>% nrow %>% equals(0)
	
}

check_if_column_missing = function(.data, .sample, .transcript, .abundance) {
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	# Check that the intersection is length good
	my_cols = c(
		quo_names(.sample),
		quo_names(.transcript),
		quo_names(.abundance)
	)
	
	.data %>% colnames %>%
		intersect(my_cols) %>%
		length %>%
		equals(length(my_cols))
}

column_type_checking = function(.data, .sample, .transcript, .abundance) {
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	.data %>% select(!!.sample) %>% map_chr(~ class(.x)) %in% c("character", "factor") %>% all &
	.data %>% select(!!.transcript) %>% map_chr(~ class(.x)) %in% c("character", "factor") %>% all &
		.data %>% select(!!.abundance) %>% map_chr(~ class(.x)) %in% c("integer", "numeric", "double") %>% all
	
}

check_if_attribute_present = function(.data) {
	"tt_internals" %in% (.data %>% attributes %>% names) &&
		"tt_columns" %in% (.data %>% attr("tt_internals")  %>% names)
}

eliminate_sparse_transcripts = function(.data, .transcript){
	# Parse column names
	.transcript = enquo(.transcript)
	
	warning("nanny says: Some transcripts have been omitted from the analysis because not present in every sample.")
	
	.data %>%
		add_count(!!.transcript, name = "my_n") %>%
		filter(my_n == max(my_n)) %>%
		select(-my_n)
}

check_if_data_rectangular = function(.data, .sample, .transcript, .abundance, type = "hard"){
	
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	is_rectangular =
		.data %>%
		distinct(!!.sample, !!.transcript, !!.abundance) %>%
		count(!!.sample) %>%
		count(n) %>%
		nrow %>%
		equals(1)
	
	is_rectangular
	
	# if(!is_rectangular & type == "hard") stop("nanny says: the data must have the same number of transcript per sample.")
	#
	# if(!is_rectangular & type == "soft") warning("nanny says: the data should have the same number of transcript per sample.")
	
	
	# # Eliminate sparse transcripts
	# .data %>% eliminate_sparse_transcripts(!!.transcript)
	
	
}

nanny_to_tbl = function(.data) {
	.data %>%	drop_class(c("nanny", "tt"))
}

validation_default = function(.data,
															.sample,
															.transcript,
															.abundance,
															type = "hard",
															skip_dupli_check = FALSE) {
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	# Type check
	is_present = check_if_column_missing(.data,!!.sample,!!.transcript,!!.abundance)
	if (type == "hard" &
			!is_present)
		stop(
			"nanny says: One or more columns .sample .transcript or .abundance are missing from your data frame."
		)
	if (type == "soft" & !is_present) {
		warning(
			"nanny says: One or more columns .sample .transcript or .abundance are missing from your data frame. The nanny object has been converted to a `tbl`"
		)
		return(.data %>% nanny_to_tbl)
	}
	
	# Type check
	is_type = column_type_checking(.data,!!.sample,!!.transcript,!!.abundance)
	if (type == "hard" &
			!is_type)
		stop(
			"nanny says: The column provided as .sample .transcript or .abundance do not comply with the required types (<FACTOR/CHARACTER>, <FACTOR/CHARACTER>, <NUMERIC>)."
		)
	if (type == "soft" & !is_type) {
		warning(
			"nanny says: The column provided as .sample .transcript or .abundance do not comply with the required types. The nanny object has been converted to a `tbl`"
		)
		return(.data %>% nanny_to_tbl)
	}
	
	# Check if duplicated genes
	if (!skip_dupli_check) {
		is_unique = check_if_duplicated_genes(.data,!!.sample,!!.transcript,!!.abundance)
		if (type == "hard" &
				!is_unique)
			stop(
				"nanny says: Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding (e.g., aggregate_duplicates())."
			)
		if (type == "soft" & !is_unique) {
			warning(
				"nanny says: Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding (e.g., aggregate_duplicates()). The nanny object has been converted to a `tbl`"
			)
			return(.data %>% nanny_to_tbl)
		}
	}
	
	# Check if NA in counts
	is_count_good = check_if_counts_is_na(.data,!!.abundance)
	if (type == "hard" &
			!is_count_good)
		stop("nanny says: You have NA values in your counts. Please check your data frame.")
	if (type == "soft" & !is_count_good) {
		warning(
			"nanny says: You have NA values in your counts. The nanny object has been converted to a `tbl`"
		)
		return(.data %>% nanny_to_tbl)
	}
	
}

validation <- function(.data,
											 .sample = NULL,
											 .transcript = NULL,
											 .abundance = NULL,
											 type = "hard",
											 skip_dupli_check = FALSE) {
	UseMethod("validation", .data)
}

validation.default = validation_default

# validation.nanny = function(.data,
# 														 .sample = NULL,
# 														 .transcript = NULL,
# 														 .abundance = NULL,
# 														 type = "hard",
# 														 skip_dupli_check = FALSE) {
# 	# Check if attribute is present
# 	is_attr = check_if_attribute_present(.data)
# 	if (type == "hard" &
# 			!is_attr)
# 		stop(
# 			"nanny says: The object provided has nanny class but no attribute containing the column names. Insert a valid nanny object or provide `.sample`, `.transcript`, `.abundance` column names as arguments "
# 		)
# 	if (type == "soft" & !is_attr) {
# 		warning(
# 			"nanny says: The object provided has nanny class but no attribute containing the column names. The nanny object has been converted to a `tbl`"
# 		)
# 		return(.data %>% nanny_to_tbl)
# 	}
# 
# 	# Get column names
# 	.sample = enquo(.sample)
# 	.transcript = enquo(.transcript)
# 	.abundance = enquo(.abundance)
# 	col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
# 	.sample = col_names$.sample
# 	.transcript = col_names$.transcript
# 	.abundance = col_names$.abundance
# 
# 	validation_default(
# 		.data,
# 		!!.sample,
# 		!!.transcript,
# 		!!.abundance,
# 		type = type,
# 		skip_dupli_check = skip_dupli_check
# 	)
# 
# }
