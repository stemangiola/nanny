% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{subset}
\alias{subset}
\alias{subset.spec_tbl_df}
\alias{subset.tbl_df}
\alias{subset.tbl}
\title{Extract selected-column-wise information}
\usage{
subset(.data, .column)

\method{subset}{spec_tbl_df}(.data, .column)

\method{subset}{tbl_df}(.data, .column)

\method{subset}{tbl}(.data, .column)
}
\arguments{
\item{.data}{A `tbl`}

\item{.column}{The name of the column of interest}
}
\value{
A `tbl` object

A `tbl` object

A `tbl` object

A `tbl` object
}
\description{
subset() takes as input a `tbl` and returns a `tbl` with only selected-column-related columns
}
\details{
\lifecycle{maturing}

This functon extracts only selected-column-related information for downstream analysis (e.g., visualisation). It is disruptive in the sense that it cannot be passed anymore to nanny function.
}
\examples{

subset(mtcars_tidy,car_model)


}
