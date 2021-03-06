% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{reduce_dimensions}
\alias{reduce_dimensions}
\alias{reduce_dimensions.spec_tbl_df}
\alias{reduce_dimensions.tbl_df}
\title{Dimension reduction of the feature value data}
\usage{
reduce_dimensions(
  .data,
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
  ...
)

\method{reduce_dimensions}{spec_tbl_df}(
  .data,
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
  ...
)

\method{reduce_dimensions}{tbl_df}(
  .data,
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
  ...
)
}
\arguments{
\item{.data}{A `tbl` formatted as | <element> | <feature> | <value> | <...> |}

\item{.element}{The name of the element column (normally elements).}

\item{.feature}{The name of the feature column (normally features)}

\item{.value}{The name of the column including the numerical value the clustering is based on (normally feature value)}

\item{method}{A character string. The dimension reduction algorithm to use (PCA, MDS, tSNE).}

\item{.dims}{A list of integer vectors corresponding to principal components of interest (e.g., list(1:2, 3:4, 5:6))}

\item{top}{An integer. How many top genes to select for dimensionality reduction}

\item{of_elements}{A boolean. In case the input is a nanny object, it indicates Whether the element column will be element or feature column}

\item{transform}{A function to use to tranforma the data internalli (e.g., log1p)}

\item{scale}{A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.}

\item{action}{A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).}

\item{...}{Further parameters passed to the function prcomp if you choose method="PCA" or Rtsne if you choose method="tSNE"}
}
\value{
A tbl object with additional columns for the reduced dimensions

A tbl object with additional columns for the reduced dimensions

A tbl object with additional columns for the reduced dimensions
}
\description{
reduce_dimensions() takes as input a `tbl` formatted as | <element> | <feature> | <value> | <...> | and calculates the reduced dimensional space of the feature value. The functions available are PCA, MDS (Robinson et al., 2010, <doi:10.1093/bioinformatics/btp616>), tSNE (Laurens van der Maaten, 2009)
}
\details{
\lifecycle{maturing}

This function reduces the dimensions of the feature values.
It can use multi-dimensional scaling (MDS) of principal component analysis (PCA).
}
\examples{


  reduce_dimensions(mtcars_tidy, car_model, feature, value, method="PCA")
  
  reduce_dimensions(mtcars_tidy, car_model, feature, value, method="MDS")
  
  reduce_dimensions(mtcars_tidy, car_model, feature, value, method="tSNE")



}
