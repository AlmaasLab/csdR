% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{partial_argsort}
\alias{partial_argsort}
\title{Extract indecies corresponding to the largest elements}
\usage{
partial_argsort(x, n_elements)
}
\arguments{
\item{x}{Numeric vector, the vector containing the numbers to sort.}

\item{n_elements}{Integer scalar, the number of indecies to return.}
}
\value{
Numeric vector, the indecies of the largest elements (in sorted order) in
\code{x}.
}
\description{
Extracts the indecies of the \eqn{n} largest elements of the input
This procedure is equivalent to \code{order(x, decreasing = TRUE)[1:n_elements]},
but is much faster and avoids the overhead of sorting discarded elements.
This function is useful for extracting the rows in a data frame having the
largest values in one of the columns.
}
\examples{
x <- c(10L,5L,-2L,12L,15L)
max_indecies <- partial_argsort(x,3L)
max_indecies
x[max_indecies]
order(x)[1:3]
mtcars[partial_argsort(mtcars$hp,5L),]
}
