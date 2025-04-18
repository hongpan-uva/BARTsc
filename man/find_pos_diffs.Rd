% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bartsc_find_differential.R
\name{find_pos_diffs}
\alias{find_pos_diffs}
\title{find differentially expressed features}
\usage{
find_pos_diffs(
  object,
  mod,
  fore_celltype,
  back_celltype,
  fore_cells,
  back_cells,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  random.seed = 1
)
}
\arguments{
\item{object}{A bartsc object.}
}
\value{
A character vector.
}
\description{
find differentially expressed features
}
