% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/id_mapping.R
\name{ambiguity}
\alias{ambiguity}
\title{Inspect the ambiguity of a mapping}
\usage{
ambiguity(
  d,
  from_col,
  to_col,
  groups = NULL,
  quantify = TRUE,
  qualify = TRUE,
  expand = NULL,
  global = FALSE,
  summary = FALSE
)
}
\arguments{
\item{d}{Data frame: a data frame with two columns to be inspected. It might
contain arbitrary other columns. Existing grouping will be removed.}

\item{from_col}{Character: column name of the "from" side of the mapping.}

\item{to_col}{Character: column name of the "to" side of the mapping.}

\item{groups}{Character vector of column names. Inspect ambiguity within
these groups; by default, ambiguity is determined across all rows.}

\item{quantify}{Logical or character: inspect the mappings for each
ID for ambiguity. If TRUE, for each translated column, two new columns
will be created with numeric values, representing the ambiguity of the
mapping on the "from" and "to" side of the translation, respectively.
If a character value provided, it will be used as a column name suffix
for the new columns.}

\item{qualify}{Logical or character: inspect the mappings for each
ID for ambiguity. If TRUE, for each translated column, a new column
will be inculded with values `one-to-one`, `one-to-many`, `many-to-one`
or `many-to-many`. If a character value provided, it will be used as a
column name suffix for the new column.}

\item{expand}{Logical: override the expansion of target columns, including
`to_col`: by default, this function expands data into multiple rows if
the `to_col` has already been expanded. Using this argument, the
`to_col` and other target columns will be lists of vectors for `expand =
FALSE`, and simple vectors for `expand = TRUE`.}

\item{global}{Logical or character: if `groups` are provided, analyse
ambiguity also globally, across the whole data frame. Character value
provides a custom suffix for the columns quantifying and qualifying
global ambiguity.}

\item{summary}{Logical: generate a summary about the ambiguity of the
translation and make it available as an attribute.}
}
\value{
A data frame (tibble) with ambiguity information added in new
    columns, as described at the "quantify" and "qualify" arguments.
}
\description{
Inspect the ambiguity of a mapping
}
