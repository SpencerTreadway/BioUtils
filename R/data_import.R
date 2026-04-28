#' Load GEO Dataset from SOFT File
#'
#' Imports a GEO dataset from a .soft file and converts it into
#' a structured ExpressionSet object for downstream analysis.
#'
#' This function abstracts the workflow of using GEOquery to read
#' raw GEO data and prepares it for consistent use within BioUtils.
#'
#' @param file.path Character. Path to a GEO .soft file.
#' @param log.transform Logical. Whether to apply log2 transformation to the
#'   expression matrix during import. Default is \code{FALSE}. Set to
#'   \code{TRUE} if the raw data is not already on a log scale, which is
#'   common for Affymetrix and Agilent arrays.
#'
#' @return An \code{ExpressionSet} object containing expression data,
#'   phenotype data, and feature data. This object is the required input
#'   for \code{extract.expression()}.
#'
#' @details
#' Internally delegates to \code{GEOquery::getGEO()} to parse the .soft file
#' and \code{GEOquery::GDS2eSet()} to construct the \code{ExpressionSet}.
#' The returned object is passed directly to \code{extract.expression()} to
#' decompose it into the list format used throughout BioUtils.
#'
#' @examples
#' \dontrun{
#' eset <- load.geo.soft("GDS3268.soft")
#' eset <- load.geo.soft("GDS3268.soft", log.transform = TRUE)
#' }
#'
#' @export
load.geo.soft <- function(file.path, log.transform=FALSE)
{
  if(!file.exists(file.path))
  {
    print(paste(file.path, "not found", sep=" "))
    return()
  }

  return(GEOquery::GDS2eSet(GEOquery::getGEO(filename=file.path), do.log2=log.transform))
}

#' Extract Expression and Metadata from GEO Dataset
#'
#' Decomposes an \code{ExpressionSet} object into its core components,
#' returning a named list that is used as the primary data structure
#' throughout BioUtils.
#'
#' @param eset An \code{ExpressionSet} object as returned by
#'   \code{load.geo.soft()}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{expression}{Numeric matrix of gene expression values. Rows are
#'     probe IDs, columns are sample IDs.}
#'   \item{phenotype}{Data frame of sample metadata (e.g., disease state,
#'     age, tissue). Row names are sample IDs matching the column names of
#'     \code{expression}.}
#'   \item{gene}{Data frame of probe-level gene annotations (e.g., gene
#'     symbol, gene title, chromosomal location). Row names are probe IDs
#'     matching the row names of \code{expression}.}
#' }
#'
#' @details
#' The three returned components form the basis of all downstream BioUtils
#' workflows. The \code{gene} data frame is passed to
#' \code{find.probe.by.gene()} for probe lookup, \code{expression} is
#' passed to \code{get.gene.expression()}, and \code{phenotype} is passed
#' to \code{build.analysis.df()} for group labeling.
#'
#' @examples
#' \dontrun{
#' eset <- load.geo.soft("GDS3268.soft")
#' geo <- extract.expression(eset)
#' head(geo$phenotype)
#' dim(geo$expression)
#' }
#'
#' @export
extract.expression <- function(eset)
{
  if(!inherits(eset, "ExpressionSet"))
  {
    print(paste("Invalid Argument of Class:", class(eset), sep=" "))
    return()
  }

  gene.expression <- Biobase::exprs(eset)
  phenotype <- Biobase::pData(eset)
  genes <- Biobase::fData(eset)

  return(list(
    expression = gene.expression,
    phenotype = phenotype,
    gene = genes
  ))
}
