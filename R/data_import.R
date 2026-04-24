#' Load GEO Dataset from SOFT File
#'
#' Imports a GEO dataset from a .soft file and converts it into
#' a structured ExpressionSet object for downstream analysis.
#'
#' This function abstracts the workflow of using GEOquery to read
#' raw GEO data and prepares it for consistent use within BioUtils.
#'
#' @param file.path Character. Path to a GEO .soft file.
#' @param log.transform Logical. Whether to apply log2 transformation.
#'
#' @return An ExpressionSet object containing expression data,
#' phenotype data, and feature data.
#'
#' @details
#' Internally uses GEOquery::getGEO() and GEOquery::GDS2eSet().
#'
#' @examples
#' eset <- load.geo.soft("GDS3268.soft")
#'
#' @export

library(GEOquery)

load.geo.soft <- function(file.path, log.transform=FALSE)
{
  if(!file.exists(file.path))
  {
    print(paste(filepath, "not found", sep=" "))
    return()
  }

  return(GDS2eSet(getGEO(filename=file.path), do.log2=log.transform))
}

#' Extract Expression and Metadata from GEO Dataset
#'
#' Decomposes an ExpressionSet object into its core components.
#'
#' @param eset ExpressionSet object.
#'
#' @return A list containing:
#' \describe{
#'   \item{expression}{Matrix of gene expression values}
#'   \item{phenotype}{Data frame of sample metadata}
#'   \item{features}{Data frame of gene/probe annotations}
#' }
#'
#' @export
extract.expression <- function(eset)
{
  if(class(eset) != "ExpressionSet")
  {
    print(paste("Invalid Argument of Class:", class(eset), sep=" "))
    return()
  }

  gene.expression <- exprs(eset)
  phenotype <- pData(eset)
  genes <- fData(eset)

  return(list(
    expression=gene.expression,
    phenotype=phenotype,
    gene=genes
  ))
}
