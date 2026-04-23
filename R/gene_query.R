#' Find Probe ID by Gene Name
#'
#' Searches feature data to find the probe ID corresponding
#' to a given gene name or description.
#'
#' @param features Data frame of gene annotations.
#' @param gene_name Character. Target gene name.
#'
#' @return Character. Probe ID.
#'
#' @export
find.probe.by.gene <- function(features, gene.name)
{
  return(as.integer(rownames(features$gene[(features$gene[,2]==gene.name),])[1]))
}

#' Extract Gene Expression by Probe ID
#'
#' Retrieves expression values for a specific gene using its probe ID.
#'
#' @param expression Matrix of gene expression values.
#' @param probe_id Character. Probe identifier.
#'
#' @return Numeric vector of expression values.
#'
#' @export
get_gene_expression <- function(expression, probe.id)
{
  if(is.na(expression[probe.id]))
  {
    return()
  }

  return(expression[probe.id,])
}
