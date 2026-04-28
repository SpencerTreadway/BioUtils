#' Find Probe IDs by Gene Name
#'
#' Searches the gene annotation data frame to find the probe IDs corresponding
#' to one or more gene names or descriptions.
#'
#' @param genes Data frame of gene annotations as returned by
#'   \code{extract.expression()$gene}. The second column (index 2) is expected
#'   to contain gene titles or descriptions.
#' @param gene.names Character vector of one or more gene names or descriptions
#'   to search for. Matching is exact and case-sensitive.
#'
#' @return Integer vector of probe IDs corresponding to the matched genes.
#'   The length of the vector equals the number of matches found. Returns an
#'   empty integer vector if no matches are found.
#'
#' @details
#' Probe IDs are the integer row names of the \code{genes} annotation data
#' frame. These IDs are used directly to index rows of the expression matrix
#' in \code{get.gene.expression()}. If multiple gene names are supplied, all
#' matching probe IDs are returned as a vector, which is then passed as-is to
#' \code{get.gene.expression()} to retrieve a multi-row expression matrix.
#'
#' @examples
#' \dontrun{
#' geo <- extract.expression(load.geo.soft("GDS3268.soft"))
#'
#' # Single gene
#' probe <- find.probe.by.gene(geo$gene, "mucin 20, cell surface associated")
#'
#' # Multiple genes
#' probes <- find.probe.by.gene(geo$gene, c(
#'   "mucin 20, cell surface associated",
#'   "alcohol dehydrogenase 1A (class I), alpha polypeptide"
#' ))
#' }
#'
#' @export
find.probe.by.gene <- function(genes, gene.names)
{
  return(as.integer(rownames(genes[which(genes[, 2] %in% gene.names), ])))
}

#' Extract Gene Expression by Probe ID
#'
#' Retrieves expression values for one or more genes from the expression matrix
#' using their integer probe IDs.
#'
#' @param expression Numeric matrix of gene expression values as returned by
#'   \code{extract.expression()$expression}. Rows are probe IDs, columns are
#'   sample IDs.
#' @param probe.id Integer vector of one or more probe IDs as returned by
#'   \code{find.probe.by.gene()}.
#'
#' @return A numeric matrix with one row per probe and one column per sample.
#'   Row names are the probe IDs; column names are the sample IDs. When a
#'   single probe ID is supplied the result is still a matrix (not a vector),
#'   ensuring that \code{build.analysis.df()} receives a consistent input type.
#'
#' @details
#' The returned matrix is the direct input for \code{build.analysis.df()}.
#' Using \code{as.matrix()} ensures the single-probe case returns a one-row
#' matrix rather than a named vector, which keeps the \code{pivot_longer}
#' step in \code{build.analysis.df()} consistent regardless of how many genes
#' are queried.
#'
#' @examples
#' \dontrun{
#' geo <- extract.expression(load.geo.soft("GDS3268.soft"))
#' probes <- find.probe.by.gene(geo$gene, c(
#'   "mucin 20, cell surface associated",
#'   "alcohol dehydrogenase 1A (class I), alpha polypeptide"
#' ))
#' expr <- get.gene.expression(geo$expression, probes)
#' dim(expr)  # probes x samples
#' }
#'
#' @export
get.gene.expression <- function(expression, probe.id)
{
  return(as.matrix(expression[probe.id, ]))
}
