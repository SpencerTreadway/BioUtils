#' Build Analysis Data Frame
#'
#' Combines gene expression values with phenotype group labels
#' into a clean analysis-ready data frame.
#'
#' @param expression.vector Numeric vector of gene expression values.
#' @param phenotype Data frame of sample metadata.
#'
#' @return Data frame with Expression and Group columns.
#'
#' @export
build.analysis.df <- function(expression.vector, phenotype)
{
  if(class(expression.vector) != "numeric")
  {
    print("Expression vector must be numeric")
    return()
  }
  analysis.df <- data.frame(expression=expression.vector, group=phenotype$disease.state)
  analysis.df <- analysis.df[!is.nan(analysis.df$expression),]
  analysis.df <- na.omit(analysis.df)
  return(analysis.df)
}
