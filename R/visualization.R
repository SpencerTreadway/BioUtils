#' Plot Expression by Group
#'
#' Generates scatter and box plots of expression values grouped by phenotype.
#'
#' @param df Analysis data frame.
#' @param test.results test object to pull results from.
#'
#' @return Plot object(s).
#'
#' @export

library(ggplot2)

plot.expression <- function(df)
{
  p <- ggplot(df, aes(x=group, y=expression)) +
    geom_violin(trim=FALSE, alpha=0.6) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_jitter(width=0.15, alpha=0.4) +
    theme(
      panel.background=element_blank(), # Blank background
      panel.border=element_rect(color="black"), # Black border
      panel.grid.major=element_blank(), # No grid lines
      axis.line=element_line(color="black",linewidth=1.2), # Thick black axes
      title=element_text(size=20), # Increase title size
      axis.title=element_text(size=16), # Slightly increase axis title size
      plot.title=element_text(hjust=0.5) # Specifically justify plot title
    )
  return(p)
}
