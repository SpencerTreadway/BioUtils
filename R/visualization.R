#' Plot Gene Expression with Automatic Statistical Analysis
#'
#' Generates an annotated visualization of gene expression between two groups.
#' This function internally performs statistical analysis using \code{analyze_gene()}
#' and overlays key results on the plot.
#'
#' @param df A data frame containing:
#' \describe{
#'   \item{Expression}{Numeric vector of gene expression values}
#'   \item{Group}{Factor or character vector with exactly two groups}
#' }
#' @param alpha Numeric. Significance level. Default is 0.05.
#' @param n.boot Integer. Number of bootstrap samples for confidence interval.
#' @param show.points Logical. Whether to overlay jittered points.
#'
#' @return A ggplot object with annotated statistical results.
#'
#' @examples
#' \dontrun{
#' plot.gene.analysis(df)
#' }
#'
#' @export
plot.gene.analysis <- function(df, alpha=0.05, n.boot=1000, show.points=TRUE)
{
  result <- analyze.gene(df, alpha=alpha, n.boot=n.boot)
  ci <- result$confidence.interval

  label <- paste(
    "p = ", signif(result$p.value, 3),
    "\nd = ", round(result$effect.size, 2),
    " (", result$effect.size.class, ")",
    "\n", 1 - alpha, "% CI: [", signif(ci$lower, 2), ", ", signif(ci$upper, 2), "]"
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x=group, y=expression)) +
    ggplot2::geom_violin(trim=FALSE, alpha=0.6) +
    ggplot2::geom_boxplot(width=0.1, outlier.shape=NA)
  if(show.points)
  {
    p <- p + ggplot2::geom_jitter(width=0.15, alpha=0.4)
  }
  y.pos <- max(df$expression, na.rm=TRUE)

  p <- p +
    ggplot2::annotate("text",
             x=1.5,
             y=y.pos - 0.2,
             label=label,
             hjust=0.5,
             vjust=-0.5,
             size=4
            )
  p <- p +
    ggplot2::labs(
      title="Gene Expression Comparison",
      subtitle=result$interpretation,
      x="Group",
      y="Expression"
    )
  p <- p + ggplot2::theme(
    panel.background=ggplot2::element_blank(), # Blank background
    panel.border=ggplot2::element_rect(color="black"), # Black border
    panel.grid.major=ggplot2::element_blank(), # No grid lines
    axis.line=ggplot2::element_line(color="black",linewidth=1.2), # Thick black axes
    title=ggplot2::element_text(size=20), # Increase title size
    axis.title=ggplot2::element_text(size=16), # Slightly increase axis title size
    plot.title=ggplot2::element_text(hjust=0.5) # Specifically justify plot title
  )
  return(p)
}
