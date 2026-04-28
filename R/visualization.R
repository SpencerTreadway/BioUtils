#' Plot Gene Expression with Automatic Statistical Analysis
#'
#' Generates an annotated visualization of gene expression between two groups.
#' For a single-gene data frame, the plot is annotated with statistical results
#' from \code{analyze.gene()}. For a multi-gene data frame, the plot is
#' faceted by gene and annotations are omitted — use \code{analyze.gene()} on
#' each gene subset for per-gene statistics.
#'
#' @param df A data frame as produced by \code{build.analysis.df()}, containing:
#' \describe{
#'   \item{expression}{Numeric vector of gene expression values.}
#'   \item{group}{Character or factor vector with exactly two group labels.}
#'   \item{gene}{Character vector of gene names. Optional; when present and
#'     containing more than one unique value, enables faceted multi-gene
#'     plotting.}
#' }
#'
#' @param alpha Numeric. Significance level passed to \code{analyze.gene()}.
#'   Default is \code{0.05}. Only used in single-gene mode.
#' @param n.boot Integer. Number of bootstrap resamples for the confidence
#'   interval, passed to \code{analyze.gene()}. Default is \code{1000}.
#'   Only used in single-gene mode.
#' @param show.points Logical. Whether to overlay jittered individual sample
#'   points on the violin/boxplot. Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object. In single-gene mode the plot includes a
#'   text annotation with p-value, Cohen's d, and bootstrapped CI, and a
#'   subtitle with the full interpretation string. In multi-gene mode the
#'   plot is faceted by gene with no statistical annotation.
#'
#' @details
#' Single-gene mode is triggered when the \code{gene} column is absent or
#' contains exactly one unique value. Multi-gene mode is triggered when
#' \code{gene} contains more than one unique value, producing a
#' \code{facet_wrap(~ gene)} layout. The plot theme is consistent across
#' both modes and matches the style conventions of the BioUtils package.
#'
#' @examples
#' \dontrun{
#' geo <- extract.expression(load.geo.soft("GDS3268.soft"))
#'
#' # Single-gene plot with statistical annotation
#' probe <- find.probe.by.gene(geo$gene, "mucin 20, cell surface associated")
#' expr <- get.gene.expression(geo$expression, probe)
#' df <- build.analysis.df(expr, geo$phenotype, geo$gene)
#' plot.gene.analysis(df)
#'
#' # Multi-gene faceted plot
#' probes <- find.probe.by.gene(geo$gene, c(
#'   "mucin 20, cell surface associated",
#'   "alcohol dehydrogenase 1A (class I), alpha polypeptide"
#' ))
#' expr.multi <- get.gene.expression(geo$expression, probes)
#' df.multi <- build.analysis.df(expr.multi, geo$phenotype, geo$gene)
#' plot.gene.analysis(df.multi)
#' }
#'
#' @export
plot.gene.analysis <- function(df, alpha=0.05, n.boot=1000, show.points=TRUE)
{
  multi.gene <- "gene" %in% colnames(df) && length(unique(df$gene)) > 1

  # Base plot shared by both modes
  p <- ggplot2::ggplot(df, ggplot2::aes(x=group, y=expression)) +
    ggplot2::geom_violin(trim=FALSE, alpha=0.6) +
    ggplot2::geom_boxplot(width=0.1, outlier.shape=NA)

  if(show.points)
  {
    p <- p + ggplot2::geom_jitter(width=0.15, alpha=0.4)
  }

  if(multi.gene)
  {
    p <- p +
      ggplot2::facet_wrap(~ gene, scales="free_y") +
      ggplot2::labs(
        title = "Gene Expression Comparison",
        x = "Group",
        y = "Expression"
      )
  }
  else
  {
    # Single-gene mode: run full statistical analysis and annotate
    result <- analyze.gene(df, alpha=alpha, n.boot=n.boot)
    ci <- result$confidence.interval

    label <- paste(
      "p = ", signif(result$p.value, 3),
      "\nd = ", round(result$effect.size, 2),
      " (", result$effect.size.class, ")",
      "\n", 1 - alpha, "% CI: [", signif(ci$lower, 2), ", ", signif(ci$upper, 2), "]"
    )

    y.pos <- max(df$expression, na.rm=TRUE)

    p <- p +
      ggplot2::annotate(
        "text",
        x = 1.5,
        y = y.pos - 0.2,
        label = label,
        hjust = 0.5,
        vjust = -0.5,
        size = 4
      ) +
      ggplot2::labs(
        title = "Gene Expression Comparison",
        subtitle = result$interpretation,
        x = "Group",
        y = "Expression"
      )
  }

  p <- p + ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color="black"),
    panel.grid.major = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color="black", linewidth=1.2),
    title = ggplot2::element_text(size=20),
    axis.title = ggplot2::element_text(size=16),
    plot.title = ggplot2::element_text(hjust=0.5)
  )

  return(p)
}
