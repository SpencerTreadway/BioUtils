#' Compute Effect Size (Cohen's d)
#'
#' Calculates the standardized mean difference between two groups.
#'
#' @param df Data frame with Expression and Group columns.
#'
#' @return Numeric. Cohen's d value.
#' @export
compute.effect.size <- function(df)
{
  groups <- split(df$expression, df$group)

  g1 <- groups[[1]]
  g2 <- groups[[2]]

  mean1 <- mean(g1)
  mean2 <- mean(g2)

  sd1 <- sd(g1)
  sd2 <- sd(g2)

  n1 <- length(g1)
  n2 <- length(g2)

  pooled.sd <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2))
  d <- (mean1 - mean2) / pooled.sd

  return(d)
}

#' Compute Confidence Interval for Effect Size
#'
#' Estimates a confidence interval for Cohen's d using bootstrap resampling.
#' This provides a measure of uncertainty around the effect size estimate,
#' which is especially useful when sample sizes are small or distributions
#' deviate from normality.
#'
#' @param df A data frame containing at least two columns:
#' \describe{
#'   \item{expression}{Numeric vector of gene expression values}
#'   \item{group}{Factor or character vector indicating group membership}
#' }
#' @param n.boot Integer. Number of bootstrap resamples to perform.
#' Default is 1000.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{lower}{Lower bound of the confidence interval}
#'   \item{upper}{Upper bound of the confidence interval}
#' }
#'
#' @details
#' The function repeatedly resamples the input data with replacement and
#' recomputes Cohen's d for each resample. The confidence interval is
#' derived from the empirical distribution of bootstrapped effect sizes.
#'
#' @examples
#' \dontrun{
#' ci <- compute.ci(df, n.boot = 1000)
#' }
#'
#' @export
compute.ci <- function(df, alpha=0.05, n.boot=1000)
{
  boot.ds <- numeric(n.boot)

  for(i in 1:n.boot)
  {
    sample.df <- df[sample(nrow(df), replace=TRUE),]
    boot.ds[i] <- compute.effect.size(sample.df)
  }

  ci <- quantile(boot.ds, c(alpha / 2, 1- (alpha / 2)))

  return(list(
    lower = ci[1],
    upper = ci[2]
  ))
}

#' Perform Nonparametric Group Comparison
#'
#' Conducts a Wilcoxon rank-sum test (Mann-Whitney U test) to compare
#' expression values between two groups without assuming normality.
#'
#' @param df A data frame containing:
#' \describe{
#'   \item{expression}{Numeric vector of gene expression values}
#'   \item{group}{Factor or character vector with exactly two groups}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{p.value}{P-value from the Wilcoxon test}
#'   \item{test}{Character string identifying the test used}
#' }
#'
#' @details
#' This test is useful as a robustness check against parametric methods
#' such as the t-test, especially when data are skewed or contain outliers.
#'
#' @examples
#' \dontrun{
#' res <- nonparametric.test(df)
#' }
#'
#' @export
nonparametric.test <- function(df)
{
  res <- wilcox.test(expression ~ group, data = df)

  return(list(
    p.value = res$p.value,
    test = "Wilcoxon rank-sum"
  ))
}

#' Classify Effect Size
#'
#' Categorizes Cohen's d into qualitative magnitude levels.
#'
#' @param d Numeric. Effect size.
#'
#' @return Character description.
#' @export
classify.effect.size <- function(d) {

  abs.d <- abs(d)

  if (abs.d < 0.2)
  {
    return("negligible")
  }
  else if (abs.d < 0.5)
  {
    return("small")
  }
  else if (abs.d < 0.8)
  {
    return("moderate")
  }
  else
  {
    return("large")
  }
}

#' Assess Biological Relevance of Gene Expression Differences
#'
#' Provides a heuristic interpretation of whether a statistically
#' significant difference in gene expression is likely to be
#' biologically meaningful based on effect size and p-value.
#'
#' @param effect.size Numeric. Cohen's d or other standardized effect size.
#' @param p.value Numeric. P-value from a statistical test.
#' @param alpha Numeric. Significance threshold. Default is 0.05.
#'
#' @return Character string describing the inferred biological relevance.
#'
#' @details
#' This function applies simple thresholds to distinguish between:
#' \itemize{
#'   \item Statistically significant but small effects
#'   \item Potentially meaningful biological differences
#'   \item Lack of strong evidence
#' }
#'
#' These rules are heuristic and intended as guidance rather than
#' definitive biological conclusions.
#'
#' @examples
#' \dontrun{
#' flag <- flag.biological.relevance(effect.size = 0.6, p.value = 0.01)
#' }
#'
#' @export
flag.biological.relevance <- function(effect.size, p.value)
{

  if (p.value < 0.05 && abs(effect.size) > 0.5)
  {
    return("Potentially biologically meaningful")
  }

  if (p.value < 0.05 && abs(effect.size) <= 0.5)
  {
    return("Statistically significant but small effect")
  }

  return("No strong evidence of biological relevance")
}

#' Perform Adaptive T-Test
#'
#' Automatically determines whether to use a Welch or Student
#' t-test based on variance equality.
#'
#' @param df Analysis data frame with Expression and Group columns.
#' @param alpha Significance level for variance test.
#'
#' @return List containing:
#' \describe{
#'   \item{variance.test}{Result of variance test}
#'   \item{t.test}{Result of t-test}
#'   \item{test.type}{Character indicating test used}
#' }
#'
#' @export
adaptive.t.test <- function(df, alpha=0.05)
{
  variance.result <- var.test(expression ~ group, data=df)
  p.value <- variance.result$p.value

  if(p.value < 0.05)
  {
    t.result <- t.test(expression ~ group, data=df, var.equal=FALSE, conf.level=1-alpha)
    return(t.result)
  }
  else
  {
    t.result <- t.test(expression ~ group, data=df, var.equal=TRUE, conf.level=1-alpha)
    anova.result <- oneway.test(expression ~ group, data=df, var.equal=TRUE)
    return(list(
      variance.test=anova.result,
      t.test=t.result,
      p.value=t.result$p.value
    ))
  }
}

#' Analyze Gene Expression Between Groups
#'
#' Performs a comprehensive statistical analysis comparing gene expression
#' between two groups. This includes parametric and nonparametric testing,
#' effect size estimation, confidence intervals, and interpretation.
#'
#' @param df A data frame containing:
#' \describe{
#'   \item{Expression}{Numeric vector of gene expression values}
#'   \item{Group}{Factor or character vector with exactly two groups}
#' }
#' @param alpha Numeric. Significance threshold. Default is 0.05.
#' @param n_boot Integer. Number of bootstrap samples for confidence interval.
#' Default is 1000.
#'
#' @return A named list containing:
#' \describe{
#'   \item{p.value}{P-value from t-test}
#'   \item{effect.size}{Cohen's d}
#'   \item{effect.size_class}{Qualitative magnitude of effect size}
#'   \item{confidence.interval}{List with lower and upper bounds}
#'   \item{test.type}{Type of t-test used}
#'   \item{nonparametric.p}{P-value from Wilcoxon test}
#'   \item{robustness}{Agreement between parametric and nonparametric tests}
#'   \item{biological.relevance}{Heuristic interpretation}
#'   \item{interpretation}{Human-readable summary}
#'   \item{raw}{Raw statistical test outputs}
#' }
#'
#' @details
#' This function integrates multiple statistical perspectives to provide
#' a more complete understanding of group differences. It goes beyond
#' simple hypothesis testing by including effect size estimation,
#' uncertainty quantification, and robustness checks.
#'
#' @examples
#' \dontrun{
#' result <- analyze.gene(df)
#' result$interpretation
#' }
#'
#' @export
analyze.gene <- function(df, alpha=0.05, n.boot=1000)
{
  if (!all(c("expression", "group") %in% colnames(df)))
  {
    stop("Data frame must contain 'Expression' and 'Group' columns.")
  }

  if (length(unique(df$group)) != 2)
  {
    stop("Group column must contain exactly two groups.")
  }

  test.results <- adaptive.t.test(df, alpha)
  p.value <- test.results$p.value

  d <- compute.effect.size(df)
  d.class <- classify.effect.size(d)

  ci <- compute.ci(df, alpha, n.boot)

  np.result <- nonparametric.test(df)
  np.p <- np.result$p.value

  if (p.value < alpha && np.p < alpha)
  {
    robustness <- "robust across parametric and nonparametric tests"
  }
  else if (p.value < alpha && np.p >= alpha)
  {
    robustness <- "sensitive to distributional assumptions"
  }
  else
  {
    robustness <- "no consistent evidence of difference"
  }

  bio.flag <- flag.biological.relevance(d, p.value)

  interpretation <- paste(
    "The difference is ",
    ifelse(p.value < alpha, "statistically significant ", "not statistically significant "),
    "with a ", d.class, " effect size.",
    "\nResult is ", robustness, ".\n",
    bio.flag, ".", sep=""
  )

  return(list(
    p.value = p.value,
    effect.size = d,
    effect.size.class = d.class,
    confidence.interval = ci,
    nonparametric.p = np.p,
    robustness = robustness,
    biological.relevance = bio.flag,
    interpretation = interpretation,
    raw = list(
      parametric = test.results,
      nonparametric = np.result
    )
  ))
}
