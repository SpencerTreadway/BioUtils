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

#' Interpret Statistical Results
#'
#' Combines p-value and effect size into a human-readable interpretation.
#'
#' @param p.value Numeric.
#' @param effect.size Numeric.
#' @param alpha Significance threshold.
#'
#' @return Character interpretation.
#' @export
interpret.results <- function(p.value, effect.size, alpha = 0.05) {

  significance <- ifelse(p.value < alpha, "statistically significant", "not statistically significant")

  magnitude <- classify.effect.size(effect.size)

  interpretation <- paste(
    "The difference is", significance,
    "with a", magnitude, "effect size."
  )

  return(interpretation)
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
#' Performs statistical testing, computes effect size,
#' and generates an interpretation of results.
#'
#' @param df Data frame with Expression and Group columns.
#' @param alpha Significance level.
#'
#' @return List containing full analysis results.
#' @export
analyze.gene <- function(df, alpha=0.05) {

  test.results <- adaptive.t.test(df, alpha)

  d <- compute.effect.size(df)

  interpretation <- interpret.results(
    p.value = test.results$p.value,
    effect.size = d,
    alpha = alpha
  )

  return(list(
    p.value = test.results$p.value,
    effect.size = d,
    effect.size.class = classify.effect.size(d),
    interpretation = interpretation,
    raw = test.results
  ))
}
