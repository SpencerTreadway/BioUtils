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
#'   \item{variance_test}{Result of variance test}
#'   \item{t_test}{Result of t-test}
#'   \item{test_type}{Character indicating test used}
#' }
#'
#' @export
adaptive.t.test <- function(df, alpha=0.05)
{
  variance_result <- var.test(expression ~ group, data=df)
  p.value <- variance_result$p.value

  if(p.value < 0.05)
  {
    t.result <- t.test(expression ~ group, data=df, var.equal=FALSE, conf.level=1-alpha)
    t.result
    return(t.result)
  }
  else
  {
    t.result <- t.test(expression ~ group, data=df, var.equal=TRUE, conf.level=1-alpha)
    anova.result <- oneway.test(expression ~ group, data=df, var.equal=TRUE)
    t.result
    anova.result
    return(list(t=t.result, anova=anova.result))
  }
}
