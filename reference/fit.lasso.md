# Fit LASSO Regression for Multi-Gene Biomarker Discovery

Applies cross-validated LASSO logistic regression to identify a sparse
set of genes that jointly predict a binary phenotype. Unlike univariate
tests, LASSO accounts for the joint contribution of all genes
simultaneously, automatically penalizing redundant predictors.

## Usage

``` r
fit.lasso(expression.matrix, phenotype.vector, alpha = 1, nfolds = 10)
```

## Arguments

- expression.matrix:

  Numeric matrix of gene expression values as returned by
  `extract.expression()$expression`. Rows are probes, columns are
  samples. The matrix is transposed internally so samples become rows as
  required by `glmnet`.

- phenotype.vector:

  Factor or integer vector of binary phenotype labels, one per sample
  (e.g., `0` for control, `1` for disease). Must be the same length as
  the number of columns in `expression.matrix`.

- alpha:

  Numeric. Elastic net mixing parameter. `1` (default) gives pure LASSO;
  `0` gives Ridge; values between 0 and 1 give elastic net.

- nfolds:

  Integer. Number of cross-validation folds. Default is `10`.

## Value

A `cv.glmnet` object. Key elements:

- lambda.min:

  Lambda with minimum cross-validated error.

- lambda.1se:

  Largest lambda within 1 SE of the minimum (sparser model).

- glmnet.fit:

  The full sequence of fitted models across lambda values.

Extract selected genes with `coef(fit, s = "lambda.1se")`; probes with
non-zero coefficients are the LASSO-selected biomarkers.

## Details

LASSO is suited to high-dimensional genomic data where the number of
genes far exceeds samples. It shrinks most coefficients to exactly zero,
retaining only genes with independent predictive value. This complements
[`run.limma.de()`](https://github.com/SpencerTreadway/BioUtils/reference/run.limma.de.md) -
while limma ranks genes by individual differential expression, LASSO
identifies the minimal subset with non-redundant joint predictive
signal.

## Examples

``` r
# \donttest{
# Synthetic example — small expression matrix with binary outcome
set.seed(42)
expr.mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(expr.mat) <- paste0("probe", 1:20)
colnames(expr.mat) <- paste0("sample", 1:10)
phenotype.binary <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
lasso.fit <- fit.lasso(expr.mat, phenotype.binary)
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
coef(lasso.fit, s = "lambda.1se")
#> 21 x 1 sparse Matrix of class "dgCMatrix"
#>             lambda.1se
#> (Intercept)          .
#> probe1               .
#> probe2               .
#> probe3               .
#> probe4               .
#> probe5               .
#> probe6               .
#> probe7               .
#> probe8               .
#> probe9               .
#> probe10              .
#> probe11              .
#> probe12              .
#> probe13              .
#> probe14              .
#> probe15              .
#> probe16              .
#> probe17              .
#> probe18              .
#> probe19              .
#> probe20              .
# }
```
