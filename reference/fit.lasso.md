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
geo <- extract.expression(load.geo.soft(accession = "GDS507", log.transform = TRUE))
#> GDS507 not found locally, downloading from NCBI GEO...
#> Using locally cached version of GDS507 found here:
#> /tmp/RtmpYU9fuS/GDS507.soft.gz 
#> Using locally cached version of GPL97 found here:
#> /tmp/RtmpYU9fuS/GPL97.annot.gz 
phenotype.binary <- ifelse(geo$phenotype$disease.state == "disease", 1, 0)
lasso.fit <- fit.lasso(geo$expression, phenotype.binary)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'drop': non-conformable arguments
selected <- coef(lasso.fit, s = "lambda.1se")
#> Error: object 'lasso.fit' not found
selected[selected[, 1] != 0, ]
#> Error: object 'selected' not found
# }
```
