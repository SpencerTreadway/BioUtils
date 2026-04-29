# Adjust P-Values for Multiple Comparisons

Applies a multiple testing correction to a vector of raw p-values,
reducing the false discovery rate when many genes are tested
simultaneously.

## Usage

``` r
adjust.pvalues(p.values, method = "BH")
```

## Arguments

- p.values:

  Numeric vector of raw p-values, one per gene or probe.

- method:

  Character. Correction method passed to
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Common options
  are `"BH"` (Benjamini-Hochberg), `"bonferroni"`, and `"holm"`. Default
  is `"BH"`.

## Value

A numeric vector of adjusted p-values, the same length and order as the
input. For the BH method values represent the estimated false discovery
rate (FDR); for Bonferroni and Holm they represent the family-wise error
rate.

## Details

When
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md)
is applied across many probes in a loop, each test is performed
independently and p-values are not corrected for multiplicity. This
function should be applied to the resulting p-value vector before
interpreting significance. BH correction is recommended for exploratory
genomic analyses; Bonferroni is more conservative and suited to
confirmatory settings. For genome-wide analysis, prefer
[`run.limma.de()`](https://github.com/SpencerTreadway/BioUtils/reference/run.limma.de.md)
which handles correction internally.

## Examples

``` r
if (FALSE) { # \dontrun{
geo <- extract.expression(load.geo.soft("GDS3268.soft"))
probe.ids <- find.probe.by.gene(geo$gene, c("GENE1", "GENE2", "GENE3"))

raw.pvals <- sapply(probe.ids, function(id) {
  expr <- get.gene.expression(geo$expression, id)
  df <- build.analysis.df(expr, geo$phenotype, geo$gene)
  analyze.gene(df)$p.value
})
adj.pvals <- adjust.pvalues(raw.pvals, method = "BH")
} # }
```
