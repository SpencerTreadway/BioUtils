# Plot Volcano Plot from Differential Expression Results

Generates a volcano plot visualizing the relationship between
statistical significance and fold-change magnitude for all probes in a
differential expression analysis. Points are colored by whether they
exceed user-defined significance and fold-change thresholds.

## Usage

``` r
volcano.plot(de.results, fc.threshold = 1, fdr.threshold = 0.05)
```

## Arguments

- de.results:

  Data frame as returned by
  [`run.limma.de()`](https://github.com/SpencerTreadway/BioUtils/reference/run.limma.de.md),
  containing at minimum the columns `logFC` and `adj.P.Val`.

- fc.threshold:

  Numeric. Absolute log2 fold-change cutoff for labeling genes as
  differentially expressed. Default is `1` (i.e., 2-fold).

- fdr.threshold:

  Numeric. Adjusted p-value cutoff. Default is `0.05`.

## Value

A ggplot object showing each probe as a point with `logFC` on the x-axis
and `-log10(adj.P.Val)` on the y-axis. Points are colored as
upregulated, downregulated, or not significant. Dashed threshold lines
are drawn at the specified cutoffs.

## Details

The volcano plot is the standard summary visualization for a TopTable
from
[`run.limma.de()`](https://github.com/SpencerTreadway/BioUtils/reference/run.limma.de.md)
and serves as the genome-wide complement to
[`gene.analysis.plot()`](https://github.com/SpencerTreadway/BioUtils/reference/gene.analysis.plot.md),
which visualizes a single gene. Probes that appear in the upper corners
(high fold-change, low FDR) are the strongest candidates for follow-up
with
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md)
or inclusion in a gene signature for
[`run.gsea()`](https://github.com/SpencerTreadway/BioUtils/reference/run.gsea.md).

## Examples

``` r
if (FALSE) { # \dontrun{
de.results <- run.limma.de(eset, condition.col = "disease.state")
volcano.plot(de.results, fc.threshold = 1, fdr.threshold = 0.05)
} # }
```
