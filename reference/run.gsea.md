# Run Gene Set Enrichment Analysis

Tests whether predefined gene sets (e.g., pathways or GO terms) are
systematically enriched among genes ranked by differential expression,
using the fgseaMultilevel algorithm from `fgsea`.

## Usage

``` r
run.gsea(de.results, genes, pathways, min.size = 15, max.size = 500)
```

## Arguments

- de.results:

  Data frame as returned by
  [`run.limma.de()`](https://github.com/SpencerTreadway/BioUtils/reference/run.limma.de.md).
  Must contain a `logFC` column. Row names are probe IDs.

- genes:

  Data frame of gene annotations as returned by
  `extract.expression()$gene`. Must contain columns `"ID"` for probe
  identifiers and `"Gene symbol"` for gene symbols matching the format
  used by pathway databases such as MSigDB.

- pathways:

  Named list of character vectors, where each element is a gene set and
  names are pathway labels. MSigDB or GO gene sets in this format can be
  loaded via `msigdbr`.

- min.size:

  Integer. Minimum number of genes required in a gene set for it to be
  tested. Default is `15`.

- max.size:

  Integer. Maximum gene set size. Default is `500`.

## Value

A data frame with one row per tested pathway, containing:

- pathway:

  Gene set name.

- pval:

  Nominal p-value.

- padj:

  BH-adjusted p-value across all tested pathways.

- NES:

  Normalized enrichment score. Positive values indicate enrichment among
  upregulated genes; negative values indicate enrichment among
  downregulated genes.

- size:

  Number of genes from the pathway present in the ranked list.

- leadingEdge:

  List-column of genes driving the enrichment score.

## Details

Uses the `fgseaMultilevel` algorithm, which estimates p-values using an
adaptive multilevel splitting Monte Carlo approach. This is more
accurate than the original permutation-based `fgseaSimple` and does not
require a fixed permutation count.

GSEA operates on the full ranked gene list from
[`run.limma.de()`](https://github.com/SpencerTreadway/BioUtils/reference/run.limma.de.md)
rather than a filtered subset, preserving information from all genes.
Probe IDs from the TopTable row names are resolved to gene symbols via
the `genes` annotation data frame before matching against pathway gene
sets. Probes with no gene annotation are silently dropped from the
ranked list, as they cannot be matched to any pathway. The `leadingEdge`
genes are strong candidates for follow-up with
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md)
or
[`gene.analysis.plot()`](https://github.com/SpencerTreadway/BioUtils/reference/gene.analysis.plot.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(msigdbr)
geo <- extract.expression(load.geo.soft(accession = "GDS507",
                                                 log.transform = TRUE))
hallmark.df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark.df$gene_symbol, hallmark.df$gs_name)
de.results <- run.limma.de(geo)
gsea.results <- run.gsea(de.results, geo$gene, pathways)
head(gsea.results[order(gsea.results$padj), ])
} # }
```
