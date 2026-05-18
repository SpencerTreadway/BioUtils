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
# \donttest{
library(msigdbr)
geo <- extract.expression(load.geo.soft(accession = "GDS507",
                                                 log.transform = TRUE))
#> GDS507 not found locally, downloading from NCBI GEO...
#> Using locally cached version of GDS507 found here:
#> /tmp/RtmpxRZSjV/GDS507.soft.gz 
#> Using locally cached version of GPL97 found here:
#> /tmp/RtmpxRZSjV/GPL97.annot.gz 
hallmark.df <- msigdbr(species = "Homo sapiens", category = "H")
#> Warning: The `category` argument of `msigdbr()` is deprecated as of msigdbr 10.0.0.
#> ℹ Please use the `collection` argument instead.
pathways <- split(hallmark.df$gene_symbol, hallmark.df$gs_name)
de.results <- run.limma.de(geo)
#> Removing intercept from test coefficients
gsea.results <- run.gsea(de.results, geo$gene, pathways)
head(gsea.results[order(gsea.results$padj), ])
#>                                       pathway         pval         padj
#>                                        <char>        <num>        <num>
#> 1:         HALLMARK_INTERFERON_GAMMA_RESPONSE 8.505567e-11 3.827505e-09
#> 2:               HALLMARK_ALLOGRAFT_REJECTION 2.282585e-07 5.135817e-06
#> 3:           HALLMARK_TNFA_SIGNALING_VIA_NFKB 2.844309e-06 4.266464e-05
#> 4:         HALLMARK_INTERFERON_ALPHA_RESPONSE 1.862513e-04 2.095328e-03
#> 5: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION 9.521288e-04 8.569159e-03
#> 6:                           HALLMARK_HYPOXIA 1.516218e-03 8.952530e-03
#>      log2err         ES       NES  size
#>        <num>      <num>     <num> <int>
#> 1: 0.8390889 -0.6451016 -2.551362    67
#> 2: 0.6901325 -0.6559719 -2.299225    38
#> 3: 0.6272567 -0.5518351 -2.126387    59
#> 4: 0.5188481 -0.5941171 -1.983752    32
#> 5: 0.4772708 -0.4610456 -1.814634    64
#> 6: 0.4550599 -0.4486926 -1.759598    63
#>                                     leadingEdge
#>                                          <list>
#> 1: ST8SIA4,IRF1,SAMHD1,SSPN,GBP4,SLAMF7,...[42]
#> 2:    ST8SIA4,TLR3,CTSS,EGFR,GBP2,BCAT1,...[15]
#> 3: IRF1,DUSP4,INHBA,EDN1,MARCKS,SLC16A6,...[36]
#> 4:    GBP2,IRF1,GBP4,RSAD2,CMPK2,SAMD9L,...[18]
#> 5:  NNMT,CTHRC1,ACTA2,EDIL3,INHBA,PRRX1,...[21]
#> 6:     ANGPTL4,CP,PDK1,NR3C1,ALDOA,EGFR,...[22]
# }
```
