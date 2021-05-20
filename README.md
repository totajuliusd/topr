
<!-- README.md is generated from README.Rmd. Please edit that file -->

# topR

See full documentation at <https://wuxi-nextcode.github.io/topR/>

### Installing from github using `devtools`

``` r
devtools::install_github("wuxi-nextcode/topR")
```

## Example

In this example, we’re going to demonstrate the basic usage of the topR
library.

### Load packages

First load the `gorr` package, the `tidyverse` package is recommended in
general, but not required for this example

``` r
library(topR)
#> 
#> Attaching package: 'topR'
#> The following object is masked from 'package:stats':
#> 
#>     qqplot
library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
#> ✓ tibble  3.1.2     ✓ dplyr   1.0.4
#> ✓ tidyr   1.1.2     ✓ stringr 1.4.0
#> ✓ readr   1.4.0     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
```

### Loading and exploring prebuilt datasets

Load the `gwas_CD` and checking how many variants are in the datasets,
since we do not want to plot with too many datapoints

``` r
data(gwas_CD)
head(gwas_CD)
#>   CHROM     POS          ID       REF ALT        SE           P       OR
#> 1  chr1 1006415 rs145588482 TGGCAGCTC   T 0.1540620 0.000468758 0.583384
#> 2  chr1 1006415 rs145588482 TGGCAGCTC   T 0.1540620 0.000468758 0.583384
#> 3  chr1 1007256  rs76233940         G   A 0.1540250 0.000401567 0.579783
#> 4  chr1 1007256  rs76233940         G   A 0.1540250 0.000401567 0.579783
#> 5  chr1 1007256  rs76233940         G   A 0.1540250 0.000401567 0.579783
#> 6  chr1 1341559 rs376494450         C   T 0.0732974 0.000151216 1.320130
#>          AF Gene_Symbol Max_Impact         max_consequence
#> 1 0.0129317  AL645608.3     LOWEST downstream_gene_variant
#> 2 0.0129317       ISG15        LOW          intron_variant
#> 3 0.0130091  AL645608.1     LOWEST downstream_gene_variant
#> 4 0.0130091  AL645608.3     LOWEST downstream_gene_variant
#> 5 0.0130091       ISG15        LOW          intron_variant
#> 6 0.0270627        DVL1        LOW          intron_variant
```

``` r
paste("Number of SNPs in rhw dataset: [", length(gwas_CD$POS),"]", sep = "")
#> [1] "Number of SNPs in rhw dataset: [26821]"
```

### Manhattan plots

Get an overview of the crohn’s disease (CD) results in a Manhattan plot

``` r
manhattan(gwas_CD)
#> Scale for 'x' is already present. Adding another scale for 'x', which will
#> replace the existing scale.
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

get the top SNP per 10 MB

``` r
snps_CD <- get_best_snp_per_MB(gwas_CD, thresh = 1e-09, region = 10000000)
manhattan(gwas_CD, variants=snps_CD, annotation_thresh = 1e-09)
#> Scale for 'x' is already present. Adding another scale for 'x', which will
#> replace the existing scale.
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Chromplot

Take a closer look at the results by chromosome

``` r
CHR="chr16"
chromplot(gwas_CD,chr=CHR, variants = snps_CD, annotation_thresh = 5e-09)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

### Regionplot

Currently the regionplot uses the `gorr` library to fetch genes from
Genuity Science’ API to get the genes, so platform connection is needed,
see Article `Regionplot Example` for an example.

### Manhattan Multi

Display the output from more than one GWAS on the same manhattan plot

``` r
data(gwas_CD_IIBD)
manhattan(list(gwas_CD_IIBD,gwas_CD), legend_labels=  c("Crohns IIBD", "Crohns UKBB"), alpha = c(1, 0.5))
#> Scale for 'x' is already present. Adding another scale for 'x', which will
#> replace the existing scale.
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />
