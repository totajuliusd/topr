
# topr

# This repo is deprecated and has been moved to https://github.com/GenuityScience/topr and released on CRAN.

See full documentation at <https://wuxi-nextcode.github.io/topr/>

### Installing from github using devtools

``` r
devtools::install_github("wuxi-nextcode/topr")
```

## Example

In this example we demonstrate the basic usage of the topr library.

### Load packages

First load the topr package, the tidyverse package is recommended in
general, but not required for this example

``` r
library(topr)
#> 
#> Attaching package: 'topr'
#> The following object is masked from 'package:stats':
#> 
#>     qqplot
library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
#> ✔ ggplot2 3.3.3     ✔ purrr   0.3.4
#> ✔ tibble  3.1.2     ✔ dplyr   1.0.6
#> ✔ tidyr   1.1.3     ✔ stringr 1.4.0
#> ✔ readr   1.4.0     ✔ forcats 0.5.0
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
library(ggrepel)
```

### Loading and exploring prebuilt datasets

Load the gwas\_CD dataset, which is a subset of association results
(SNPs with P&lt;1e-03) for Crohn´s disease from the UK biobank.

It is highly recommended to theck the number of datapoints in your
dataset before you plot, since a very large dataset will take a long
time to plot.

``` r
paste("Number of SNPs in the dataset: [", length(CD_UKBB$POS),"]", sep = "")
#> [1] "Number of SNPs in the dataset: [26821]"
```

### Manhattan plots

Get an overview of association results for crohn’s disease (CD) in a
Manhattan plot

``` r
manhattan(CD_UKBB)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="120%" />

### QQ plots

``` r
qqtopr(CD_UKBB,n_variants=length(CD_UKBB$POS))
```

### Label the top SNPs with the name of their nearest gene

Use the annotate argument in the manhattan function to label the top
SNPs with p-values below the annotate threshold with their nearest gene

``` r
manhattan(CD_UKBB, annotate = 1e-09)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="120%" />

## Highlight genes of interest

``` r
manhattan(CD_UKBB, annotate = 1e-09, highlight_genes = c("NOD2","IL23R","JAK2"))
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="120%" />

### View one chromsome only

Take a closer look at the results by chromosome. Here we plot the
results on chromosome 7 only.

``` r
manhattan(CD_UKBB, annotate = 1e-09, chr = "7")
#> [1] "7"
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="120%" />

### Regionplot

Zoom in further on the chromosome plot with the regionplot function.

Zoom in on a gene of interest, e.g IKZF1:

``` r
regionplot(CD_UKBB, gene="IKZF1", annotate= 1e-09)
#> [1] "7"
#> [1] "Zoomed to region:  chr7:50204067-50505101"
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="120%" />

Zoom in on the top hit on a chromosome

``` r
CHR <- "chr1"
top_hit <- get_top_hit(CD_UKBB,chr=CHR)
regionplot(CD_UKBB, chr = CHR, xmin=top_hit$POS-250000 ,xmax= top_hit$POS+250000)
#> [1] "1"
#> [1] "Zoomed to region:  chr1:66966513-67466513"
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="120%" />

### Display multiple phenotypes/datasets on the same plot

Display the output from more than one GWAS on the same plot

### Manhattan multiple phenotypes

``` r
manhattan(list(CD_UKBB,CD_FINNGEN,UC_UKBB),legend_labels = c("CD UKBB","CD Finngen","UC UKBB"),title="IBD")
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="120%" />

``` r
regionplot(list(CD_UKBB,CD_FINNGEN),gene="NOD2", annotate=c(1e-15, 1e-09), legend_labels = c("UKBB", "Finngen"), title="Crohn's disease (CD)")
#> [1] "16"
#> [1] "16"
#> [1] "Zoomed to region:  chr16:50593587-50834041"
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="120%" />

``` r
manhattan(list(CD_UKBB,CD_FINNGEN,UC_UKBB),legend_labels = c("CD UKBB","CD Finngen","UC UKBB"),title="IBD")
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="120%" />

#### The ntop argument

Use the ntop argument to set the number of datasets displayed at the top
(default value is 3)

``` r
manhattan(list(CD_UKBB,CD_FINNGEN,UC_UKBB),ntop=2,legend_labels = c("CD UKBB","CD Finngen","UC UKBB"),title="IBD")
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="120%" />

## Useful functions

``` r
get_top_hit(CD_UKBB, chr="chr16")
dat1 <- get_best_snp_per_MB(CD_UKBB,thresh = 1e-07, region=1000000)
#get overlapping SNPS overlapping in two datasets
overlapping_snps <- dat1 %>% get_overlapping_snps_by_pos(CD_FINNGEN)

overlapping_snps_matched <- overlapping_snps %>% match_alleles()
overl_snps_matched_pos_allele_dat1 <- overlapping_snps_matched %>% flip_to_positive_allele_for_dat1()
snpset1 <- overl_snps_matched_pos_allele_dat1 %>% annotate_with_nearest_gene(protein_coding_only = T)

#or do all this in one go, by calling the create_snpset functon
snpset1 <- create_snpset(CD_FINNGEN, CD_UKBB, thresh = 1e-06)
snpset2 <- create_snpset(CD_UKBB, CD_FINNGEN, thresh= 1e-06)
```

``` r
e1 <- effect_plot(snpset1, pheno_x="CD Finngen", pheno_y="CD UKBB",color=get_topr_colors()[1], gene_label_thresh = 1)
e2 <- effect_plot(snpset2, pheno_x="CD UKBB", pheno_y="CD Finngen", color=get_topr_colors()[2], gene_label_thresh=1,annotate_with = "ID")
grid.arrange(e1,e2)
```

### Plot apperance: setting text sizes

``` r
manhattan(list(CD_UKBB,CD_FINNGEN), annotate=1e-09,axis_title_size = 20,axis_text_size = 16,label_size = 5, title_text_size = 16, legend_text_size = 20)
#> [1] "Use the legend_labels argument to change the legend labels from color names to meaningful labels! "
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="120%" />

``` r
regionplot(list(CD_UKBB,CD_FINNGEN), gene="IKZF1", vline=50274703,title="CD UKBB", title_text_size = 16,axis_title_size = 20,axis_text_size = 20,legend_text_size = 20)
#> [1] "7"
#> [1] "7"
#> [1] "Use the legend_labels argument to change the legend labels from color names to meaningful labels! "
#> Warning in min(dat[[i]]$log10p): no non-missing arguments to min; returning Inf
#> Warning in max(dat[[i]]$log10p): no non-missing arguments to max; returning -Inf
#> [1] "Zoomed to region:  chr7:50204067-50505101"
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="120%" />

Setting alpha, size and shape
