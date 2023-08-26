
# *topr*: an R package for viewing and annotating genetic association results

<img src="man/figures/manhattan_wide_1080_high.gif" alt="topr GIF" width="100%">

## Citation

Please cite the following paper if you use *topr* in a publication:

Juliusdottir, T. *topr*: an R package for viewing and annotating genetic association results. BMC Bioinformatics 24, 268 (2023). https://doi.org/10.1186/s12859-023-05301-4

## Installation
<hr>

Install from CRAN:

``` r
install.packages("topr")
```

Or from github:

``` r
devtools::install_github("totajuliusd/topr")
```

And then load the package:

``` r
library(topr)
```


## Main features and functionality
<hr>

*topr* is written in the R programming language and utilises the ggplot2 and ggrepel R graphics libraries for plotting.

*topr's* two main plotting functions are <code>manhattan()</code> and <code>regionplot()</code>. 

The manhattan() function returns a **ggplot object**. The regionplot() function draws three ggplotGrobs aligned using egg::gtable_frame, however it can be called with the <code>extract_plots</code> argument set to TRUE to return a list of three **ggplot objects** instead.


### Example input datasets 

*See the <a href="https://totajuliusd.github.io/topr/articles/input_datasets.html">Input datasets vignette</a> for more detailed information.*

Input datasets must include least three columns (<code>CHROM, POS</code> and <code>P</code>), where naming of the columns is flexible (i.e the chr label can be either chr or chrom and is case insensitive).

*topr* has 3 inbuilt GWAS results (<code>CD_UKBB, CD_FINNGEN</code> and <code>UC_UKBB</code>). To get information on them, do:

``` r
?CD_UKBB
?CD_FINNGEN
?UC_UKBB
```

The chromosome in the <code>CHROM</code> column can be represented with or without the <i>chr</i> suffix, e.g (chr1 or 1)


### Basic usage 
<hr>

Basic usage of *topr's* key functions is as follows.
<br>

#### Manhattan
<hr>

See the <a href="https:///totajuliusd.github.io/topr/articles/manhattan.html">Manhattan vignette</a> for more detailed examples of how to use the manhattan plot function.

View the whole genome association results on a Manhattan plot:

``` r
manhattan(CD_UKBB)
```

Annotate the lead/index variants (with p-values below 5e-9) with their nearest gene:

``` r
manhattan(CD_UKBB, annotate=5e-9)
```

View multiple GWAS results on the same plot

``` r
manhattan(list(CD_UKBB, CD_FINNGEN), legend_labels = c("UKBB", FinnGen"))
```

<br>

#### Regionplot
<hr>
See the <a href="https:///totajuliusd.github.io/topr/articles/regionplot.html">Regionplot vignette</a> for more detailed examples of how to use the regionplot function.

Further zoom-in on a genetic region by gene name (*IL23R*):

``` r
regionplot(CD_UKBB, gene="IL23R")
```

View the correlation pattern between the variants within the region in a locuszoom like plot.
Note that the variant correlation (<code>R2</code>) has to be pre-calculated and included in the input dataframe.

``` r
locuszoom(R2_CD_UKBB)
```

Display multiple GWAS results zoomed in on the *IL23R* gene

``` r
regionplot(list(UC_UKBB, CD_UKBB), gene="IL23R")
```

<br>

#### Other useful functions
<hr>

Extract lead/index variants from the GWAS dataset (<code>CD_UKBB</code>):

```{r}
get_lead_snps(CD_UKBB)
```

Annotate the lead/index variants with their nearest gene:

```{r}
get_lead_snps(CD_UKBB) %>% annotate_with_nearest_gene()
```

Get genomic coordinates for a gene (*topr* uses genome build GRCh38.p13 by default):

```{r}
get_gene_coords("IL23R")
```
Get genomic coordinates for a gene using genome build GRCh37 instead.

```{r}
get_gene_coords("IL23R", build="37")
```

Get snps within a region:

```{r}
get_snps_within_region(CD_UKBB, region = "chr1:67138906-67259979")
```

Get the top variant on a chromosome:
```{r}
get_top_snp(CD_UKBB, chr="chr1")
```

Create a snpset by extracting the top/lead variants from the CD_UKBB dataset and overlapping variants (by position) in the CD_FINNGEN dataset. 

```{r}
get_snpset(CD_UKBB, CD_FINNGEN)
```

Create an effecplot by plotting the effect sizes of top/lead variants overlapping in two datasets.
```{r}
effectplot(list(CD_UKBB, CD_FINNGEN), annotate = 1e-08)
```

<br>

#### Get help:

``` r
?manhattan()
?regionplot()
?locuszoom()
```
