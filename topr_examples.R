
library(gorr)

setwd("/Users/thorhildur/IdeaProjects/topr")
source("/Users/thorhildur/IdeaProjects/topr/assoc_plots_v1.R")


#Get the input data
gwas_CD=read.table(file="inst/extdata/crohns_K50_UKBB_p1e-03.tsv", as.is=T,header=T)

#Start by checking how many variants are in the datasets, since we do not want to plot with too many datapoints
paste("Number of SNPs in rhw dataset: [",length(gwas_CD$POS),"]", sep="")

## MANHATTAN
#Get an overview of the crohn's disease (CD) results in a Manhattan plot
manhattan(gwas_CD)

#get the top SNP per 10 MB
snps_CD=get_best_snp_per_MB(gwas_CD,thresh=1e-09,region=10000000)

manhattan(gwas_CD, variants=snps_CD,annotation_thresh = 1e-09)

## CHROMPLOT
#Take a closer look at the results by chromosome
CHR="chr16"
chromplot(gwas_CD,chr=CHR, variants=snps_CD,annotation_thresh=5e-09)

## REGIONPLOT
# Zoom into a 1 MB region centered on the top variant on the chromosome
top_snp=gwas_CD %>% filter(CHROM == CHR) %>% arrange(P) %>% head(n=1)
#set the 1 MB range centered on the top snp
xmin=top_snp$POS-500000
xmax=top_snp$POS+500000
xmin=ifelse(xmin<0, 0,xmin)
#currently the regionplot uses gor to get the genes, so platform connection is needed
#regionplot(gwas_results,chr=CHR,xmin=xmin,xmax=xmax,variants=lead_snps,annotation_thresh=1e-08,vline=top_snp$POS)

## MANHATTAN MULTI
# Display the output from more than one GWAS on the same manhattan plot

gwas_CD_IIBD=read.table(file="inst/extdata/EUR.CD.gwas_info03_filtered.assoc.cutoff_p1e_03_hg38_filt.gor", as.is=T,header=T)
manhattan(list(gwas_CD_IIBD,gwas_CD),legend_labels=c("Crohns IIBD", "Crohns UKBB"),alpha=c(1,0.5))
