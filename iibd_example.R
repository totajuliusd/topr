
library(gorr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library("egg")
library("grid")
#devtools::install_github("wuxi-nextcode/gorr")

setwd("/Users/thorhildur/IdeaProjects/topR")
source("/Users/thorhildur/IdeaProjects/topR/assoc_plots_v1.R")

#sessionInfo()

#https://platform.wuxinextcodedev.com/
#loggaðu þig inn, og þar neðst er hlekku á api-key-service
#https://platform.wuxinextcodedev.com/api-key-service/token
api_key <- "eyJhbGciOiJIUzI1NiIsInR5cCIgOiAiSldUIiwia2lkIiA6ICI1YzkxMjYyOC01NGRkLTQ3MTctODRmNi04NDM3ZTcyMDIyMzQifQ.eyJpYXQiOjE2MDYzMzg2MzgsImp0aSI6ImUzZjYwMDE4LWMzM2YtNGM2Mi05NDYwLThhOTQ5NjBhOTFkYiIsImlzcyI6Imh0dHBzOi8vcGxhdGZvcm0ud3V4aW5leHRjb2RlZGV2LmNvbS9hdXRoL3JlYWxtcy93dXhpbmV4dGNvZGUuY29tIiwiYXVkIjoiaHR0cHM6Ly9wbGF0Zm9ybS53dXhpbmV4dGNvZGVkZXYuY29tL2F1dGgvcmVhbG1zL3d1eGluZXh0Y29kZS5jb20iLCJzdWIiOiJlYmEwYjQ4OS1mZGQxLTQ3YjgtOGMxNy04OWFjMjg4YWM1YWQiLCJ0eXAiOiJPZmZsaW5lIiwiYXpwIjoiYXBpLWtleS1jbGllbnQiLCJub25jZSI6IjI1MWZlYmY4LTIxMGUtNGVhNC04NTczLTM2NTU1ODk3NDE2NCIsInNlc3Npb25fc3RhdGUiOiJmMjEyMTE1Mi01YzdiLTQ0YWEtOTgyZC1hZWVhN2QzYzAxZDkiLCJzY29wZSI6Im9wZW5pZCBvZmZsaW5lX2FjY2VzcyJ9.6H6-k6wBzZnRl99MCssbxfXpt-MtyE6gaPEMO266zb0"
project = "ukbb_hg38"
conn <- gor_connect(api_key, project)

query <- function(query) {gor_query(query, conn)} # Defining a simple function so we don't need to re-enter the connection every time we run a query

CD_f="/Users/thorhildur/work/genuityscience/public_data/EUR.CD.gwas_info03_filtered.assoc.cutoff_p_0.05.tsv"
UC_f="/Users/thorhildur/work/genuityscience/public_data/EUR.UC.gwas_info03_filtered.assoc.cutoff_p_0.05.tsv"
IBD_f="/Users/thorhildur/work/genuityscience/public_data/EUR.IBD.gwas_info03_filtered.assoc.cutoff_p_0.05.tsv"

crohns=read.table(CD_f,as.is=T,header=T)
uc=read.table(UC_f, as.is=T,header=T)
ibd=read.table(IBD_f,as.is=T,header=T)

install.packages("showtext")


CD_filt=crohns %>% filter(P<1e-02)
UC_filt=uc %>% filter(P<1e-02)
IBD_filt=ibd %>% filter(P<1e-02)
CD_filt=thin_data(CD_filt)
UC_filt=thin_data(UC_filt)
IBD_filt=thin_data(IBD_filt)

color=c("#E69F00","darkblue","#00AFBB","#FC4E07","darkorange1") # "#E69F00",
manhattan(list(CD_filt,UC_filt,IBD_filt),size=c(1,0.6, 0.6),color=color,legend_labels=c("CD","UC","IBD"),ntop=2,alpha=c(1,0.5,1))



CD_snps=get_best_snp_per_MB(CD_filt)
UC_snps=get_best_snp_per_MB(UC_filt)
IBD_snps=get_best_snp_per_MB(IBD_filt)

length(CD_snps$POS)
length(UC_snps$POS)
length(IBD_snps$POS)
manhattan(UC_filt)
manhattan(CD_filt)
#"#999999"

colnames(UC_filt)

dat=CD_filt %>% filter(P<1e-09)

if(is.data.frame(dat)) dat=list(dat)
dat=dat_column_check_and_set(dat)
dat=dat_chr_check(dat)
incl_chrX=include_chrX(dat)
offsets=get_chr_offsets(incl_chrX)
for(i in 1:length(dat)){
  dat[[i]]=dat[[i]] %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
}

dat[[1]]$tmp=NA

if(annotate_with_lead_snps){
  lead_snps=get_best_per_MB_snps(dat)
}
CD_filt=crohns %>% filter(P<1e-03)
dat=CD_filt %>% filter(P<1e-06)
dat
lead_snps=get_best_per_MB_snps(dat)

dat
length(lead_snps$POS)

get_best_per_MB_snps=function(df,thresh=1e-09){
    if(! "CHROM" %in% colnames(df) || (! "POS" %in% colnames(df)) || (! "P" %in% colnames(df))){
       if(is.data.frame(dat)) dat=list(df)
       dat=dat_column_check_and_set(dat)
       dat=dat_chr_check(dat)
       df=dat[[1]]
    }
    df=df %>% filter(P<thresh)
    df$tmp=NA
    for(row in 1:nrow(df)){
      df$tmp=round(df$POS/1000000)
    }
    lead_snps=df %>% group_by(CHROM,tmp) %>% arrange(P) %>% filter(P== min(P))
    if(length(lead_snps)== 0){
      print(paste("There are no SNPs with a p-value below ",thresh, " in the input. Use the [thresh] argument to adjust the threshold.",sep=""))
    }
    return(lead_snps)
}
