
#' Match the variants in the snpset by their alleles
#'
#' @description
#'
#' \code{match_alleles()}
#'
#' This method is deprecated and will be removed in future versions. use \code{\link{match_by_alleles}} instead.
#'
#' @param df A dataframe that is in the snpset format (like returned by the get_snpset() function)
#' @param verbose A logical scalar (default: FALSE). Assign to TRUE to get information on which alleles are matched and which are not.
#'
#' @return The input dataframe containing only those variants whith matched alleles in the snpset
#' @export
#'
#' @examples
#' \dontrun{
#' match_alleles(df)
#' }
#'

match_alleles <- function(df, verbose=F){
  .Deprecated("match_by_alleles")
  df$ID_tmp <- paste(df$CHROM, df$POS, sep="_")
  matched_snps <- df %>% dplyr::filter(REF1 == REF2 & ALT1 == ALT2)
  #check whether ref and alt are reversed
  if(length(matched_snps$POS) != length(df$POS)){
    swapped_allele <-  df %>% dplyr::filter(REF1 == ALT2 & REF2 == ALT1)
    swapped_allele$E2 <- swapped_allele$E2*(-1)
    swapped_allele$REF2tmp <- swapped_allele$REF2
    swapped_allele$REF2 <- swapped_allele$ALT2
    swapped_allele$ALT2 <- swapped_allele$REF2tmp
    swapped_allele <- swapped_allele %>% dplyr::select(-REF2tmp)
    matched_snps <- rbind(matched_snps, swapped_allele)
  }
  #check whether all variants were mached
  not_matched <- df %>% dplyr::filter(! ID_tmp %in% matched_snps$ID_tmp)
  if(verbose){
    if(length(not_matched$POS)>0){
      print(paste("Could not match ",length(not_matched$POS), " snps on their alleles."))
      print(not_matched)
      print(paste(length(matched_snps$POS), " SNPs have matching REF and ALT alleles. ", sep=""))
    }else{
      print("All SNPs were matched by their REF and ALT alleles")
    }
  }
  return(matched_snps %>% dplyr::select(-ID_tmp))
}

#' Get variants that overlap between two datasets
#'
#' @description
#'
#' \code{get_overlapping_snps_by_pos()}
#' 
#' This method is deprecated and will be removed in future versions. use \code{\link{match_by_pos}} instead.
#'
#' @param df1 A dataframe of variants, has to contain CHROM and POS
#' @param df2 A dataframe of variants, has to contain CHROM and POS
#' @param verbose A logical scalar (default: FALSE). Assign to TRUE to get information on which alleles are matched and which are not.
#'
#'
#' @return The input dataframe containing only those variants with matched alleles in the snpset
#' @export
#'
#' @examples
#' \dontrun{
#' get_overlapping_snps_by_pos(dat1, dat2)
#' }
#'
get_overlapping_snps_by_pos <- function(df1, df2,verbose=F){
  .Deprecated("match_by_pos")
  dat2 <- dat_check(df2)
  dat1 <- dat_check(df1)
  df1 <- dat1[[1]]
  df1_orig <- df1
  df2 <- dat2[[1]]
  df2$ID_tmp <- paste(df2$CHROM, df2$POS, sep="_")
  df1$ID_tmp <- paste(df1$CHROM, df1$POS, sep="_")
  if((! "BETA" %in% colnames(df2)) & ("OR" %in% colnames(df2))){
    df2$BETA <- log(df2$OR)
  }
  if((! "BETA" %in% colnames(df1)) & ("OR" %in% colnames(df1))){
    df1$BETA <- log(df1$OR)
  }
  if((! "REF" %in% colnames(df1)) & (! "REF" %in% colnames(df2))){
    warning("The input datasets have to include REF and ALT columns to be able to use this function")
  }
  df2 <- df2 %>% dplyr::select(CHROM,POS,REF,ALT,P,BETA) %>% dplyr::rename(P2=P,E2=BETA,ALT2=ALT,REF2=REF)
  if("Gene_Symbol" %in% colnames(df1)){
    df1 <- df1 %>% dplyr::select(CHROM,POS,REF,ALT,P,BETA,Gene_Symbol,ID_tmp) %>% dplyr::rename(P1=P,E1=BETA,ALT1=ALT,REF1=REF)
  }
  else{
    df1 <- df1 %>% dplyr::select(CHROM,POS,REF,ALT,P,BETA, ID_tmp ) %>% dplyr::rename(P1=P,E1=BETA,ALT1=ALT,REF1=REF)
  }
  snpset <- df1 %>% dplyr::inner_join(df2,by=c("CHROM","POS")) %>% dplyr::arrange(CHROM,POS,P1) %>% dplyr::distinct(CHROM,POS, .keep_all = TRUE)
  if("ID" %in% colnames(df1_orig)){
    #add the ID column back to the dataframe
    snpset <- snpset %>% dplyr::inner_join(df1_orig %>% dplyr::select(CHROM,POS,ID), by=c("CHROM","POS")) %>% dplyr::distinct(CHROM,POS, .keep_all = TRUE)
  }
  else{
    print("ID is not in colnames df1")
    print(colnames(df1))
  }
  if(verbose){
    print("Overlapping SNPs: ")
    not_found <- df1 %>% dplyr::filter(! ID_tmp %in% snpset$ID_tmp)
    print(paste("There are a total of ",length(df1$POS), " SNPs in the first dataset. Thereof [",length(snpset$POS), "] were are also the second dataset (dat2) and [", length(not_found$POS), "] are not.", sep=""))
    print(paste("SNPs FOUND in dat2: ",length(snpset$POS), sep="" ))
    print(snpset)
    if(length(not_found$POS)>0){
      print(paste("SNPs NOT FOUND in dat2: " ,length(not_found$POS), sep=""))
      print(not_found)
    }
  }
  return(snpset %>% dplyr::select(-ID_tmp))
}

#' Create a dataframe that can be used as input for making effect plots
#'
#' @description
#'
#' \code{create_snpset()}
#' 
#' This method is deprecated and will be removed in future versions. use \code{\link{get_snpset}} instead.
#'
#' @param df1 The dataframe to extract the top snps from (with p-value below thresh)
#' @param df2 The dataframe in which to search for overlapping SNPs from dataframe1
#' @param thresh Numeric, the p-value threshold used for extracting the top snps from dataset 1
#' @param region_size Integer, the size of the interval which to extract the top snps from
#' @param protein_coding_only Logical, set this variable to TRUE to only use protein_coding genes for the annotation
#' @param verbose Logical, (default: FALSE). Assign to TRUE to get information on which alleles are matched and which are not.

#' @return Dataframe containing the top hit
#' @export
#'
#' @examples
#' \dontrun{
#' create_snpset(CD_UKBB,CD_FINNGEN, thresh=1e-09)
#' }
#'

create_snpset <- function(df1, df2, thresh=1e-08,protein_coding_only=TRUE, region_size=1000000, verbose=F){
  .Deprecated("get_snpset")
  snpset <- df1 %>% get_best_snp_per_MB(thresh = thresh,region_size=region_size) %>% get_overlapping_snps_by_pos(df2,verbose=verbose) %>% match_alleles(verbose=verbose) %>% flip_to_positive_allele_for_dat1()
  if(! "Gene_Symbol" %in% colnames(df1)  ||  ! "gene_symbol" %in% colnames(df1)){
     snpset <- snpset %>% annotate_with_nearest_gene(protein_coding_only = protein_coding_only)
  }
  return(snpset)
}

#' Show the code/functions used to create a snpset
#'
#' @description
#' 
#' This method is deprecated and will be removed in future versions. use \code{\link{get_snpset_code}} instead.
#' 
#' \code{create_snpset_code()}
#'
#' @return Dataframe containing the top hit
#' @export
#'
#' @examples
#' \dontrun{
#' create_snpset_code()
#' }
#'
create_snpset_code <-function(){
  .Deprecated("get_snpset_code")
  get_snpset_code()
}


#' Create a plot comparing effects within two datasets
#'
#' @description
#'
#' \code{effect_plot()}
#' 
#' This method is deprecated and will be removed in future versions. use \code{\link{effectplot}} instead.
#'
#' @param dat The input dataframe (snpset) containing one row per variant and P values (P1 and P2) and effects (E1 and E2) from two datasets/phenotypes
#' @param pheno_x A string representing the name of the phenotype whose effect is plotted on the x axis
#' @param pheno_y A string representing the name of the phenotype whose effect is plotted on the y axis
#' @param annotate_with  A string, The name of the column that contains the label for the datapoints (default value is Gene_Symbol)
#' @param thresh A number. Threshold cutoff, datapoints with P2 below this threshold are shown as filled circles whereas datapoints with P2 above this threshold are shown as open circles
#' @param ci_thresh A number.Show the confidence intervals if the P-value is below this threshold
#' @param gene_label_thresh A string, label datapoints with P2 below this threshold
#' @param color A string, default value is the first of the topr colors
#' @param scale A number, to change the size of the title and axes labels and ticks at the same time (default = 1)
#' @export
#'
#' @examples
#' \dontrun{
#' effect_plot(dat)
#' }
#'
#'
effect_plot <- function(dat,pheno_x="pheno_x", pheno_y="pheno_", annotate_with="Gene_Symbol", thresh=1e-08, ci_thresh=1,gene_label_thresh = 1e-08, color=get_topr_colors()[1],scale=1){
  .Deprecated("effectplot")
  return(effectplot(dat, pheno_x=pheno_x, pheno_y=pheno_y, annotate_with=annotate_with, thresh=thresh, ci_thresh=ci_thresh, gene_label_thresh = gene_label_thresh, color=color, scale=scale))
}


#' Get the genetic position of a gene by gene name
#'
#' @description
#'
#' \code{get_gene()} Get the gene coordinates for a gene
#' Required parameter is gene name
#' 
#' This method is deprecated and will be removed in future versions. use \code{\link{get_gene_coords}} instead.
#'
#' @param gene_name A string representing a gene name (e.g. "FTO")
#' @param chr A string, search for the genes on this chromosome only, (e.g chr="chr1")
#' @param build A string, genome build, choose between builds 37 (GRCh37) and 38 (GRCh38) (default is 38)
#' @return Dataframe with the gene name and its genetic coordinates
#' @export
#'

get_gene <- function(gene_name,chr=NULL, build=38){
  .Deprecated("get_gene")
  return(get_gene_coords(gene_name=gene_name,chr=chr,build=build))
}


#' Get the index/lead variants
#'
#' @description
#'
#' \code{get_best_snp_per_MB()} Get the top variants within 1 MB windows of the genome with association p-values below the given threshold 
#'
#' This method is deprecated and will be removed in future versions. use \code{\link{get_lead_snps}} instead.
#' 
#' @param df Dataframe
#' @param thresh A number. P-value threshold, only extract variants with p-values below this threshold (5e-09 by default)
#' @param protein_coding_only Logical, set this variable to TRUE to only use protein_coding genes for annotation
#' @param chr String, get the top variants from one chromosome only, e.g. chr="chr1"
#' @param .checked Logical, if the input data has already been checked, this can be set to TRUE so it wont be checked again (FALSE by default)
#' @param verbose Logical, set to TRUE to get printed information on number of SNPs extracted
#' @return Dataframe of lead variants. Returns the best variant per MB (by default, change the region size with the region argument) with p-values below the input threshold (thresh=5e-09 by default)
#' @export
#' @inheritParams regionplot


get_best_snp_per_MB <- function(df, thresh=5e-09,region_size=1000000,protein_coding_only=FALSE,chr=NULL, .checked=FALSE, verbose=FALSE){
  .Deprecated("get_lead_snps")
  return(get_lead_snps(df=df, thresh=thresh,region_size=region_size,protein_coding_only=protein_coding_only,chr=chr, .checked=.checked, verbose=verbose))
}


# Temporary deprecated argument handler
# msg of form "[arg_name] argument deprecated [- optional custom message]"

deprecated_argument_msg <- function(arg, custom=NULL) {
  deparse(substitute(arg)) %>%
    paste("argument deprecated", if (!is.null(custom)) paste("-", custom))
}
