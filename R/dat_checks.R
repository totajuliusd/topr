dat_chr_check=function(dat,chr=NULL){
  for(i in 1:length(dat)){
    df=as.data.frame(dat[[i]])
    #remove chr from CHROM if its there, and set chrX to 23
    df=df %>% dplyr::mutate(CHROM=gsub('chr','',CHROM))
    df[df$CHROM=='X', 'CHROM']="23"
    df$CHROM = as.integer(df$CHROM)
    df = df %>% filter(CHROM<24)
    df$POS=as.integer(df$POS)
    df=df %>% arrange(CHROM,POS)
    dat[[i]]=df
  }
  return(dat)
}


#' dat_column_check_and_set
#'
#' @description
#'
#' This function is used to standardize the column names in the input dataframe
#'
#'
#' @param dat  The data input, a data frame or a list of data frames
dat_column_check_and_set=function(dat){
  for(i in 1:length(dat)){
    df=as.data.frame(dat[[i]])
    if("pos" %in% colnames(df)) df=df %>% rename(POS=pos)
    if("BP" %in% colnames(df)) df=df %>% rename(POS=BP)
    if("chrom" %in% colnames(df)) df=df %>% rename(CHROM=chrom)
    if("CHR" %in% colnames(df)) df=df %>% rename(CHROM=CHR)
    if("chr" %in% colnames(df)) df=df %in% rename(CHROM=chr)

    if("rsid" %in% colnames(df)) df=df %>% rename(ID=rsid)
    if("rsId" %in% colnames(df)) df=df %>% rename(ID=rsId)
    if("rsName" %in% colnames(df)) df=df %>% rename(ID=rsName)
    if("SNP" %in% colnames(df)) df=df %>% rename(ID=SNP)
    if("snp" %in% colnames(df)) df=df %>% rename(ID=snp)

    if("gene_symbol" %in% colnames(df)) df=df %>% rename("Gene_Symbol" = "gene_symbol")
    if("gene" %in% colnames(df)) df=df %>% rename("Gene_Symbol" = "gene")
    if("geneName" %in% colnames(df)) df=df %>% rename("Gene_Symbol" = "geneName")

    if("pval" %in% colnames(df)) df=df %>% rename("P" = "pval")
    if("p" %in% colnames(df)) df=df %>% rename("P" = "p")
    if("pvalue" %in% colnames(df)) df=df %>% rename("P" = "pvalue")

    if(! "ID" %in% colnames(df)){
      df$ID=paste(df$CHROM, df$POS, sep="_")
    }
    if(! "Gene_Symbol" %in% colnames(df)) df$Gene_Symbol=df$ID
    if("Gene_Symbol" %in% colnames(df)){
      df$Gene_Symbol=sub(",.*","", df$Gene_Symbol)
      df$Gene_Symbol=ifelse(df$Gene_Symbol == ".",df$ID, df$Gene_Symbol)
    }
    if("max_impact" %in% colnames(df)) df=df %>% rename("Max_Impact" = "max_impact")
    if(! "Max_Impact" %in% colnames(df)){
      if("Impact" %in% colnames(df)) df=df %>% rename("Max_Impact" = "Impact")
      if("impact" %in% colnames(df)) df=df %>% rename("Max_Impact" = "impact")
    }
    if( (!"CHROM" %in% colnames(df)) || (!"POS" %in% colnames(df)) || (!"P" %in% colnames(df) ))
      stop("Some columns are missing from the dataset. Required columns are CHROM,POS and P. Add the required columns and try again, or rename existing columns, e.g. df=df %>% dplyr::rename(CHROM=yourColname)")
    dat[[i]]=df
  }
  return(dat)
}
