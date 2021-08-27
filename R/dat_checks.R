
dat_check <- function(dat){
  if(is.data.frame(dat)) dat <- list(dat)
  dat <- dat_column_check_and_set(dat) %>% dat_chr_check()
  return(dat)
}

dat_chr_check <- function(dat){
  for(i in seq_along(dat)){
    df <- as.data.frame(dat[[i]])
    #remove chr from CHROM if its there, and set chrX to 23
    df <- df %>% dplyr::mutate(CHROM=gsub('chr','',CHROM))
    df[df$CHROM=='X', 'CHROM'] <- "23"
    df$CHROM <-  as.integer(df$CHROM)
    df <- df %>% dplyr::filter(CHROM<24)
    df$POS <- as.integer(df$POS)
    df <- df %>% dplyr::arrange(CHROM,POS)
    dat[[i]] <- df
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
#'

dat_column_check_and_set <- function(dat){
  for(i in seq_along(dat)){
    df <- as.data.frame(dat[[i]])
    if("pos" %in% colnames(df)) df <- df %>%  dplyr::rename(POS=pos)
    else if("BP" %in% colnames(df)) df <- df %>%  dplyr::rename(POS=BP)

    if("chrom" %in% colnames(df)) df <- df %>%  dplyr::rename(CHROM=chrom)
    else if("CHR" %in% colnames(df)) df <- df %>%  dplyr::rename(CHROM=CHR)
    else if("chr" %in% colnames(df)) df <- df %in%  dplyr::rename(CHROM=chr)

    if("rsid" %in% colnames(df)) df <- df %>%  dplyr::rename(ID=rsid)
    else if("rsId" %in% colnames(df)) df <- df %>%  dplyr::rename(ID=rsId)
    else if("rsName" %in% colnames(df)) df <- df %>%  dplyr::rename(ID=rsName)
    else if("SNP" %in% colnames(df)) df <- df %>%  dplyr::rename(ID=SNP)
    else if("snp" %in% colnames(df)) df <- df %>%  dplyr::rename(ID=snp)

    if("gene_symbol" %in% colnames(df)) df <- df %>%  dplyr::rename("Gene_Symbol" = "gene_symbol")
    else if("gene" %in% colnames(df)) df <- df %>%  dplyr::rename("Gene_Symbol" = "gene")
    else if("geneName" %in% colnames(df)) df <- df %>%  dplyr::rename("Gene_Symbol" = "geneName")

    if(! "P" %in% colnames(df)){
      if("pval" %in% colnames(df)) df <- df %>%  dplyr::rename("P" = "pval")
      else if("p" %in% colnames(df)) df <- df %>%  dplyr::rename("P" = "p")
      else if("pvalue" %in% colnames(df)) df <- df %>%  dplyr::rename("P" = "pvalue")
    }
    if(! "ID" %in% colnames(df)){
      if("CHROM" %in% colnames(df)){
        df$ID <- paste(df$CHROM, df$POS, sep="_")
      }
      else{
        stop("Could not find the chromosome column (CHROM) in the input data")
      }
    }
   # if(! "Gene_Symbol" %in% colnames(df)) df$Gene_Symbol=df$ID
    if("Gene_Symbol" %in% colnames(df)){
      df$Gene_Symbol <- sub(",.*","", df$Gene_Symbol)
      df$Gene_Symbol <- ifelse(df$Gene_Symbol == ".",df$ID, df$Gene_Symbol)
    }
    if("max_impact" %in% colnames(df)) df <- df %>%  dplyr::rename("Max_Impact" = "max_impact")
    if(! "Max_Impact" %in% colnames(df)){
      if("Impact" %in% colnames(df)) df <- df %>%  dplyr::rename("Max_Impact" = "Impact")
      if("impact" %in% colnames(df)) df <- df %>%  dplyr::rename("Max_Impact" = "impact")
    }
    #only needed for the effec plot
    if("beta" %in% colnames(df)) df <- df %>%  dplyr::rename("BETA" = "beta")
    if("ref" %in% colnames(df)) df <- df %>%  dplyr::rename("REF" = "ref")
    if("alt" %in% colnames(df)) df <- df %>%  dplyr::rename("ALT" = "alt")
    if( (!"CHROM" %in% colnames(df)) || (!"POS" %in% colnames(df)) || (!"P" %in% colnames(df) ))
      stop("Some columns are missing from the dataset. Required columns are CHROM,POS and P. Add the required columns and try again, or rename existing columns, e.g. df=df %>% dplyr::rename(CHROM=yourColname)")
    dat[[i]] <- df
  }

  return(dat)
}
