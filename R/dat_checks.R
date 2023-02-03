
dat_check <- function(dat, verbose=TRUE,locuszoomplot=FALSE){
  if(is.data.frame(dat)) dat <- list(dat)
  for(i in seq_along(dat)){
    if(nrow(dat[[i]])==0){
      warning(paste0("Dataset [",i, "] is empty with 0 rows. Please remove it from the input list and re-run!"))
    }
  }
  dat <- dat_column_check_and_set(dat, verbose=verbose,locuszoomplot=locuszoomplot) %>% dat_chr_check()
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

convert_region_size <- function(region_size){
  region_size <- gsub("(\\d+)(kb|KB|Kb)","\\1000", region_size)
  region_size <- gsub("(\\d+)(mb|MB|Mb)","\\1000000", region_size)
  if(! is.numeric(region_size)){
    region_size <- as.numeric(region_size)
    
  }
  return(region_size)
}

dat_column_chr_pos_check <- function(df){
  df <- df %>% 
    rename_with(~ rename_value(.x, "CHROM"), matches(c("^chrom$","^chr$","^chromosome$"), ignore.case = TRUE)) %>% 
    rename_with(~ rename_value(.x, "POS"), matches(c("^pos$","^bp$","^base_pair_location$"), ignore.case = TRUE))
  return(df)
}

dat_column_check_and_set <- function(dat, verbose=TRUE,locuszoomplot=FALSE){
  for(i in seq_along(dat)){
    df <- as.data.frame(dat[[i]])
    df<- dat_column_chr_pos_check(df)
    if(locuszoomplot)
      dfwithcolor=df
    df <- df %>% 
      dplyr::rename_with(~ rename_value(.x, "ID"), matches(c("^rsid$","^rsname$","^snp$"), ignore.case = TRUE)) %>% 
      dplyr::rename_with(~ rename_value(.x, "Gene_Symbol"), matches(c("^gene_symbol$","^gene$","^genename$","^gene_name$"), ignore.case = TRUE)) %>%
      dplyr::rename_with(~ rename_value(.x, "Max_Impact"), matches(c("^max_impact$","^impact$"), ignore.case = TRUE)) %>% 
      dplyr::rename_with(~ rename_value(.x, "OR"), matches(c("^odds_ratio$","^or$"), ignore.case = TRUE)) %>% 
      dplyr::rename_with(~ rename_value(.x, "BETA"), matches(c("^beta$"), ignore.case=TRUE)) %>% 
      dplyr::rename_with(~ rename_value(.x, "REF"), matches(c("^ref$"), ignore.case=TRUE)) %>% 
      dplyr::rename_with(~ rename_value(.x, "ALT"), matches(c("^alt$"), ignore.case=TRUE)) 
    if(locuszoomplot)
        df$color <- dfwithcolor$color
   
    if(! "P" %in% colnames(df)){
      df <- df %>% 
        dplyr::rename_with(~ rename_value(.x, "P"), matches(c("^pval$","^pvalue$","^p_value$","^p$"), ignore.case = TRUE)) 
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
    if( (!"CHROM" %in% colnames(df)) || (!"POS" %in% colnames(df)) || (!"P" %in% colnames(df) ))
      stop("Some columns are missing from the dataset. Required columns are CHROM,POS and P. Add the required columns and try again, or rename existing columns, e.g. df=df %>% dplyr::rename(CHROM=yourColname)")
     zero_pvals <- dplyr::filter(df, P <= .Machine$double.xmin )
    caption <- ""
    if (nrow(zero_pvals) > 0) {
      exclusion_list <- with(zero_pvals, paste(CHROM, POS, REF, ALT, sep = ":", collapse = ", "))
      exclusion_list <- stringr::str_trunc(exclusion_list, width = 100)
      caption <- stringr::str_glue("{nrow(zero_pvals)} zero-value p-value{if(nrow(zero_pvals) > 1) 's' else ''} found in the input dataset and removed from plot: {exclusion_list}")
      if(is.null(verbose) || verbose==TRUE ){ print(caption) }
    }
    df <- df %>% dplyr::filter(P > .Machine$double.xmin) 
    df$P <- as.numeric(df$P)
    df <- df[!is.na(df$P), ]
    dat[[i]] <- df
  }
  return(dat)
}
