
dat_check <- function(dat, verbose=TRUE,locuszoomplot=FALSE){
  if(is.data.frame(dat)) dat <- list(dat)
  for(i in seq_along(dat)){
    if(nrow(dat[[i]])==0){
      warning(paste0("Dataset [",i, "] is empty with 0 rows. Please remove it from the input list and re-run!"))
    }
    
  }
  dat <- dat_column_check_and_set(dat, verbose=verbose,locuszoomplot=locuszoomplot) %>% rm_chr_prefix_and_sort()
  return(dat)
}

downsample <- function(dat, downsample_cutoff=0.05, downsample_prop=0.1){
  for(i in seq_along(dat)){
    df <- as.data.frame(dat[[i]])
    if(max(df$P) > downsample_cutoff){
      sig <- df |>
      subset(P < downsample_cutoff)
      notsig <- df |>
      subset(P >= downsample_cutoff) |>
      group_by(CHROM) |>
      slice_sample(prop=downsample_prop)
      dat[[i]]  <- bind_rows(sig, notsig)
    }
  }
  return(dat)
}

rm_chr_prefix_and_sort <- function(dat){
  for(i in seq_along(dat)){
    df <- as.data.frame(dat[[i]])
    df <- df %>% dplyr::mutate(CHROM=gsub('chr','',CHROM))
    df$POS <- as.integer(df$POS)
    df <- df %>% dplyr::arrange(CHROM,POS)
    dat[[i]] <- df
  }
  return(dat)
}


convert_chrs_to_numeric <- function(dat, get_chr_lengths_from_data){
    chrs <- NULL
    if(get_chr_lengths_from_data)
      chrs <- get_chrs_from_data(dat)
    else{
      chrs <- chr_lengths$V1[2:24] 
      chrs <- gsub("chr", "", chrs)
    }
    numeric_chrs <- chrs[grepl('^-?[0-9.]+$', chrs)] %>% as.numeric() %>% sort()
    non_numeric_chrs <- chrs[!grepl('^-?[0-9.]+$', chrs)]
    chr_order <- append(c(numeric_chrs), c(non_numeric_chrs))
    chr_number=c(1:length(chr_order))
    chr_map <- stats::setNames(as.list(chr_number), chr_order)
   for(i in seq_along(dat)){
     df <- dat[[i]]
     for(chr in chr_order){
      df[df$CHROM==chr, 'CHROM'] <- as.numeric(unlist(unname(chr_map[chr])))
    }
    df$CRHOM <- as.integer(df$CHROM)
    df <- df %>% dplyr::arrange(CHROM) 
    dat[[i]] <- df
   }
  return(list("dat"=dat, "chr_map"=chr_map))
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
        dplyr::rename_with(~ rename_value(.x, "P"), matches(c("^pval$","^pvalue$","^p_value$","^p-value$","^p$"), ignore.case = TRUE)) 
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
      exclusion_list <- with(zero_pvals, paste(CHROM, POS, sep = ":", collapse = ", "))
      exclusion_list <- stringr::str_trunc(exclusion_list, width = 100)
      caption <- stringr::str_glue("{nrow(zero_pvals)} zero-value p-value{if(nrow(zero_pvals) > 1) 's' else ''} (P<",format_P(.Machine$double.xmin)," (.Machine$double.xmin)) found in the input dastaset and excluded from the plot: {exclusion_list}")
      if(is.null(verbose) || verbose==TRUE ){ print(caption) }
    }
    df <- df %>% dplyr::filter(P > .Machine$double.xmin) 
    df$P <- as.numeric(df$P)
    df <- df[!is.na(df$P), ]
    dat[[i]] <- df
  }
  return(dat)
}
