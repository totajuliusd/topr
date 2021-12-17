
get_genes_in_region <- function(chr=chr, xmin=xmin,xmax=xmax,protein_coding_only=F, show_exons=F,show_genes=T){
  if(show_genes){
    genes <- get_genes(chr,xmin,xmax,protein_coding_only = protein_coding_only )
  }else if(show_exons){
    genes <- get_exons(chr,xmin,xmax,protein_coding_only = protein_coding_only)
  }else{
    if(xmax-xmin < 1000001){
      genes <- get_exons(chr,xmin,xmax,protein_coding_only = protein_coding_only )
    }
    else{
      genes <- get_genes(chr,xmin,xmax,protein_coding_only = protein_coding_only )
    }
  }
  return(genes)
}

get_main_LD_snp <- function(dat){
  #df <- dat[[1]]
  label_cols <-  c("CHROM","POS","P","ID","log10p")
  top_snps <-  data.frame(matrix(nrow = 0, ncol = length(label_cols)))
  colnames(top_snps) <- label_cols
  for(i in seq_along(dat)){
      top_snps <- rbind(top_snps, dat[[i]] %>% dplyr::filter("R2" >= 1) %>% dplyr::distinct(ID, .keep_all=T) %>% dplyr::arrange(-R2) %>% utils::head(n=1) %>% dplyr::select("CHROM","POS","P","ID","log10p") )
  }
  return(top_snps)
}

get_annotation <- function(dat, annotate=1e-09, region_size=1000000,distinct_gene_labels=FALSE,protein_coding_only=FALSE, verbose=FALSE){
  if(is.data.frame(dat)){dat <- list(dat)}
  if("log10p" %in% colnames(dat[[1]])){
    label_cols <-  c("CHROM","POS","P","ID","Gene_Symbol","biotype", "log10p", "color","alpha","size","shape")
  }
  else{
    label_cols <-  c("CHROM","POS","P","ID","Gene_Symbol","biotype")
  }
  plot_labels <-  data.frame(matrix(nrow = 0, ncol = length(label_cols)))
  colnames(plot_labels) <- label_cols
  #retrieve the top variants

  for(i in seq_along(dat)){
    df <- as.data.frame(dat[[i]])
    if(is.vector(annotate)){
      annot_thresh <- ifelse(i <= length(annotate),  annotate[i], annotate[length(annotate)])
    }
    else{ annot_thresh <- annotate}
    tmp_labels <- get_best_snp_per_MB(df, thresh = annot_thresh, region_size=region_size, .checked=TRUE, protein_coding_only = protein_coding_only, verbose=FALSE)

    if(nrow(tmp_labels) > 0){
      if(!"biotype" %in% tmp_labels){tmp_labels$biotype <- "unknown"}
      if(! "Gene_Symbol" %in% colnames(tmp_labels)){
        tmp_labels <- annotate_with_nearest_gene(tmp_labels, protein_coding_only=protein_coding_only)
      }
      if("log10p" %in% colnames(dat[[1]])){
      tmp_labels <- tmp_labels %>% dplyr::select("CHROM","POS","P","ID","Gene_Symbol","biotype", "log10p", "color","alpha","size","shape")
      }
      else{
        tmp_labels <- tmp_labels %>% dplyr::select("CHROM","POS","P","ID","Gene_Symbol","biotype")
      }
    }
    plot_labels <- rbind(plot_labels,tmp_labels)
  }
  if(distinct_gene_labels){
    plot_labels <- plot_labels %>% dplyr::arrange(CHROM,POS,P) %>% dplyr::distinct(CHROM,POS,Gene_Symbol, .keep_all = TRUE)
  }
  return(plot_labels)
}


#' Get the nearest gene for one or more snps
#'
#' @description
#'
#' \code{annotate_with_nearest_gene()} Annotate the variant/snp with their nearest gene
#' Required parameters is a dataframe of SNPs (with the columns CHROM and POS)
#'
#' @param variants a dataframe of variant positions (CHROM and POS)
#' @param protein_coding_only Logical, if set to TRUE only annotate with protein coding genes (the default value is FALSE)
#' @return the input dataframe with Gene_Symbol as an additional column
#' @export
#'
#' @examples
#' \dontrun{
#' annotate_with_nearest_gene(variants)
#' }
annotate_with_nearest_gene <- function(variants, protein_coding_only=F){
  if("POS" %in% colnames(variants) & "CHROM" %in% colnames(variants)){
    if(length(variants$POS) > 1000){
      warning(paste("The dataset includes [",length(variants$POS),"] variants. This may take a while...", sep=""))
    }
    for(i in seq_along(variants$POS)){
      if(length(variants$POS) > 1000){
        if(i %% 1000==0) {
          # Print on the screen some message
          print(paste(i," variants annotated", sep=""))
        }
      }
      nearest_gene <- NULL
      variant <- variants[i,]
      chr <- gsub("chr", "", variant$CHROM)
      chr <- paste("chr",chr,sep="")
      genes_on_chr <- topr::ENSGENES %>% dplyr::filter(chrom == chr) %>% dplyr::arrange(gene_start)
      if(protein_coding_only){
        genes_on_chr <- genes_on_chr %>% dplyr::filter(biotype == "protein_coding")
      }
      within_gene <-  genes_on_chr %>% dplyr::filter(gene_end >= variant$POS & gene_start <= variant$POS)
      if(length(within_gene$gene_symbol) > 0 ){  #TODO: order the genes by their biotype, and pull out the top one
        if(length(within_gene) == 1){ nearest_gene <- within_gene$gene_symbol }
        else{
          prot_coding <- within_gene %>% dplyr::filter(biotype=="protein_coding")
          if(length(prot_coding$gene_symbol) > 0){  nearest_gene <- prot_coding %>% utils::head(n=1)}
          else{  nearest_gene <- within_gene %>% utils::head(n=1) }
        }
      }else{
        genes_left <-  genes_on_chr %>% dplyr::filter(gene_end <= variant$POS) %>% dplyr::arrange(gene_end)
        genes_right <- genes_on_chr %>% dplyr::filter(gene_start >= variant$POS) %>% dplyr::arrange(gene_start)
        if(length(genes_left$gene_symbol)>0  & length(genes_right$gene_symbol)> 0){
          gene_left <- genes_left[as.numeric(length(genes_left$gene_symbol)),]
          gene_right <- genes_right[1,]
          dist_left <- variant$POS-gene_left$gene_end
          dist_right <- gene_right$gene_start-variant$POS
          if(abs(dist_left) < abs(dist_right)){  nearest_gene <- gene_left }
          else{  nearest_gene <- gene_right }
        }
        else if(length(genes_left$gene_symbol)== 0  & length(genes_right$gene_symbol)> 0){
          nearest_gene <- genes_right[1,]
        }
        else if(length(genes_left$gene_symbol)> 0  & length(genes_right$gene_symbol)==0){
          nearest_gene <- genes_left[as.numeric(length(genes_left$gene_symbol)),]
        }

      }
      if(! is.null(nearest_gene)){
        variants[i,"Gene_Symbol"] <- nearest_gene$gene_symbol
        variants[i, "biotype"] <- nearest_gene$biotype
      }else{
        variants[i,"Gene_Symbol"] <- "not_found"
        variants[i, "biotype"] <- "."
      }
    }
  }
  else{
    stop("Cannot find the columns CHROM and POS in the input data. Add the required columns and try again, or rename existing columns, e.g. df=df %>% dplyr::rename(CHROM=yourColname)")
  }
  return(variants)
}