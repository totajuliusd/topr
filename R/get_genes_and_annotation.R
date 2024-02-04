#' Get SNPs/variants within region
#'
#' @description
#'
#' \code{get_genes_in_region()}
#'
#' @param region A string representing the genetic region (e.g chr16:50693587-50734041)
#' @param chr A string, chromosome (e.g. chr16)
#' @param xmin An integer representing genetic position 
#' @param xmax An integer representing genetic position
#' @param show_genes A logical scalar, show genes instead of exons (default show_genes=FALSE)
#' @param show_exons Deprecated : A logical scalar, show exons instead of genes (default show_exons=FALSE)
#' @inheritParams manhattan
#' @return the genes the requested region
#' @export
#'
#' @examples
#' \dontrun{
#' get_genes_in_region(region="chr16:50593587-50834041")
#'}


get_genes_in_region <- function(chr=chr, xmin=xmin,xmax=xmax,protein_coding_only=F, show_exons=F,show_genes=T, build=38, region=NULL){
  if (!missing(show_exons)) deprecated_argument_msg(show_exons)
   if(!is.null(region)){
    tmp <- unlist(stringr::str_split(region, ":"))
    chr <- tmp[1]
    tmp_pos <- unlist(stringr::str_split(tmp[2], "-"))
    xmin <- tmp_pos[1]
    xmax <- tmp_pos[2]
   }
  if(show_genes){
    genes <- get_genes(chr,xmin,xmax,protein_coding_only = protein_coding_only, build=build)
  }else{
    genes <- get_exons(chr,xmin,xmax,protein_coding_only = protein_coding_only, build=build)
  }
  return(genes)
}

get_main_LD_snp <- function(dat, nudge_x=0.1,nudge_y=0.1,angle=0,label_fontface="plain",label_family="",label_alpha=1){
  label_cols <-  c("CHROM","POS","P","ID","log10p", "alpha","shape","nudge_x","nudge_y","angle","fontface","family")
  top_snps <-  data.frame(matrix(nrow = 0, ncol = length(label_cols)))
  colnames(top_snps) <- label_cols
  for(i in seq_along(dat)){
    top_snps <- rbind(top_snps, dat[[i]] %>% dplyr::filter(R2 >= 1) %>% dplyr::distinct(ID, .keep_all=T) %>% dplyr::arrange(-R2) %>% utils::head(n=1) %>% dplyr::select("CHROM","POS","P","ID","log10p") )
    if(length(top_snps$P) == 0){
      print("Could not find the lead snps (with R2==1) in the dataset!")      
    }
  }
  top_snps$nudge_x <- nudge_x; top_snps$nudge_y=nudge_y; top_snps$fontface<- label_fontface; top_snps$family <- label_family; top_snps$angle <- angle; top_snps$alpha=label_alpha;
  return(top_snps)
}
get_annotation <- function(dat, annotate=5e-08, region_size=1000000,distinct_gene_labels=FALSE,protein_coding_only=FALSE, verbose=NULL,nudge_x=0.1,nudge_y=0.1,
                           angle=0,label_fontface="plain",label_family="",build=38, label_alpha=1, chr_map=NULL){
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
    if(is.vector(region_size)){
      region_size_tmp <- ifelse(i <= length(region_size),  region_size[i], region_size[length(region_size)])
    }else{region_size_tmp=region_size}
    if(is.vector(annotate)){
      annot_thresh <- ifelse(i <= length(annotate),  annotate[i], annotate[length(annotate)])
    } else{ annot_thresh <- annotate}
    if(is.vector(nudge_x)){
      nudge_x_tmp <- ifelse(i <= length(nudge_x),  nudge_x[i], nudge_x[length(nudge_x)])
    }else{nudge_x_tmp <- nudge_x}
    if(is.vector(nudge_y)){
      nudge_y_tmp <- ifelse(i <= length(nudge_y),  nudge_y[i], nudge_y[length(nudge_y)])
    }else{nudge_y_tmp <- nudge_y}
    
    if(is.vector(angle)){
      angle_tmp <- ifelse(i <= length(angle),  angle[i], angle[length(angle)])
    }else{angle_tmp <- angle}
    
    if(is.vector(label_fontface)){
      label_fontface_tmp <- ifelse(i <= length(label_fontface),  label_fontface[i], label_fontface[length(label_fontface)])
    }else{label_fontface_tmp <- label_fontface}
    
    if(is.vector(label_alpha)){
      alpha_tmp <- ifelse(i <= length(label_alpha),  label_alpha[i], label_alpha[length(label_alpha)])
    }else{alpha_tmp <- label_alpha}
    
    
    if(is.vector(label_family)){
      label_family_tmp <- ifelse(i <= length(label_family),  label_family[i], label_family[length(label_family)])
    }else{label_family_tmp <- label_family}
    if(! is.null(verbose)){
      if(verbose){
        print(paste0("Dataset: [",i,"]:"))
      }
    }
    tmp_labels <- get_lead_snps(df, thresh = annot_thresh, region_size=region_size_tmp, .checked=TRUE, protein_coding_only = protein_coding_only, verbose=verbose, keep_chr=FALSE)
    if(nrow(tmp_labels) > 0){
      tmp_labels$nudge_x <- nudge_x_tmp
      tmp_labels$nudge_y <- nudge_y_tmp
      tmp_labels$angle <- angle_tmp
      tmp_labels$fontface <- label_fontface_tmp
      tmp_labels$family <- label_family_tmp
      tmp_labels$alpha <- alpha_tmp
      if(!"biotype" %in% tmp_labels){tmp_labels$biotype <- "unknown"}
      if(! "Gene_Symbol" %in% colnames(tmp_labels)){
        tmp_labels <- annotate_with_nearest_gene(tmp_labels, protein_coding_only=protein_coding_only, build=build, .chr_map=chr_map)
      }
      if("log10p" %in% colnames(dat[[1]])){
        tmp_labels <- tmp_labels %>% dplyr::select("CHROM","POS","P","ID","Gene_Symbol","biotype", "log10p", "color","alpha","size","shape","nudge_x","nudge_y","angle","fontface","family")
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
#' @param build A number representing the genome build. Set to 37 to change to build (GRCh37). The default is build 38 (GRCh38).
#' @param .chr_map An internally used list which maps chromosome names to numbers.
#' @return the input dataframe with Gene_Symbol as an additional column
#' @export
#'
#' @examples
#' \dontrun{
#' variants <-get_lead_snps(CD_UKBB)
#' annotate_with_nearest_gene(variants)
#' }
#' 
annotate_with_nearest_gene <- function(variants, protein_coding_only=FALSE, build=38, .chr_map=NULL){
  if("POS" %in% colnames(variants) & "CHROM" %in% colnames(variants)){
    if(length(variants$POS) > 1000){
      print(paste("The dataset includes [",length(variants$POS),"] variants. This may take a while...", sep=""))
    }
    for(i in seq_along(variants$POS)){
      if(length(variants$POS) > 1000){
        if(i %% 1000==0) {
          print(paste(i," variants annotated", sep=""))
        }
      }
      nearest_gene <- NULL
      variant <- variants[i,]
      chr <- gsub("chr", "", variant$CHROM)
      if(! is.null(.chr_map))
         chr <- names(.chr_map[as.numeric(chr)])
      genes_on_chr <- get_db(build) %>% dplyr::filter(chrom == chr) %>% dplyr::arrange(gene_start)
      
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
