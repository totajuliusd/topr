
set_size_shape_alpha <- function(df,size,shape,alpha){
  df <- set_size(df,size)
  df <- set_shape(df,shape)
  df <- set_alpha(df,alpha)
  df <- set_shape_by_impact(df)
  return(df)
}

set_shape_by_impact <- function(dat){
  for(i in seq_along(dat)){
    if("Max_Impact" %in% colnames(dat[[i]])){
      dat[[i]]$shape0 <- ifelse(dat[[i]]$Max_Impact == "MODERATE", 17, dat[[i]]$shape)
      dat[[i]]$shape <- ifelse(dat[[i]]$Max_Impact =="HIGH", 8, dat[[i]]$shape0)
    }
  }
  return(dat)
}

set_annotate <- function(dat,annotate){
  if(! is.null(annotate)){
    for(i in seq_along(dat)){
      if(!is.null(annotate) & nrow(dat[[i]]>0)){
        if(is.vector(annotate) & (length(annotate) >= i))
          dat[[i]]$annotate <- annotate[i]
        else
          dat[[i]]$annotate <- annotate
      }
    }
  }
  return(dat)
}

set_color <- function(dat,color){
  if(length(color) < length(dat))
    stop(paste("There are ",length(dat), " datasets, but only ",length(color), " color. Add more colors with the color argument, eg. color=c(\"blue\",\"yellow\",\"green\", etc)"))
  for(i in seq_along(dat)){
    if(!is.null(color) & nrow(dat[[i]]>0)){
      if(is.vector(color) & (length(color) >= i))
        dat[[i]]$color <- color[i]
      else
        dat[[i]]$color <- color
    }
  }
  return(dat)
}

set_size <- function(dat,size){
  if(is.null(size))
    size <- 1
  for(i in seq_along(dat)){
    if(!is.null(size) & nrow(dat[[i]]>0)){
      if(is.vector(size) &  length(size) >= i )
        dat[[i]]$size <- size[i]
      else
        dat[[i]]$size <- size
    }
  }
  return(dat)
}

set_alpha <- function(dat,alpha){
  if(is.null(alpha))
    alpha <- 0.7
  for(i in seq_along(dat)){
    if(!is.null(alpha) & nrow(dat[[i]]>0)){
      if(is.vector(alpha) & (length(alpha) >= i))
        dat[[i]]$alpha <- alpha[i]
      else
        dat[[i]]$alpha <- alpha
    }
  }
  return(dat)
}

set_shape <- function(dat,shape){
  if(is.null(shape))
    shape <- 19
  for(i in seq_along(dat)){
    if(!is.null(shape) & nrow(dat[[i]]>0)){
      if(is.vector(shape) & (length(shape) >= i))
        dat[[i]]$shape <- shape[i]
      else
        dat[[i]]$shape <- shape
    }
  }
  return(dat)
}

set_log10p <- function(dat, ntop){
  for(i in seq_along(dat)){
    df <- dat[[i]]
    if(ntop > length(dat)){ ntop <- length(dat)}
    if(i <= ntop){
      dat[[i]] <- df %>% dplyr::mutate(log10p=-log10(P))

    }else{
      dat[[i]] <- df %>% dplyr::mutate(log10p=log10(P))
    }
  }
  return(dat)
}

## Getters

get_genes <- function(chr,xmin=0,xmax=NULL,protein_coding_only=FALSE){
  chr <- gsub("chr", "", chr)
  chr <- paste("chr",chr,sep="")
  if(is.null(xmax)){
    xmax <- chr_lengths[chr_lengths$V1 == chr,]$V2
  }
  if(protein_coding_only)
    genes <- topr::ENSGENES %>% dplyr::filter(chrom==chr) %>% dplyr::filter(biotype=="protein_coding") %>% dplyr::filter(gene_start>xmin & gene_start<xmax | gene_end > xmin & gene_end < xmax | gene_start < xmin & gene_end>xmax)
  else
    genes <- topr::ENSGENES %>% dplyr::filter(chrom==chr) %>% dplyr::filter(gene_start>xmin & gene_start<xmax | gene_end > xmin & gene_end < xmax | gene_start < xmin & gene_end>xmax)
  return(genes)
}

get_exons <- function(chr,xmin=0,xmax=NULL,protein_coding_only=FALSE){
  genes <- get_genes(chr,xmin,xmax)
  if(protein_coding_only)
    exons <- merge(topr::ENSEXONS, genes %>% dplyr::select(gene_symbol,biotype,chrom,gene_start,gene_end) %>% dplyr::filter(biotype=="protein_coding"), by=c("gene_symbol","chrom","gene_start","gene_end"))
  else
    exons <- merge(topr::ENSEXONS, genes %>% dplyr::select(gene_symbol,biotype,chrom,gene_start,gene_end), by=c("gene_symbol","chrom","gene_start","gene_end"))
    #rename chromstart exon_chromstart | rename chromend exon_chromend)",sep=""))
  return(exons)
}


get_chr_from_df <- function(df){
  #take the chromosome from the dataframe
  if(length(unique(df$CHROM)) == 1)
    chr <- unique(df$CHROM)
  else
    stop("More than one chromosome in the input data, please select which chromosome to display by setting the chr argument (e.g. chr='chr10') or use input data with just one chromosome")
  return(chr)
}


get_shades <- function(offsets,dat,ntop=ntop,include_chrX=FALSE,ymin=NULL,ymax=NULL){
  n_offsets <- 11
  if(!(include_chrX)) n_offsets <-10
  y1 <- c(rep(0, n_offsets))  #if there is no bottom plot
  if(is.null(ymin)){
    ymin <- get_ymin(dat)
  }
  if(length(dat) > ntop){
    y1 <- c(rep(ymin, n_offsets))
  }
  if(is.null(ymax)){
    ymax <- get_ymax(dat)
  }
  ##### here we need to accoutn for whether we have chr X or not
  if(include_chrX){
    shades <- data.frame(x1=c(offsets[[2]],offsets[[4]],offsets[[6]],offsets[[8]],offsets[[10]],offsets[[12]],offsets[[14]],offsets[[16]],offsets[[18]],offsets[[20]],offsets[[22]]),
                      x2=c(offsets[[3]],offsets[[5]],offsets[[7]], offsets[[9]],offsets[[11]],offsets[[13]],offsets[[15]],offsets[[17]],offsets[[19]],offsets[[21]],offsets[[23]]),
                      y1=y1,
                      y2=c(rep(ymax, n_offsets)))
  }
  else{
    shades <- data.frame(x1=c(offsets[[2]],offsets[[4]],offsets[[6]],offsets[[8]],offsets[[10]],offsets[[12]],offsets[[14]],offsets[[16]],offsets[[18]],offsets[[20]]),
                      x2=c(offsets[[3]],offsets[[5]],offsets[[7]], offsets[[9]],offsets[[11]],offsets[[13]],offsets[[15]],offsets[[17]],offsets[[19]],offsets[[21]]),
                      y1=y1,
                      y2=c(rep(ymax, n_offsets)))
  }
  return(shades)
}
get_ymax <- function(dat){
  ymax <- 0
  if(is.data.frame(dat)) dat <- list(dat)
  for(i in seq_along(dat)){
    tmp <- max(dat[[i]]$log10p)
    if(tmp>ymax){
      ymax <- tmp
    }
  }
  return(ymax)
}

get_ymin <- function(dat){
  ymin <- 5
  if(is.data.frame(dat)) dat <- list(dat)
    for(i in seq_along(dat)){
      tmp <- min(dat[[i]]$log10p)
      if(tmp<ymin){
        ymin <- tmp
      }
    }
  return(ymin)
}

get_chr_offsets <- function(include_chrX=F){
  #create the offsets from the internal chr_lengths data
  chr_lengths$CHROM <- gsub("chr","", chr_lengths$V1)
  chr_lengths <- chr_lengths %>% dplyr::filter(! V1 %in% c("chrM") )
  chr_lengths[chr_lengths$CHROM=="X",'CHROM'] <- "23"
  chr_lengths[chr_lengths$CHROM=="Y",'CHROM'] <- "24"
  chr_lengths$CHROM <- as.integer(chr_lengths$CHROM)
  no_chrs <- 22
  if(include_chrX)
    no_chrs <- 23
  chr_lengths <- chr_lengths %>% dplyr::filter(CHROM< no_chrs +1 )
  tmp <- chr_lengths %>% dplyr::group_by(CHROM) %>% dplyr::summarize(m=V2) %>% dplyr::mutate(offset=cumsum(lag(m, default=0)))
  offsets <- stats::setNames(tmp$offset,tmp$CHROM)
  return(offsets)
}

get_ticknames <- function(df){
  no_chrs <- ifelse("chrX" %in% df$CHROM || "X" %in% df$CHROM || "chr23" %in% df$CHROM || "23" %in% df$CHROM, 23, 22)
  if(no_chrs == 23){
    ticknames <- c(1:16, '',18, '',20, '',22, 'X')
  }else{
    ticknames <- c(1:16, '',18, '',20, '',22)
  }
  tickpos <-df %>% dplyr::group_by(CHROM)%>%
                                dplyr::summarize(pm=mean(POS))%>%
                                dplyr::pull(pm)

  names(tickpos) <- NULL
  return(list(names=ticknames, pos=tickpos))
}

get_ticks <- function(dat){
  df <- dat[[1]]
  for(i in seq_along(dat)){ if(length(unique(dat[[i]]$CHROM))  > length(unique(df$CHROM))){ df <- dat[[i]] } }
  no_chrs <- ifelse("chrX" %in% df$CHROM || "X" %in% df$CHROM || "chr23" %in% df$CHROM || "23" %in% df$CHROM, 23, 22)
  if(no_chrs == 23){
    ticknames <- c(1:16, '',18, '',20, '',22, 'X')
  }else{
    ticknames <- c(1:16, '',18, '',20, '',22)
  }
  tickpos <- df %>%
                                dplyr::group_by(CHROM)%>%
                                dplyr::summarize(pm=mean(POS))%>%
                                dplyr::pull(pm)

  names(tickpos) <- NULL
  return(list(names=ticknames, pos=tickpos))
}

get_pos_with_offset4list <- function(dat, offsets){
  for(i in seq_along(dat)){
    dat[[i]] <- get_pos_with_offset(dat[[i]], offsets)
  }
  return(dat)
}

get_pos_with_offset <- function(df,offsets){
  df$POS_orig <- df$POS
  df <- df %>% dplyr::mutate(POS=POS_orig+offsets[CHROM])
  return(df)
}


#' Get the top SNP per 10 MB
#'
#' @description
#'
#' \code{get_best_snp_per_MB()} Get the top SNP per 10 MB
#' All other input parameters are optional
#'
#' @param df Dataframe
#' @param thresh P-value threshold, only consider variants with p-values below this threshold (1e-09 by default)
#' @param region Get the top/best variant (with p-value below thresh) within this region (region=1000000 by default)
#' @param protein_coding_only Set this variable to TRUE to only use protein_coding genes for annotation
#' @param chr Use this argument to get the top variants from one chromosome only
#' @param .checked If the input data has already been checked, this can be set to TRUE so it wont be checked again (FALSE by default)
#' @return Dataframe of lead variants. Returns the best variant per MB (by default, change the region size with the region argument) with p-values below the input threshold (thresh=1e-09 by default)
#' @export
#'
#' @examples
#' \dontrun{
#' data(gwas_CD)
#' get_best_snp_per_MB(gwas_CD, thresh = 1e-09, region = 10000000)
#' }
get_best_snp_per_MB <- function(df, thresh=1e-09,region=1000000,protein_coding_only=FALSE,chr=NULL, .checked=FALSE){
  variants <- df[0,]
  dat <- df
  if(is.data.frame(dat)) dat <- list(dat)
  if(! .checked){
    dat <- dat_check(dat)
  }
  df <- dat[[1]]
 # df$CHROM=gsub("chr","", df$CHROM)
  if(! is.null(chr)){
    chr <- gsub("chr","", chr)
    df <- df %>% dplyr::filter(CHROM==chr)
  }
  df <- df %>% dplyr::filter(P<thresh)
  if(protein_coding_only & ("biotype" %in% colnames(df))){
       df <- df %>% dplyr::filter(biotype == "protein_coding")
  }
  if(nrow(df) > 0){
    df$tmp <- NA
    for(row in seq_len(nrow(df))){
      df$tmp <- base::round(df$POS/region)
    }
    lead_snps <- df %>% dplyr::group_by(CHROM,tmp) %>% dplyr::arrange(P) %>% dplyr::filter(P== min(P))
    if(length(lead_snps)== 0){
        print(paste("There are no SNPs with a p-value below ",thresh, " in the input. Use the [thresh] argument to adjust the threshold.",sep=""))
    }
    variants <- dplyr::ungroup(lead_snps)%>% dplyr::distinct(CHROM,POS, .keep_all = T)
  }
  return(variants)
}

#' Get the genetic position of a gene or genes by their gene name
#'
#' @description
#'
#' \code{get_genes_by_Gene_Symbol()} Get genes by their gene symbol/name
#' Required parameters is on gene name or a vector of gene names
#'
#' @param genes A gene name (e.g. "FTO") or a vector of gene names ( c("FTO","NOD2"))
#' @param chr Search for the genes on this chromsome only
#' @return Dataframe of genes
#' @export
#'
#' @examples
#' \dontrun{
#' get_genes_by_Gene_Symbol(c("FTO","THADA"))
#' }
#'
get_genes_by_Gene_Symbol <- function(genes, chr=NULL){
    genes_df <- topr::ENSGENES %>% dplyr::filter(gene_symbol %in% genes)
    if(! is.null(chr)){
      chr <- gsub('chr','',chr)
      genes_df <- genes_df %>% dplyr::filter(chrom == paste("chr",chr,sep=""))
  }
  genes_df <- genes_df %>% dplyr::mutate(POS = gene_start + round((gene_end-gene_start)/2)) %>% dplyr::rename(Gene_Symbol=gene_symbol,CHROM=chrom)
  dist_genes <- genes_df %>% dplyr::distinct(CHROM,gene_start,gene_end, .keep_all = T)
  return(dist_genes)
}

#' Get the genetic position of a gene or genes by their gene name
#'
#' @description
#'
#' \code{get_gene()} Get gene coordinates by its gene symbol/name
#' Required parameters is on gene name or a vector of gene names
#'
#' @param gene_name string A gene name (e.g. "FTO")
#' @param chr Search for the gene on this chromosome
#' @return Dataframe of genes
#' @export
#'
#' @examples
#' \dontrun{
#' get_gene("FTO")
#' }
#'

get_gene  <- function(gene_name,chr=NULL){
  if(! is.null(chr)){
    chr <- gsub("chr", "", chr)
    chr <- paste("chr",chr,sep="")
    gene <- topr::ENSGENES %>% dplyr::filter(chrom == chr) %>% dplyr::filter(gene_symbol == gene_name)
  }
  else{
    gene <- topr::ENSGENES %>% dplyr::filter(gene_symbol == gene_name)
  }
  if(is.null(gene)){
    print(paste("Could not find a gene with the gene name: [", gene_name, "]", sep=""))
  }
  return(gene)
}

get_legend<-function(p1){
  tmp <- ggplot_gtable(ggplot_build(p1))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



#' Get the top hit from the dataframe
#'
#' @description
#'
#' \code{get_top_hit()} Get the top hit from the dataframe
#' All other input parameters are optional
#'
#' @param df Dataframe containing association results
#' @param chr Get the top hit in the data frame for this chromosome. If chromosome is not provided, the top hit from the entire dataset is returned.

#' @return Dataframe containing the top hit
#' @export
#'
#' @examples
#' \dontrun{
#' data(gwas_CD)
#' get_top_hit(gwas_CD, chr="chr1")
#' }
#'
get_top_hit <- function(df, chr=NULL){
  if (!is.null(chr)) {
    chr_tmp <- gsub('chr','',chr)
    chr <- paste("chr",chr_tmp, sep="")
    df <- df %>%
      filter(CHROM == chr)
  } else {
    print("Returning the top hit over the entire genome, since chromosome was not provided as an argument")
  }
  top_hit <- df %>%
    dplyr::arrange(P) %>%
    head(n = 1)

  return(top_hit)
}



#' Get the top hit from the dataframe
#'
#' @description
#'
#' \code{get_topr_colors()} Get the top hit from the dataframe
#' All other input parameters are optional
#'
#' @return Vector of colors used for plotting
#' @export
#'
#' @examples
#' \dontrun{
#' get_topr_colors()
#' }
#'
get_topr_colors <- function(){
  return(c("darkblue","#E69F00","#00AFBB","#999999","#FC4E07","darkorange1"))
}

get_manhattan_options <- function(){
  print("options(repr.plot.width=12, repr.plot.height=4.5,repr.plot.res = 170)")
}

get_regionplot_options <- function(){
  print("options(repr.plot.width=11, repr.plot.height=6.5,repr.plot.res = 170)")
}


set_genes_pos_adj <- function(genes, offsets){
  #CHROM,POS,Gene_Symbol
  genes <- genes %>% dplyr::mutate(CHROM=gsub('chr','',CHROM))
  genes[genes$CHROM=='X', 'CHROM'] <- "23"
  genes$CHROM <- as.integer(genes$CHROM)
  genes <- genes %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
  return(genes)
}