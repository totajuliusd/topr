
set_size_shape_alpha <- function(df,size,shape,alpha,locuszoomplot=F){
  df <- set_size(df,size, locuszoomplot = locuszoomplot)
  df <- set_shape(df,shape,locuszoomplot = locuszoomplot)
  df <- set_alpha(df,alpha)
  if(!locuszoomplot){
    df <- set_shape_by_impact(df)
  }
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

set_size <- function(dat,size, locuszoomplot = locuszoomplot){
  if(is.null(size))
    size <- 1
  for(i in seq_along(dat)){
    if(!is.null(size) & nrow(dat[[i]]>0)){
      if(is.vector(size) &  length(size) >= i )
        dat[[i]]$size <- size[i]
      else
        dat[[i]]$size <- size
    }
    if(locuszoomplot & ("R2" %in% colnames(dat[[i]]))){
      dat[[i]]$size <- ifelse(dat[[i]]$R2 == 1, dat[[i]]$size*1.5, dat[[i]]$size)
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

set_shape <- function(dat,shape,locuszoomplot=F){
  if(is.null(shape))
    shape <- 19
  for(i in seq_along(dat)){
    if(!is.null(shape) & nrow(dat[[i]]>0)){
      if(is.vector(shape) & (length(shape) >= i))
        dat[[i]]$shape <- shape[i]
      else
        dat[[i]]$shape <- shape
    }
    if(locuszoomplot & ("R2" %in% colnames(dat[[i]]))){
        dat[[i]]$shape <- ifelse(dat[[i]]$R2 == 1, 18, dat[[i]]$shape)
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

get_genes <- function(chr,xmin=0,xmax=NULL,protein_coding_only=FALSE, build=38){
  chr <- gsub("chr", "", chr)
  chr <- paste("chr",chr,sep="")
  if(is.null(xmax)){
    xmax <- chr_lengths[chr_lengths$V1 == chr,]$V2
  }
  if(chr == "chr23"){ chr="chrX"}
  if(chr == "23"){ chr="X"}
  if(build == "38")
    genes <- toprdata::ENSGENES %>% dplyr::filter(chrom==chr)  %>% dplyr::filter(gene_start>xmin & gene_start<xmax | gene_end > xmin & gene_end < xmax | gene_start < xmin & gene_end>xmax | gene_start == xmin | gene_end == xmax)
  else if(build == "37")
    genes <- toprdata::ENSGENES_37 %>% dplyr::filter(chrom==chr)  %>% dplyr::filter(gene_start>xmin & gene_start<xmax | gene_end > xmin & gene_end < xmax | gene_start < xmin & gene_end>xmax | gene_start == xmin | gene_end == xmax)
  else
    print("The requested build does not exist. Available genome builds are: 37 and 38")
  if(protein_coding_only)  
      genes <- genes %>% dplyr::filter(biotype=="protein_coding")
  return(genes)
}

get_exons <- function(chr,xmin=0,xmax=NULL,protein_coding_only=FALSE, build=38){
  genes <- get_genes(chr,xmin,xmax, build=build)

  if(protein_coding_only){
    if(build == "38")
      exons <- merge(toprdata::ENSEXONS, genes %>% dplyr::select(gene_symbol,biotype,chrom,gene_start,gene_end) %>% dplyr::filter(biotype=="protein_coding"), by=c("gene_symbol","chrom","gene_start","gene_end"))
    else if(build == "37")
      exons <- merge(toprdata::ENSEXONS_37, genes %>% dplyr::select(gene_symbol,biotype,chrom,gene_start,gene_end) %>% dplyr::filter(biotype=="protein_coding"), by=c("gene_symbol","chrom","gene_start","gene_end"))
  }
  else{
    if(build == "38")
      exons <- merge(toprdata::ENSEXONS, genes %>% dplyr::select(gene_symbol,biotype,chrom,gene_start,gene_end), by=c("gene_symbol","chrom","gene_start","gene_end"))
    else if(build == "37")
      exons <- merge(toprdata::ENSEXONS_37, genes %>% dplyr::select(gene_symbol,biotype,chrom,gene_start,gene_end), by=c("gene_symbol","chrom","gene_start","gene_end"))
  }
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
  y1 <- c(rep(0, n_offsets))  #if there is no bottom plot
  if(is.null(ymin))
    ymin <- get_ymin(dat)
  if(length(dat) > ntop)
    y1 <- c(rep(ymin, n_offsets))
  if(is.null(ymax)){
    ymax <- get_ymax(dat)
  }
  offset23 <- offsets[[23]]
   if(! include_chrX)
    offset23 <- offset23+4000000
  shades <- data.frame(x1=c(offsets[[2]],offsets[[4]],offsets[[6]],offsets[[8]],offsets[[10]],offsets[[12]],offsets[[14]],offsets[[16]],offsets[[18]],offsets[[20]],offsets[[22]]),
                      x2=c(offsets[[3]],offsets[[5]],offsets[[7]], offsets[[9]],offsets[[11]],offsets[[13]],offsets[[15]],offsets[[17]],offsets[[19]],offsets[[21]],offset23),
                      y1=y1,
                      y2=c(rep(ymax, n_offsets)))
  return(shades)
}


get_ymax <- function(dat){
  ymax <- 0
  if(is.data.frame(dat)) dat <- list(dat)
  for(i in seq_along(dat)){
    if(length(dat[[i]]$log10p) > 0){
      tmp <- max(dat[[i]]$log10p)
      if(tmp>ymax){
        ymax <- tmp
      }
    }
  }
  return(ymax)
}

get_ymin <- function(dat){
  ymin <- 5
  if(is.data.frame(dat)) dat <- list(dat)
    for(i in seq_along(dat)){
      if(length(dat[[i]]$log10p) > 0){
        tmp <- min(dat[[i]]$log10p)
        if(tmp<ymin){
          ymin <- tmp
      }
    }
    }
  return(ymin)
}

get_chr_lengths_and_offsets <- function(include_chrX=F){
  #create the offsets from the internal chr_lengths data
  chr_lengths$CHROM <- gsub("chr","", chr_lengths$V1)
  chr_lengths <- chr_lengths %>% dplyr::filter(! V1 %in% c("chrM") )
  chr_lengths[chr_lengths$CHROM=="X",'CHROM'] <- "23"
  chr_lengths[chr_lengths$CHROM=="Y",'CHROM'] <- "24"
  chr_lengths$CHROM <- as.integer(chr_lengths$CHROM)
  no_chrs <- 23
  #no_chrs <-ifelse(include_chrX, 23,22)
  chr_lengths <- chr_lengths %>% dplyr::filter(CHROM< no_chrs +1 )
  chr_lengths_and_offsets <- chr_lengths %>% dplyr::group_by(CHROM) %>% dplyr::summarize(m=V2) %>% dplyr::mutate(offset=cumsum(as.numeric(lag(m, default=0))))
  return(chr_lengths_and_offsets)
}

get_chr_offsets <- function(include_chrX=F){
  tmp <- get_chr_lengths_and_offsets(include_chrX)
  offsets <- stats::setNames(tmp$offset,tmp$CHROM)
  return(offsets)
}

get_ticknames <- function(df){
  no_chrs <- ifelse("chrX" %in% df$CHROM || "X" %in% df$CHROM || "chr23" %in% df$CHROM || "23" %in% df$CHROM, 23, 22)
  include_chrX <- T
  if(no_chrs == 23){
    ticknames <- c(1:16, '',18, '',20, '',22, 'X')
  }else{
    ticknames <- c(1:16, '',18, '',20, '',22)
    include_chrX <- F
  }
  chr_lengts_and_offsets <- get_chr_lengths_and_offsets(include_chrX)
  tickpos <-chr_lengts_and_offsets %>% dplyr::group_by(CHROM)%>%
                  dplyr::summarize(pm=offset+(m/2))%>%
    #    dplyr::summarize(pm=mean(POS))
                                dplyr::pull(pm)

  names(tickpos) <- NULL
  return(list(names=ticknames, pos=tickpos))
}

get_ticks <- function(dat){
  df <- dat[[1]]
  for(i in seq_along(dat)){ if(length(unique(dat[[i]]$CHROM))  > length(unique(df$CHROM))){ df <- dat[[i]] } }
  ticknames <- c(1:16, '',18, '',20, '',22, 'X')
  chr_lengts_and_offsets <- get_chr_lengths_and_offsets()
  tickpos <-chr_lengts_and_offsets %>% dplyr::group_by(CHROM)%>%
    dplyr::summarize(pm=offset+(m/2))%>%
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


#' Get the index/lead variants
#'
#' @description
#'
#' \code{get_lead_snps()} Get the top variants within 1 MB windows of the genome with association p-values below the given threshold 
#'
#'
#' @param df Dataframe
#' @param thresh A number. P-value threshold, only extract variants with p-values below this threshold (5e-08 by default)
#' @param protein_coding_only Logical, set this variable to TRUE to only use protein_coding genes for annotation
#' @param chr String, get the top variants from one chromosome only, e.g. chr="chr1"
#' @param .checked Logical, if the input data has already been checked, this can be set to TRUE so it wont be checked again (FALSE by default)
#' @param verbose Logical, set to TRUE to get printed information on number of SNPs extracted
#' @param keep_chr Logical, set to FALSE to remove the "chr" prefix before each chromosome if present (TRUE by default) 
#' @return Dataframe of lead variants. Returns the best variant per MB (by default, change the region size with the region argument) with p-values below the input threshold (thresh=5e-08 by default)
#' @export
#' @inheritParams regionplot
#' @examples
#' \dontrun{
#' get_lead_snps(CD_UKBB)
#' }
#'
get_lead_snps <- function(df, thresh=5e-08,region_size=1000000,protein_coding_only=FALSE,chr=NULL, .checked=FALSE, verbose=NULL, keep_chr=TRUE){
  variants <- df[0,]
  dat <- df
  if(is.data.frame(dat)) dat <- list(dat)
  chrPrefix=0
  
  dat[[1]] <- dat_column_chr_pos_check(dat[[1]]) 
  
  if("CHROM"  %in% colnames(dat[[1]])){
    if(!is.integer(dat[[1]][1,"CHROM"])){chrPrefix=1}
  }
  if(! is.numeric(region_size)){
    region_size <- convert_region_size(region_size)
  }
  if(! .checked){
    dat <- dat_check(dat, verbose=verbose)
  }
 
  df <- dat[[1]]
  if(! is.null(chr)){
    chr <- gsub("chr","", chr)
    df <- df %>% dplyr::filter(CHROM == chr)
  }
  
  df <- df %>% dplyr::filter(P<thresh)
  if(protein_coding_only & ("biotype" %in% colnames(df))){
    df <- df %>% dplyr::filter(biotype == "protein_coding")
  }
  if(nrow(df) > 0){
    df$tmp <- NA
    if(max(df$POS)-min(df$POS) <= region_size){
      df$tmp<- 1
    }else{
      for(row in seq_len(nrow(df))){
        df$tmp <- base::round(df$POS/region_size)
      }
    }
    lead_snps <- df %>% dplyr::group_by(CHROM,tmp) %>% dplyr::arrange(P) %>% distinct(P, .keep_all=T) %>% dplyr::filter(P == min(P))
    variants <- dplyr::ungroup(lead_snps)%>% dplyr::distinct(CHROM,POS, .keep_all = T)
    if(! is.null(verbose)){
      if(verbose){
        print(paste("Found ",length(variants$POS), " index/lead variants with a p-value below ", thresh, sep=""))
        print(paste("Index/lead variant: the top variant within a",format(region_size, scientific = F),"bp region_size"), sep="")
      }
    }
  }else{
    if(is.null(verbose)){
      print(paste("There are no SNPs with p-values below ",thresh, " in the input dataset. Use the [thresh] argument to lower the threshold.",sep=""))
    }else if(verbose){
      print(paste("There are no SNPs with p-values below ",thresh, " in the input dataset. Use the [thresh] argument to lower the threshold.",sep=""))
    }
    
  }
  if("tmp" %in% colnames(variants)){
    variants <- variants %>% dplyr::select(-tmp)
  }
  if(keep_chr & chrPrefix){
    variants$CHROM <- gsub("(chr|CHR|Chr)","", variants$CHROM)
    variants$CHROM  <- paste0("chr", variants$CHROM)
  }
  return(variants)
}


#' Get the genetic position of a gene by its gene name
#'
#' @description
#'
#' \code{get_genes_by_Gene_Symbol()} Get genes by their gene symbol/name
#' Required parameters is on gene name or a vector of gene names
#'
#' @param genes A string or vector of strings representing gene names, (e.g. "FTO") or  (c("FTO","NOD2"))
#' @param chr A string, search for the genes on this chromosome only, (e.g chr="chr1")
#' @param build A string, genome build, choose between builds 37 (GRCh37) and 38 (GRCh38) (default is 38)
#' @return Dataframe of genes
#' @export
#'
#' @examples
#' \dontrun{
#'   get_genes_by_Gene_Symbol(c("FTO","THADA"))
#' }
#'
get_genes_by_Gene_Symbol <- function(genes, chr=NULL, build=38){
  if(build == "38"){
    genes_df <- toprdata::ENSGENES %>% dplyr::filter(gene_symbol %in% genes)
  }else if(build == "37"){ 
    genes_df <- toprdata::ENSGENES_37 %>% dplyr::filter(gene_symbol %in% genes)
  }else{
    warning(paste("Build [",build,"] not found!!!!!!  Using build 38 GRCh38 instead ", sep=""))
    genes_df <- toprdata::ENSGENES %>% dplyr::filter(gene_symbol %in% genes)
  }
  
  if(! is.null(chr)){
    chr <- gsub('chr','',chr)
    genes_df <- genes_df %>% dplyr::filter(chrom == paste("chr",chr,sep=""))
  }
  genes_df <- genes_df %>% dplyr::mutate(POS = gene_start + round((gene_end-gene_start)/2)) %>% dplyr::rename(Gene_Symbol=gene_symbol,CHROM=chrom)
  dist_genes <- genes_df %>% dplyr::distinct(CHROM,gene_start,gene_end, .keep_all = T)
  return(dist_genes)
}


#' Get the genetic position of a gene by gene name
#'
#' @description
#'
#' \code{get_gene_coords()} Get the gene coordinates for a gene
#' Required parameter is gene name
#'
#' @param gene_name A string representing a gene name (e.g. "FTO")
#' @param chr A string, search for the genes on this chromosome only, (e.g chr="chr1")
#' @param build A string, genome build, choose between builds 37 (GRCh37) and 38 (GRCh38) (default is 38)
#' @return Dataframe with the gene name and its genetic coordinates
#' @export
#'
#' @examples
#' \dontrun{
#' get_gene_coords("FTO")
#' }
#'

get_gene_coords  <- function(gene_name,chr=NULL, build=38){
  if(! is.null(chr)){
    chr <- gsub("chr", "", chr)
    chr <- paste("chr",chr,sep="")
    if(build=="38")
      gene <- toprdata::ENSGENES %>% dplyr::filter(chrom == chr) %>% dplyr::filter(gene_symbol == gene_name)
    else if(build=="37")
      gene <- toprdata::ENSGENES_37 %>% dplyr::filter(chrom == chr) %>% dplyr::filter(gene_symbol == gene_name)
  }
  else{
    if(build=="38")
      gene <- toprdata::ENSGENES %>% dplyr::filter(gene_symbol == gene_name)
    else if(build=="37")
     gene <- toprdata::ENSGENES_37 %>% dplyr::filter(gene_symbol == gene_name)
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
#' \code{get_top_snp()} Get the top hit from the dataframe
#' All other input parameters are optional
#'
#' @param df Dataframe containing association results
#' @param chr String, get the top hit in the data frame for this chromosome. If chromosome is not provided, the top hit from the entire dataset is returned.

#' @return Dataframe containing the top hit
#' @export
#'
#' @examples
#' \dontrun{
#' get_top_snp(CD_UKBB, chr="chr1")
#' }
#'
get_top_snp <- function(df, chr=NULL){
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
  return(c("darkblue","#E69F00","#00AFBB","#999999","#FC4E07","darkorange1","darkgreen","blue","red","magenta","skyblue","grey40","grey60","yellow","black","purple","orange","pink","green","cyan"))
}

set_genes_pos_adj <- function(genes, offsets){
  #CHROM,POS,Gene_Symbol
  genes <- genes %>% dplyr::mutate(CHROM=gsub('chr','',CHROM))
  genes[genes$CHROM=='X', 'CHROM'] <- "23"
  genes$CHROM <- as.integer(genes$CHROM)
  genes <- genes %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
  return(genes)
}

#' Get SNPs/variants within region
#'
#' @description
#'
#' \code{get_snps_within_region()}
#'
#' @param df data frame of association results with the columns CHR and POS
#' @param region A string representing the genetic region (e.g chr16:50693587-50734041)
#' @param chr A string, chromosome (e.g. chr16)
#' @param xmin An integer, include variants with POS larger than xmin
#' @param xmax An integer, include variants with POS smaller than xmax
#' @param keep_chr Deprecated: Logical, set to FALSE to remove the "chr" prefix before each chromosome if present (TRUE by default) 
#' @return the variants within the requested region
#' @export
#'
#' @examples
#' \dontrun{
#' get_snps_within_region(CD_UKBB, "chr16:50593587-50834041")
#' }
#'

get_snps_within_region <- function(df, region, chr=NULL, xmin=NULL, xmax=NULL,keep_chr=NULL){
  snps <- NULL
  if (!missing(keep_chr)) deprecated_argument_msg(keep_chr)
  df <- dat_column_chr_pos_check(df)
  if(!is.null(region)){
    tmp <- unlist(stringr::str_split(region, ":"))
    chr <- tmp[1]
    tmp_pos <- unlist(stringr::str_split(tmp[2], "-"))
    xmin <- as.numeric(tmp_pos[1])
    xmax <- as.numeric(tmp_pos[2])
  }
  if(!is.null(chr) & !is.null(xmin) & !is.null(xmax)){
    snps <- df %>% filter(CHROM == chr & POS >= xmin & POS <= xmax)
    if(length(snps$POS) < 1){
      print(paste("There are no SNPs within the region: ", chr, ":", xmin, "-", xmax, sep=""))
    }
  }
  return(snps)
}
