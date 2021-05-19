
set_size_shape_alpha=function(df,size,shape,alpha){
  df=set_size(df,size)
  df=set_shape(df,shape)
  df=set_alpha(df,alpha)
  df=set_shape_by_impact(df)
  return(df)
}

set_shape_by_impact=function(dat){
  for(i in 1:length(dat)){
    if("Max_Impact" %in% colnames(dat[[i]])){
      dat[[i]]$shape0=ifelse(dat[[i]]$Max_Impact == "MODERATE", 17, dat[[i]]$shape)
      dat[[i]]$shape=ifelse(dat[[i]]$Max_Impact =="HIGH", 8, dat[[i]]$shape0)
    }
  }
  return(dat)
}

set_annotation_thresh=function(dat,annotation_thresh){
  if(! is.null(annotation_thresh)){
    for(i in 1:length(dat)){
      if(!is.null(annotation_thresh) & nrow(dat[[i]]>0)){
        if(is.vector(annotation_thresh) & (length(annotation_thresh) >= i))
          dat[[i]]$annotation_thresh=annotation_thresh[i]
        else
          dat[[i]]$annotation_thresh=annotation_thresh
      }
    }
  }
  return(dat)
}

set_color=function(dat,color){
  if(length(color) < length(dat))
    stop(paste("There are ",length(dat), " datasets, but only ",length(color), " color. Add more colors with the color argument, eg. color=c(\"blue\",\"yellow\",\"green\", etc)"))
  for(i in 1:length(dat)){
    if(!is.null(color) & nrow(dat[[i]]>0)){
      if(is.vector(color) & (length(color) >= i))
        dat[[i]]$color=color[i]
      else
        dat[[i]]$color=color
    }
  }
  return(dat)
}

set_size=function(dat,size){
  if(is.null(size))
    size=1
  for(i in 1:length(dat)){
    if(!is.null(size) & nrow(dat[[i]]>0)){
      if(is.vector(size) &  length(size) >= i )
        dat[[i]]$size=size[i]
      else
        dat[[i]]$size=size
    }
  }
  return(dat)
}

set_alpha=function(dat,alpha){
  if(is.null(alpha))
    alpha=0.7
  for(i in 1:length(dat)){
    if(!is.null(alpha) & nrow(dat[[i]]>0)){
      if(is.vector(alpha) & (length(alpha) >= i))
        dat[[i]]$alpha=alpha[i]
      else
        dat[[i]]$alpha=alpha
    }
  }
  return(dat)
}

set_shape=function(dat,shape){
  if(is.null(shape))
    shape=19
  for(i in 1:length(dat)){
    if(!is.null(shape) & nrow(dat[[i]]>0)){
      if(is.vector(shape) & (length(shape) >= i))
        dat[[i]]$shape=shape[i]
      else
        dat[[i]]$shape=shape
    }
  }
  return(dat)
}

set_log10p=function(dat,ntop){
  for(i in 1:length(dat)){
    df=dat[[i]]
    if(ntop > length(dat)){ ntop=length(dat)}
    if(i <= ntop){
      dat[[i]]=df %>% dplyr::mutate(log10p=-log10(P))

    }else{
      dat[[i]]=df %>% dplyr::mutate(log10p=log10(P))
    }
  }
  return(dat)
}

## Getters

get_exons=function(chr,xmin,xmax){
  chr=gsub("chr", "", chr)
  cand_region=paste("chr",chr,":",xmin,"-",xmax,sep="")
  exons=query(paste("gor -p ",cand_region," ref/ensgenes/ensgenes.gorz | map -c gene_symbol <(nor -h ref/ensgenes/ensgenes_exons.gorz | select gene_symbol,chromstart,chromend,exon | rename chromstart exon_chromstart | rename chromend exon_chromend)",sep=""))
  return(exons)
}

get_gene_coords=function(gene_symbol,chr=NULL){
  if(! is.null(chr)){
    chr=gsub("chr", "", chr)
    gene=query(paste("gor -p chr",chr," ref/ensgenes/ensgenes.gorz | where gene_symbol= \"",gene_symbol,"\" | select chrom,gene_start,gene_end,gene_symbol", sep=""))
  }
  else{
    gene=query(paste("pgor ref/ensgenes/ensgenes.gorz | where gene_symbol= \"",gene_symbol,"\" | select chrom,gene_start,gene_end,gene_symbol", sep=""))
  }
  return(gene)
}

get_genes=function(chr,xmin,xmax){
  chr=gsub("chr", "", chr)
  cand_region=paste("chr",chr,":",xmin,"-",xmax,sep="")
  genes=query(paste("gor -p ",cand_region," ref/ensgenes/ensgenes.gorz | select chrom,gene_start,gene_end,gene_symbol,gene_stable_id,biotype,strand", sep=""))
  return(genes)
}

get_genes_by_Gene_Symbol=function(genes){
  gene_string=paste(genes, collapse = '","')
  genes=query(paste("gor ref/ensgenes/ensgenes.gorz | where gene_symbol in (\"",gene_string,"\") | select chrom,gene_start,gene_end,gene_symbol | calc POS gene_start+ round((gene_end-gene_start)/2) | rename chrom CHROM |rename gene_symbol Gene_Symbol | select CHROM,POS,Gene_Symbol",sep=""))
  return(genes)
}

get_snps_by_ID=function(snps, file=NULL, id_col_name="ID"){
  snps_string=paste(snps, collapse = '","')
  snps_df=query(paste("gor ",file,"| where ",id_col_name," in (\"",snps_string,"\") ",sep=""))
  return(snps_df)
}

get_chr_from_df=function(df){
  #take the chromosome from the dataframe
  if(length(unique(df$CHROM)) == 1)
    chr=unique(df$CHROM)
  else
    stop("More than one chromosome in the input data, please select which chromosome to display by setting the chr argument (e.g. chr='chr10') or use input data with just one chromosome")
  return(chr)
}


get_shades=function(offsets,dat,ntop,include_chrX){
  n_offsets=11
  if(!(include_chrX)) n_offsets=10

  y1=c(rep(0, n_offsets))  #if there is no bottom plot
  ymin=get_ymin(ntop,dat)

  if(length(dat) > ntop){
    y1=c(rep(ymin, n_offsets))
  }
  ymax=get_ymax(ntop,dat)
  gene_label_size=2

  ##### here we need to accoutn for whether we have chr X or not
  if(include_chrX){
    shades=data.frame(x1=c(offsets[[2]],offsets[[4]],offsets[[6]],offsets[[8]],offsets[[10]],offsets[[12]],offsets[[14]],offsets[[16]],offsets[[18]],offsets[[20]],offsets[[22]]),
                      x2=c(offsets[[3]],offsets[[5]],offsets[[7]], offsets[[9]],offsets[[11]],offsets[[13]],offsets[[15]],offsets[[17]],offsets[[19]],offsets[[21]],offsets[[23]]),
                      y1=y1,
                      y2=c(rep(ymax, n_offsets)))
  }
  else{
    shades=data.frame(x1=c(offsets[[2]],offsets[[4]],offsets[[6]],offsets[[8]],offsets[[10]],offsets[[12]],offsets[[14]],offsets[[16]],offsets[[18]],offsets[[20]]),
                      x2=c(offsets[[3]],offsets[[5]],offsets[[7]], offsets[[9]],offsets[[11]],offsets[[13]],offsets[[15]],offsets[[17]],offsets[[19]],offsets[[21]]),
                      y1=y1,
                      y2=c(rep(ymax, n_offsets)))
  }
  return(shades)
}
get_ymax=function(ntop,dat){
  ymax=0
  if(length(dat)<ntop) ntop=length(dat)
  for(i in 1:ntop){
    tmp=max(dat[[i]]$log10p)
    if(tmp>ymax){
      ymax=tmp
    }
  }
  return(ymax)
}

get_ymin=function(ntop,dat){
  ymin=0
  if(ntop < length(dat)){
    for(i in ntop:length(dat)){
      tmp=min(dat[[i]]$log10p)
      if(tmp<ymin){
        ymin=tmp
      }
    }
  }
  return(ymin)
}

get_chr_offsets=function(include_chrX=1){
  #create the offsets from the crh lenghts file
  chr_lengths=read.table(file=chr_lengths_file,as.is=T)
  chr_lengths$CHROM=gsub("chr","", chr_lengths$V1)
  chr_lengths=chr_lengths %>% filter(! V1 %in% c("chrM") )
  chr_lengths[chr_lengths$CHROM=="X",'CHROM']="23"
  chr_lengths[chr_lengths$CHROM=="Y",'CHROM']="24"
  chr_lengths$CHROM=as.integer(chr_lengths$CHROM)
  no_chrs=22
  if(include_chrX)
    no_chrs=23
  chr_lengths = chr_lengths %>% filter(CHROM< no_chrs +1 )
  tmp=suppressMessages(chr_lengths %>% group_by(CHROM) %>% summarize(m=V2) %>% mutate(offset=cumsum(lag(m, default=0))))
  #offsets=data.frame(offsets=tmp$offset, CHROM=tmp$CHROM)
  offsets=setNames(tmp$offset,tmp$CHROM)
  return(offsets)
}

get_ticknames=function(df){
  no_chrs=ifelse("chrX" %in% df$CHROM || "X" %in% df$CHROM || "chr23" %in% df$CHROM || "23" %in% df$CHROM, 23, 22)
  if(no_chrs == 23){
    ticknames <- c(1:16, '',18, '',20, '',22, 'X')
  }else{
    ticknames <- c(1:16, '',18, '',20, '',22)
  }
  tickpos <-suppressMessages( df %>%
                                group_by(CHROM)%>%
                                summarize(pm=mean(pos_adj))%>%
                                dplyr::pull(pm))

  names(tickpos) <- NULL
  return(list(names=ticknames, pos=tickpos))
}

#' Get the top SNP per 10 MB
#'
#' @description
#'
#' \code{get_best_snp_per_MB()} Get the top SNP per 10 MB
#' All other input parameters are optional
#'
#' @param df Dataframe
#' @param thresh Threshold
#' @param region region

#' @return numeric vector
#' @export
#'
#' @examples
#' \dontrun{
#' data(gwas_CD)
#' get_best_snp_per_MB(gwas_CD, thresh = 1e-09, region = 10000000)
#' }
get_best_snp_per_MB=function(df,thresh=1e-09,region=1000000){
  if(! "CHROM" %in% colnames(df) || (! "POS" %in% colnames(df)) || (! "P" %in% colnames(df))){
    if(is.data.frame(df)) dat=list(df)
    dat=dat_column_check_and_set(dat)
    dat=dat_chr_check(dat)
    df=dat[[1]]
  }
  df$CHROM=gsub("chr","", df$CHROM)
  df=df %>% filter(P<thresh)
  df$tmp=NA
  for(row in 1:nrow(df)){
    df$tmp=round(df$POS/region)
  }
  lead_snps=df %>% group_by(CHROM,tmp) %>% arrange(P) %>% filter(P== min(P))
  if(length(lead_snps)== 0){
    print(paste("There are no SNPs with a p-value below ",thresh, " in the input. Use the [thresh] argument to adjust the threshold.",sep=""))
  }

  return(ungroup(lead_snps)%>% distinct(CHROM,POS, .keep_all = T))
}


get_legend<-function(p1){
  tmp <- ggplot_gtable(ggplot_build(p1))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  p1=p1+theme(legend.position="none")
  return(legend)
}


thin_data=function(df, breaks=100000){
  if(! "CHROM" %in% colnames(df) || (! "POS" %in% colnames(df)) || (! "P" %in% colnames(df))){
    if(is.data.frame(dat)) dat=list(df)
    dat=dat_column_check_and_set(dat)
    dat=dat_chr_check(dat)
    df=dat[[1]]
  }
  rowcount=length(df$POS)
  dat1=df %>% filter(P <= 1e-07)
  dat2=df %>% filter(P > 1e-07)
   dat2$tmp=NA
   dat1$tmp=NA
   for(row in 1:nrow(dat2)){
     dat2$tmp=round(dat2$POS/breaks)
   }
   dat_filt1=dat2 %>% group_by(CHROM,tmp) %>% arrange(P) %>% filter(P== min(P))
   #dat_filt2=dat2 %>% group_by(CHROM,tmp) %>% arrange(P) %>% filter(P== max(P))
   dat_filt=rbind(dat1,dat_filt1)
   #dat_filt=rbind(dat0,dat_filt2)
   print(paste("Data set was reduced from ",rowcount, " lines, to ",length(dat_filt$POS),sep=""))
  return(dat_filt)
}
