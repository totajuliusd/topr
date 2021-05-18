filter_on_chr=function(dat,chr){
  chr=gsub("chr", "", chr)
  for(i in 1:length(dat)){
    dat[[i]]=dat[[i]] %>% filter(CHROM==chr)
  }
  return(dat)
}

filter_on_xmin_xmax=function(dat,xmin=NULL,xmax=NULL){
  for(i in 1:length(dat)){
    if(! is.null(xmin)) {
      if(xmin < 0) xmin=0
      dat[[i]] = dat[[i]] %>% filter(POS >= xmin)
    }
    if(! is.null(xmax)) dat[[i]] = dat[[i]] %>% filter(POS <= xmax)
  }
  return(dat)
}

add_sign_thresh_to_plot=function(p1, df,sign_thresh,sign_thresh_color,xmin,xmax,ymin,ymax,overviewplot=FALSE){
  #add significance threholds to the plot
  if(is.null(xmin) || is.null(xmax)) {
    xmin=0
    xmax=max(df$POS)
  }
  xlabelpos=(xmax-xmin)*0.015
  ylabelpos=(ymax-ymin)*0.02
  #set some restrictions, since this function is used by region, chrom and manhattan functions
  if(xlabelpos > 20000) xlabelpos=20000
  if(ylabelpos > 1) ylabelpos=1

  #A thinner line is added to the overview plot
  size=ifelse(overviewplot, 0.2, 0.5)
  if(is.vector(sign_thresh)){
    for(i in 1:length(sign_thresh)){
      color="red"
      if(!is.null(sign_thresh_color) & (length(sign_thresh) == length(sign_thresh_color))){
        color=sign_thresh_color[i]
      }
      size=ifelse(overviewplot, 0.2, 0.5)

      p1=p1+geom_hline(yintercept = -log10(as.numeric(sign_thresh[[i]])),size=size, color=color, linetype="dashed")
      label=sign_thresh[[i]]
      ypos= -log10(sign_thresh[[i]])
      if(! overviewplot)
        p1=p1+geom_text(aes(x=xmin+xlabelpos,y=(ypos+ylabelpos),label=label),color=color,size=3.5)
    }
  }
  else{
    p1=p1+geom_hline(yintercept = -log10(as.numeric(sign_thresh[[i]])),size=size, colour="red", linetype="dashed")
    if( ! overviewplot)
      p1=p1+geom_text(aes(x=xmin+xlabelpos,y=(-log10(sign_thresh)+ylabelpos),label=sign_thresh),color=color,size=3.5)
  }

  return(p1)
}

include_chrX=function(dat){
  chrX=0
  for(i in 1:length(dat)){
    if("23" %in% unique(dat[[i]]$CHROM)){
      return(1)
    }
  }
  return(chrX)
}

is_df_empty=function(df, type){
  if(is.null(type)) type="main data "
  if(is.data.frame(df)) df=list(df)
  for(i in 1:length(df)){
    datrows=nrow(df[[i]])
    if(datrows == 0){
      stop(paste("One or more of the dataframes provided as input is empty!!! (", type," n=",i,")",sep=""))
    }
  }
}

color_genes=function(p1,dat,genes,genes_color,genes_ypos){
  gcol="green"
  if(! is.null(genes_color))
    gcol=genes_color[1]
  if(! is.null(genes)){
    for(i in 1:length(genes)){
      gene=genes[i]
      if(! is.null(genes_color[i]) & length(genes_color) >= i)
        gcol=genes_color[i]
      for(j in 1:length(dat)){
        df_gene=dat[[j]] %>% filter(Gene_Symbol == gene)
        if(length(df_gene$POS) < 1){
          warning(paste("NO datapoints in the dataframe for gene ",gene , sep=""))
        }
        else{
          if(j==1 & (!"log10p" %in% colnames(df_gene))){
            df_gene$log10p=-log10(df_gene$P)
          }
          if(j==1 & (!"pos_adj" %in% colnames(df_gene))){ #function called from chromplot
            df_gene$pos_adj=df_gene$POS
          }
          p1=p1+geom_point(data=df_gene, aes(x=pos_adj, y=log10p),color=gcol, size=2, shape=df_gene$shape)
          df_gene_label=df_gene %>% arrange(P) %>% head(n=1)
          if(j==1){
            p1=p1+geom_text_repel(data=df_gene_label, aes(x=pos_adj,y=log10p,label=Gene_Symbol), color="black",size=3,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5)
            #p1=p1+geom_text_repel(data=df_gene_label, aes(x=pos_adj,y=2.5,label=Gene_Symbol), force_pull=0,color="black",size=3,direction="x",angle=75,hjust=0,nudge_y =0.05,max.iter = 1e4, max.time = 1,segment.size=0.2,segment.alpha =.5)
          }
        }
      }
    }
  }
  return(p1)
}


add_genes2manhattan=function(p1,dat,offsets,genes,genes_color,genes_ypos){
  gcol="green"
  if (! (is.data.frame(genes) & ("CHROM" %in% colnames(genes)))){
    genes=get_genes_by_Gene_Symbol(genes)
  }
  #returns the genes like this:
  #CHROM,POS,Gene_Symbol
print(paste("THIIS IS THE YPOS ",genes_ypos,sep=""))
  genes=genes %>% dplyr::mutate(CHROM=gsub('chr','',CHROM))
  genes[genes$CHROM=='X', 'CHROM']="23"
  genes$CHROM = as.integer(genes$CHROM)
  genes=genes%>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
  if(is.null(genes_ypos)){
    genes_ypos=1
  }
  p1=p1+geom_point(data=genes, aes(x=pos_adj, y=genes_ypos),color=gcol, size=2, shape=15)

  p1=p1+geom_text_repel(data=genes, aes(x=pos_adj,y=genes_ypos,label=Gene_Symbol), color="black",size=3,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5)
  #p1=p1+geom_text_repel(data=df_gene_label, aes(x=pos_adj,y=2.5,label=Gene_Symbol), force_pull=0,color="black",size=3,direction="x",angle=75,hjust=0,nudge_y =0.05,max.iter = 1e4, max.time = 1,segment.size=0.2,segment.alpha =.5)

  return(p1)
}



add_variants_from_variant_list=function(p1,variant_list,annotation_thresh,rectr,show_xaxis=FALSE,offsets=NULL){
  for(i in 1:length(variant_list)){
    variants=variant_list[[i]]$v
    ##variants=set_chrs_ymax_and_shape(variants) -- TODO update
    variants=set_gene_symbol(variants)

    if(!is.null(offsets)){
      variants=variants %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
      variant_list[[i]]$alpha=1
      p1=p1+geom_point(data=variants, aes(x=pos_adj, y=log10p),  alpha=ifelse(!(is.null(variant_list[[i]]$alpha)), variant_list[[i]]$alpha,variants$alpha), color=variant_list[[i]]$col, size=ifelse(!(is.null(variant_list[[i]]$size)), variant_list[[i]]$size,variants$size), shape=ifelse(!(is.null(variant_list[[i]]$shape)), variant_list[[i]]$shape, variants$shape))

      if(i == 1){

        p1=p1+suppressWarnings(geom_text_repel(data=variants, aes(x=pos_adj,y=log10p,label=ifelse(P<annotation_thresh, Gene_Symbol,"")), color="grey40",size=3,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5))

      }
    }
    else{
      p1=p1+geom_point(data=variants, aes(x=POS, y=log10p),  alpha=ifelse(!(is.null(variant_list[[i]]$alpha)), variant_list[[i]]$alpha,variants$alpha), color=variant_list[[i]]$col, size=ifelse(!(is.null(variant_list[[i]]$size)), variant_list[[i]]$size,variants$size),
                       shape=ifelse(!(is.null(variant_list[[i]]$shape)), variant_list[[i]]$shape, variants$shape))
      if(i==1){

        if(show_xaxis){
          p1=p1+suppressWarnings(geom_text_repel(data=variants, aes(x=POS,y=log10p,label=ifelse(P<annotation_thresh, Gene_Symbol,"")), color="grey40",size=3.5,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5))
        }
        else{
          p1=p1+suppressWarnings(geom_text_repel(data=variants, aes(x=POS,y=log10p,label=ifelse(P<annotation_thresh, ID,"")), color="grey40",size=3.5,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5))
        }
      }
    }
  }
  return(p1)
}

set_gene_symbol=function(variants){
  if(! "Gene_Symbol" %in% colnames(variants)){
    if("ID" %in% colnames(variants))
      variants$Gene_Symbol=variants$ID
    else
      variants$Gene_Symbol=paste(variants$CHROM,variants$POS,sep=":")
  }
  variants$Gene_Symbol=ifelse(variants$Gene_Symbol==".", variants$ID, variants$Gene_Symbol)
  return(variants)
}

tidy_plot=function(p1,axis_text_size=12,axis_title_size=12, title_text_size=14,legend_title_size=12,legend_text_size=12){
  #remove space between axis and plots
  p1=p1+scale_y_continuous(expand=c(.02,.02))
  p1=p1+scale_x_continuous(expand=c(.01,.01),labels = scales::comma )

  #remove the vertical and horizontal lines from the plot. Set the size of the axis and title texts
  p1=p1+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              axis.text=element_text(size=axis_text_size),
              axis.title = element_text(size=axis_title_size),
              title = element_text(size=title_text_size),
              legend.title=element_text(size=legend_title_size),
              legend.text=element_text(size=legend_text_size))

  return(p1)
}

format_table=function(dat, digits_colnames, sci_colnames){
  dat=format_digits(dat,digits_colnames)
  dat=format_sci(dat,sci_colnames)
  return(dat)
}

format_digits=function(dat, colnames, ndigits=3){
  for (cname in colnames){
    dat[,cname]=as.numeric(formatC(unlist(dat[,cname]), digits=ndigits))
  }
  return(dat)
}
format_sci=function(dat, colnames,ndigits=1){
  for (cname in colnames){
    dat[,cname]=as.numeric(formatC(unlist(dat[,cname]), format="e", digits=ndigits))
  }
  return(dat)
}