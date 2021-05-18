
#' Manhattan plot
#'
#' @description
#'
#' \code{manhattan()} displays the association results for the entire genome
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase)
#'
#' All other input parameters are optional
#'
#' @param ntop Number of datasets (GWASes) to show on the top plot
#' @inheritParams chromplot

#' @return plots using egg (https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html)
#' @examples
#'   manhattan(df)

manhattan=function(df, annotation_thresh=1e-09, ntop=3, title="", color=c("darkblue","#E69F00","#00AFBB","#999999","#FC4E07","darkorange1"),
                   sign_thresh=5e-09,label_size=3, size=1,shape=19,alpha=1,highlight_genes_color="green",
                   variants_size=1.2,variants_alpha=1,variants_shape=19,axis_text_size=11,axis_title_size=12, title_text_size=13, legend_title_size=12,legend_text_size=12, rect=NULL,
                   variant_list=NULL, variants=NULL,plot_order=NULL,legend_labels=NULL,legend_name=NULL,
                   ymin=NULL,ymax=NULL,variants_color=NULL,highlight_genes=NULL,highlight_genes_ypos=NULL,legend_position=NULL){
  #do a dataframe colname check here!!
  dat=df
  if(is.data.frame(dat)) dat=list(dat)
  dat=dat_column_check_and_set(dat)
  dat=dat_chr_check(dat)
  dat=set_size_shape_alpha(dat,size,shape,alpha)
  df=dat[[1]]

  if(! is.null(variants)){
    if(is.data.frame(variants)) variants=list(variants)
    variants=dat_column_check_and_set(variants)
    variants=dat_chr_check(variants)
    variants=set_size_shape_alpha(variants,variants_size,variants_shape,variants_alpha)
    variants=set_annotation_thresh(variants,annotation_thresh)

  }
  incl_chrX=include_chrX(dat)
  offsets=get_chr_offsets(incl_chrX)

  for(i in 1:length(dat)){
    dat[[i]]=dat[[i]] %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
  }
  if(! is.null(variants)){
    for(i in 1:length(variants)){
      variants[[i]]=variants[[i]] %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
    }
  }
  #in case there are multiple dataframes, use the one with most chromosome s
  dat4ticks=dat[[1]]
  for(i in 1:length(dat)){
    if(length(unique(dat[[i]]$CHROM))  > length(unique(dat4ticks$CHROM))){
      dat4ticks=dat[[i]]
    }
  }
  ticks=get_ticknames(dat4ticks)
  if(length(dat) > 1){
    dat=set_color(dat,color)
    if(! is.null(variants)){
      variants=set_color(variants,color)
      variants=set_log10p(variants,ntop)
    }
    dat=set_log10p(dat,ntop)
    shades=get_shades(offsets,dat,ntop,incl_chrX)
    if(is.null(legend_position)) legend_position="top"
    p1=manhattan_multi(dat, color=color,shades=shades,variants=variants, ntop=ntop, label_size=label_size, ymax=ymax,ymin=ymin,legend_name=legend_name,legend_labels=legend_labels)
  }
  else{
    if(is.null(variants_color))
       variants_color=c("red","green","skyblue","orange","magenta","yellow","darkorange","brown","pink")
    if(! is.null(variants)){
      variants=set_color(variants,variants_color)
    }
    p1=manhattan_single(dat[[1]],variants=variants,color=color,legend_name=legend_name,legend_labels=legend_labels,variants_color=variants_color)
    if(is.null(legend_position)) legend_position="bottom"
  }

  p1=p1+theme(legend.background=element_blank(),legend.key=element_rect(fill="white"),legend.position = legend_position,
              legend.text=element_text(size=12),panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank())

  if(length(variants) == "1" & is.null(legend_labels)){
    p1=p1 + theme(legend.position = "none")
  }
  else{
    p1=p1+theme(legend.background=element_blank(),legend.key=element_rect(fill="white"),legend.position = legend_position,
                legend.text=element_text(size=12),panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank())
  }


  thresh=as.numeric(sign_thresh)
  # TODO add_sign_thresh_to_plot
  # p1=add_sign_thresh_to_plot(p1,df, sign_thresh, sign_thresh_color,xmin,xmax,ymin,ymax)

  p1=p1+geom_hline(yintercept = -log10(thresh), colour="red", linetype="dashed")+
    geom_text(aes(x=20000000,y=(-log10(thresh)+1),label=thresh),color="red",size=3.5)+
    labs(title=title, x='Chromosome', y=expression(-log[10](italic(p)))) +
    theme(panel.border=element_blank(), axis.line=element_line(color="grey"))

  p1=tidy_plot(p1, axis_text_size=axis_text_size,axis_title_size=axis_title_size, title_text_size=title_text_size,legend_title_size=legend_title_size,legend_text_size=legend_text_size)

  p1=p1+scale_x_continuous(breaks=ticks$pos, labels=ticks$names,expand=c(.01,.01))

  if(! is.null(highlight_genes)){
    p1=add_genes2manhattan(p1,dat, offsets,highlight_genes, highlight_genes_color,highlight_genes_ypos)
  }
  if(! is.null(ymin) & is.null(ymax))
    p1=p1+coord_cartesian(ylim=c(ymin, get_ymax(ntop,dat)))
  if(! is.null(ymin) & (! is.null(ymax)))
    p1=p1+coord_cartesian(ylim=c(ymin, ymax))
  if(is.null(ymin) & (! is.null(ymax)))
    p1=p1+coord_cartesian(ylim=c(get_ymin(ntop,dat), ymax))


  return(p1)
}



manhattan_multi=function(dat, ntop=2, shades=shades, label_size=3, ymax=NULL,ymin=NULL,plot_order=NULL,color=NULL,variants=NULL,
                         legend_labels=NULL,legend_name="Datasets:"){
  #set and get ticknames and tickpositions -requires pos_adj to be included in the df, so cannot be called before the prepare_dat function
  # colorMap=mk_colorMap(dat,colors)
  p1=ggplot()+theme_bw() #+geom_point(data=dat[[1]]$gwas, aes(dat[[1]]$gwas$pos_adj, dat[[1]]$gwas$log10p),color=colors[1] alpha=0.7,size=1)+theme_bw()

  for(i in 1:length(dat)){
    p1=p1+geom_point(data=dat[[i]], aes(x=pos_adj, y=log10p, color=color), alpha=dat[[i]]$alpha, size=dat[[i]]$size,shape=dat[[i]]$shape)
    if(! is.null(variants)){
      for(i in 1:length(variants)){
        dat_labels= variants[[i]] %>% filter(P < annotation_thresh) %>% distinct(Gene_Symbol, .keep_all = T)
        nudge_y=4
        if(i>ntop)
          nudge_y=-4
        p1=p1+geom_text_repel(data=dat_labels, aes(x=pos_adj, y=log10p, label=Gene_Symbol),nudge_y=nudge_y,size=label_size,
                              segment.size=0.2, color=dat_labels$color,segment.color = "black",
                              direction="both",angle=0, vjust=0, max.iter=5000) # fontface=gene_labels_top$fontface,
      }
    }
  }
   if(length(dat)> ntop)
    p1=p1+geom_hline(yintercept = log10(5.1e-9), color="red", linetype="dashed")
  p1=p1+geom_rect(data=shades, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1),color="#D3D3D3", alpha=0.1, size=0.1)

  #add the legend for the different datasets
  if(! is.null(legend_name)) legend_name=legend_name
  if(is.null(legend_labels)){
    legend_labels=color[1:length(dat)]
    print("Use the legend_labels argument to change the legend labels from color names to meaningful labels! ")
  }
  p1=p1+scale_color_identity(guide = "legend", name=legend_name,  breaks=color[1:length(dat)],labels=legend_labels)
  return(p1)
}


manhattan_single=function(df, variants=variants,  variant_list=variant_list,label_size=label_size,offsets=offsets,
                          legend_name="Variants:",color=NULL,
                          legend_labels=NULL,variants_color=NULL){
  #df=df %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
  #set and get ticknames and tickpositions -requires pos_adj to be included in the df
  df$col_code=df$CHROM%%2
  df$color=ifelse(df$col_code=="1", "grey60","grey80")
  p1=ggplot(df)+geom_point(data=df, aes(x=pos_adj, y=-log10(P),color=color), size=df$size, shape=df$shape)+theme_bw()

  if(! is.null(variants)){
    if(length(variants)> 0){

      for(i in 1:length(variants)){
        p1=p1+suppressWarnings(geom_text_repel(data=variants[[i]], aes(x=pos_adj,y=-log10(P),label=ifelse(P<annotation_thresh, Gene_Symbol,"")), color="grey40",size=3,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5))
        p1=p1+geom_point(data=variants[[i]], aes(x=pos_adj, y=-log10(P),color=color), size=variants[[i]]$size, shape=variants[[i]]$shape)
      }
    }
  }
  #add the legend
  if(! is.null(legend_name)) legend_name=legend_name
  if(is.null(legend_labels)) legend_labels=variants_color[1:length(variants)]

  if(length(legend_labels < length(variants))){
      #TODO - label with color
  }

  if(length(variants) < length(legend_labels)) legend_labels=legend_labels[1:length(variants)]

  #p1=p1+scale_color_identity(guide = "legend",   breaks=color[1:length(variants)],)
  #if(is.null(legend_name)) legend_name="Variants:"
  p1=p1+scale_color_identity(guide = "legend",breaks=variants_color[1:length(variants)],name=legend_name,labels=legend_labels)
  return(p1)
}


manhattan_yval = function(df, yval, use_log10=TRUE, title="", label_size=3, rect=NULL, variant_list=NULL,variants=NULL){
  #df=set_chrs_ymax_and_shape(df)
  #create the offsets using the main data frame (df)
  tmp=suppressMessages(df %>% group_by(CHROM) %>% summarize(m=max(POS)) %>% mutate(offset=cumsum(lag(m, default=0))))
  offsets=setNames(tmp$offset, tmp$CHROM)
  #set the offsets for both df and variants
  df=df %>% dplyr::mutate(pos_adj=POS+offsets[CHROM])
  #set and get ticknames and tickpositions
  ticks=get_ticknames(df)
  df$col_code=df$CHROM%%2
  df$color=ifelse(df$col_code=="1", "grey60","grey80")

  if(use_log10){
    p1=ggplot(df)+geom_point(data=df, aes(x=pos_adj, y=-log10(df %>% pull(yval))), color=df$color, size=df$size, shape=df$shape)+theme_bw()
  }else{
    p1=ggplot(df)+geom_point(data=df, aes(x=pos_adj, y=df %>% pull(yval)), color=df$color, size=df$size, shape=df$shape)+theme_bw()
  }
  #remove space between axis and plots
  p1=p1+scale_y_continuous(expand=c(.02,.02))
  p1=p1+scale_x_continuous(breaks=ticks$pos, labels=ticks$names,expand=c(.01,.01))
  #subtitle = paste("Gene name is shown for variant if P<", sign_thresh, sep="")
  p1=p1+labs(title=title, x='Chromosome', y=yval) #,theme(axis.line = element_line(),panel.border = element_blank())
  p1=p1+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  p1=p1+theme(panel.border=element_blank(),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(color="grey"))
  p1=p1+theme(axis.text.y=element_text(size=12),axis.text=element_text(size=12),
              title = element_text(size=13))
  #p1=p1+theme(legend.background=element_blank(),legend.key=element_rect(fill="white"),legend.position = "top", legend.text=element_text(size=12))
  return(p1)
}

