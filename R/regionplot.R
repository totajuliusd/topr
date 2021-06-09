
#' Region plot
#'
#' @description
#'
#' @description
#'
#' \code{regionplot()} displays the association results for a smaler region within one chromosome
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase), xmin and xmax
#'
#' All other input parameters are optional
#'
#' @param show_overview Show the overview plot (default show_overview=TRUE)
#' @param genes List of genes or exons to plot
#' @param show_genes Show genes instead of exons (default show_genes=FALSE)
#' @param show_exons Show exons instead of genees (default show_exons=FALSE)#'
#' @inheritParams chromplot
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' regionplot(df,xmin=xmin,xmax=xmax)
#' }


regionplot=function(dat, annotation_thresh=NULL, title="",label_all=0,xmin=0, size=2, shape=19, alpha=1,label_size=4,
                    color=c("darkblue","#E69F00","#00AFBB","#999999","#FC4E07","darkorange1"),
                    axis_text_size=11,axis_title_size=12,title_text_size=13,legend_title_size=12,legend_text_size=12,
                    show_genes=FALSE, show_overview=TRUE, show_exons=FALSE,sign_thresh=5.1e-9,sign_thresh_color="red", variant_list=NULL,genes=NULL,variants=NULL,variant_ids=NULL,
                    variant_ids_color=NULL,xmax=NULL,ymin=NULL,ymax=NULL,
                    chr=NULL,vline=NULL,legend_name="Data:",legend_position="right",legend_labels=NULL,gene=NULL,highlight_genes=NULL,highlight_genes_color=NULL){
  # three plots, overview_plot, main_plot and gene_plot
  #only include overview plot if df region is larger than the region between xmin and xmax
  #dat=list(gwas_ukbb,gwas_abbvie)

  #is_df_empty(dat,"main dat")

  #is_df_empty(variants, "variants")

  if(! is.null(gene)){
    gene_df=get_gene_coords(gene,chr)
    if(dim(gene_df)[1]==0){
      stop(paste("Could not find gene ",gene))
    }else{
      chr=gene_df$chrom
      xmin=gene_df$gene_start-100000
      xmax=gene_df$gene_end+100000
    }
  }
  if(is.null(xmin) || is.null(xmax))
    stop("The region to view has not been set. Set the size of the region with xmin and xmax, e.g. xmin=66716513, xmax=67716513")
  #check and set data
  if(is.data.frame(dat)) dat=list(dat)
  dat=dat_column_check_and_set(dat)
  if(is.null(chr)) chr=get_chr_from_df(dat[[1]])
  dat=dat_chr_check(dat)
  df=dat[[1]]


  dat=set_size_shape_alpha(dat,size,shape,alpha)
  dat=set_color(dat,color)
  dat=filter_on_chr(dat,chr)
  dat_full=dat
  dat=filter_on_xmin_xmax(dat,xmin,xmax)
  #check and set variants

  if(! is.null(variants)){
    if(is.data.frame(variants)) variants=list(variants)


    variants=dat_column_check_and_set(variants)
    variants=dat_chr_check(variants)
    variants=filter_on_chr(variants,chr)
    variants=filter_on_xmin_xmax(variants,xmin,xmax)
    variant_color=color
    if(length(dat) == 1) variant_color="red"
    variants=set_size_shape_alpha(variants,size,shape,alpha)
    variants=set_color(variants,variant_color)
    variants=set_annotation_thresh(variants,annotation_thresh)
  }

  main_plot=chromplot(dat, variants=variants,variant_list=variant_list,annotation_thresh=annotation_thresh,xmin=xmin,xmax=xmax,ymin=ymin, ymax=ymax,show_xaxis = F,size=size,
                      shape=shape,alpha=alpha, color=color,.checked=TRUE, label_size=label_size, variant_ids=variant_ids,variant_ids_color=variant_ids_color,
                      vline=vline,legend_name=legend_name,legend_position = legend_position,legend_labels=legend_labels,sign_thresh=sign_thresh,
                      sign_thresh_color=sign_thresh_color,chr=chr,highlight_genes=highlight_genes,highlight_genes_color=highlight_genes_color,
                      axis_text_size=axis_text_size,axis_title_size=axis_title_size,title_text_size=title_text_size,legend_title_size=legend_title_size,legend_text_size=legend_text_size)


  if(show_overview){
    if(min(df$POS) < xmin & max(df$POS) > xmax){ # showOverviewplot if TRUE
      overview_plot=mk_overview_plot(dat_full,xmin=xmin,xmax=xmax,chr=chr,sign_thresh = sign_thresh,sign_thresh_color = sign_thresh_color)
    }else{
      show_overview=FALSE
    }
  }

  if(is.null(genes)){
    if(show_genes){
      genes=get_genes(chr,xmin,xmax)
    }else if(show_exons){
      genes=get_exons(chr,xmin,xmax)
    }else{
      if(xmax-xmin < 1000001){
        show_exons=TRUE
        genes=get_exons(chr,xmin,xmax)
      }
      else{
        show_genes=TRUE
        genes=get_genes(chr,xmin,xmax)
      }
    }
  }

  if(show_exons & ("exon_chromstart" %in% colnames(genes))){
    gene_plot=exonplot(genes, xmin, as.numeric(xmax),label_size=label_size,vline=vline)
  }else{
    gene_plot=geneplot(genes, xmin, as.numeric(xmax),label_size=label_size,vline=vline)
  }
  gene_plot=tidy_plot(gene_plot,axis_text_size=axis_text_size,axis_title_size=axis_title_size,title_text_size=title_text_size,legend_title_size=legend_title_size,legend_text_size=legend_text_size)
  gene_plot=gene_plot+scale_x_continuous(expand=c(.01,.01),labels = scales::comma )

  print(paste("Zoomed to region: ",chr,":",xmin,"-",xmax,sep=""))

  #and plot using egg
  #https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.htm

  g2=ggplotGrob(main_plot)
  g3=ggplotGrob(gene_plot)
  fg2=egg::gtable_frame(g2, height=unit(7, "null"),debug=F)
  fg3=egg::gtable_frame(g3, height=unit(2, "null"),debug=F)

  if(show_overview){
    g1=ggplotGrob(overview_plot)
    fg1=egg::gtable_frame(g1, height=unit(1.25, "null"), debug=F)
    fg12=egg::gtable_frame(gridExtra::gtable_rbind(fg1,fg2,fg3),debug=FALSE)
  }else{
    fg12=egg::gtable_frame(gridExtra::gtable_rbind(fg2,fg3),debug=FALSE)
  }
  grid::grid.newpage()
  grid::grid.draw(fg12)
  #return(grid.draw(fg12))
}
