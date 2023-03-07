#' Create a regionplot
#'
#' @description
#'
#' \code{regionplot()} displays the association results for a smaller genetic regions within one chromosome.
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase) and either a variant ID, gene name or the genetic region represented as a chromosome together with start and stop positions (either as a single string or as three separate arguments).
#'
#' All other input parameters are optional
#'
#' @param show_overview A logical scalar, shows/hides the overview plot (default= TRUE)
#' @param gene A string, the name of the  gene to zoom into (e.g. gene=FTO)
#' @param show_genes A logical scalar, show genes instead of exons (default show_genes=FALSE)
#' @param show_exons Deprecated : A logical scalar, show exons instead of genes (default show_exons=FALSE)
#' @param vline A number or vector of numbers to add a vertical line to the plot at a specific chromosomal position, e.g \code{vline=204000066}. Multiple values can be provided in a vector, e.g  \code{vline=c(204000066,100500188)}
#' @param max_genes An integer, only label the genes if they are fewer than max_genes (default values is 200).
#' @param show_gene_names A logical scalar, if set to TRUE, gene names are shown even though they exceed the max_genes count
#' @param variant A string representing the variant to zoom in on. Can be either an rsid, or a dataframe (with the columns CHROM,POS,P)
#' @param gene A string representing the gene to zoom in on (e.g. gene=FTO)
#' @param region A string representing a genetic region, e.g. chr1:67038906-67359979
#' @param gene_padding An integer representing size of the region around the gene, if the gene argument was used (default = 100000)
#' @param locuszoomplot A logical scalar set to FALSE. Only set to TRUE by calling the locuszoom function
#' @param show_gene_legend A logical scalar, set to FALSE to hide the gene legend (default value is TRUE)
#' @param unit_main the height unit of the main plot (default = 7)
#' @param unit_gene the height unit of the gene plot (default= 2 )
#' @param unit_overview the height unit of the overview plot (default = 1.25)
#' @param verbose Logical, set to FALSE to get suppress printed information
#' @param gene_color A string representing a color, can be used to change the color of the genes/exons on the geneplot
#' @param unit_ratios A string of three numbers separated by ":", for the overview, main and gene plots height ratios e.g  1.25:7:2
#' @param extract_plots Logical, FALSE by default. Set to TRUE to extract the three plots separately in a list
#' @inheritParams manhattan
#'
#' @return plots within ggplotGrobs, arranged with egg::gtable_frame
#' @export
#'
#' @examples
#' \dontrun{
#' regionplot(CD_UKBB, gene="IL23R")
#' }


regionplot <- function(df, ntop=10, annotate=NULL, xmin=0, size=2, shape=19, alpha=1,label_size=4, annotate_with="ID",
                       color=get_topr_colors(), axis_text_size=11,axis_title_size=12,title_text_size=13, show_genes=NULL, show_overview=TRUE,
                       show_exons=FALSE,max_genes=200, sign_thresh=5e-08, sign_thresh_color="red", sign_thresh_label_size=3.5,
                       xmax=NULL,ymin=NULL,ymax=NULL,protein_coding_only=FALSE,region_size=1000000,gene_padding=100000,angle=0,legend_title_size=12,legend_text_size=11,
                       nudge_x=0.01,nudge_y=0.01, rsids=NULL, variant=NULL,rsids_color=NULL,legend_name="",legend_position="right",
                       chr=NULL,vline=NULL,show_gene_names=NULL,legend_labels=NULL,gene=NULL, title=NULL, label_color=NULL,locuszoomplot=FALSE,
                       region=NULL,legend_nrow=NULL,gene_label_size=NULL, scale=1, show_legend=TRUE,sign_thresh_linetype="dashed", sign_thresh_size=0.5,
                       rsids_with_vline=NULL, annotate_with_vline=NULL,show_gene_legend=TRUE, unit_main=7, unit_gene=2, unit_overview=1.25, verbose=NULL,
                       gene_color=NULL,segment.size=0.2,segment.color="black",segment.linetype="solid", max.overlaps=10, unit_ratios=NULL, 
                       extract_plots=FALSE,label_fontface="plain",label_family="",gene_label_fontface="plain",gene_label_family="",build=38,
                       label_alpha=1, vline_color="grey",vline_linetype="dashed", vline_alpha=1,vline_size=0.5){
  # three plots, overview_plot, main_plot and gene_plot
  #only include overview plot if df region is larger than the region between xmin and xmax
  if (!missing(show_exons)) {
    deprecated_argument_msg(show_exons)
  }
  annot_with_vline <- FALSE
  if(! is.numeric(region_size)){
    region_size <- convert_region_size(region_size)
  }
  if(! is.null(annotate_with_vline)){
    annotate <- annotate_with_vline
    annot_with_vline <- TRUE
  }

  if(is.null(gene) & is.null(variant) & is.null(region) & (is.null(chr)  & is.null(xmax))){
    stop("Regional argument is missing (either gene, variant, region or the 3 arguments (chr,xmin,xmax) has to be provided as input to this function")
  }
  else{
    dat <- dat_check(df,verbose=verbose,locuszoomplot=locuszoomplot) %>% set_log10p(ntop)  %>% set_size_shape_alpha(size, shape, alpha, locuszoomplot = locuszoomplot)
   if(! locuszoomplot){
      dat <- dat %>% set_color(color)
    }
  
  top_snps <- NULL
  using_ntop <- FALSE
  if(length(dat) > ntop){
    using_ntop <- TRUE
  }
  if(! is.null(gene)){
    gene_df <- get_gene_coords(gene)
    if(dim(gene_df)[1]==0){  stop(paste("Could not find gene ",gene)) }
    else{
      chr <- gene_df$chrom
      xmin <- gene_df$gene_start-gene_padding
      xmax <- gene_df$gene_end+gene_padding
    }
  }else if(! is.null(variant)){
    variant_region=1000000
    if(is.data.frame(variant)){
      if("POS" %in% colnames(variant)){
        chr <- variant$CHROM
        xmin <- variant$POS-(variant_region/4)
        xmax <- variant$POS+(variant_region/4)
      }
      else{
        stop("There is no POS column in the variant dataframe")
      }
    }
    else{
      tmp <- dat[[1]] %>% filter(ID == variant) %>% arrange(P) %>% head(n=1)
      if(length(tmp$POS > 0 )){
        chr <- tmp$CHROM
        xmin <- tmp$POS-(variant_region/4)
        xmax <- tmp$POS+(variant_region/4)
      }
      else{
        stop(paste("Could not find the variant ",variant, " in the first input dataset " ),sep="")
      }

    }
  }else if(! is.null(region)){
    tmp <- unlist(stringr::str_split(region, ":"))
    chr <- tmp[1]
    tmp_pos <- unlist(stringr::str_split(tmp[2], "-"))
    xmin <- as.numeric(tmp_pos[1])
    xmax <- as.numeric(tmp_pos[2])
  }
  #use the first chromosome in the dataframe if no chromosome was given as argument
  chr <- ifelse(is.null(chr), unique(dat[[1]]$CHROM[1]), chr)
  if(! is.null(chr)){
    dat <- dat %>% filter_on_chr(chr)
    if(! is.null(xmin) & ! is.null(xmax)){
      df_full <- dat
      dat <- dat %>% filter_on_xmin_xmax(xmin,xmax)
    }
  }
  # get the annotation
  if(! is.null(annotate)){ top_snps <- get_annotation(dat, annotate=annotate, region_size=region_size,protein_coding_only = protein_coding_only,
                                                      verbose=verbose,nudge_x=nudge_x,nudge_y=nudge_y,angle=angle,
                                                      label_fontface=label_fontface, label_family=label_family, build=build,label_alpha=label_alpha) }
  # get the main plot
  main_plot <- get_base_plot(dat,color=color,legend_labels=legend_labels, show_legend = show_legend, legend_name = legend_name,legend_position = legend_position,locuszoomplot=locuszoomplot, legend_nrow=legend_nrow,scale=scale,verbose=verbose) %>%
    set_plot_text_sizes(axis_text_size=axis_text_size,axis_title_size = axis_title_size,legend_text_size=legend_text_size, legend_title_size=legend_title_size,scale=scale)
  
  main_plot <- main_plot+theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())
  if(! is.null(annotate)){
    main_plot <- main_plot %>%  add_annotation(plot_labels = top_snps,annotate_with=annotate_with,angle=angle,label_size = label_size, nudge_y = nudge_y, 
                                               nudge_x = nudge_x, label_color=label_color,scale=scale,
                                               segment.color=segment.color,segment.size=segment.size,segment.linetype=segment.linetype,max.overlaps=max.overlaps)
    if(annot_with_vline){
      main_plot <- main_plot %>% add_vline(top_snps$POS, vline_color=vline_color, vline_linetype = vline_linetype, vline_alpha=vline_alpha, vline_size=vline_size,scale=scale)
    }
  }
  if(locuszoomplot){
    #annotate the top variant
    if(is.null(rsids) & is.null(rsids_with_vline)){
      ld_snp <- get_main_LD_snp(dat, nudge_x=nudge_x, nudge_y=nudge_y, label_fontface=label_fontface, angle=angle, label_alpha=label_alpha)
      if(is.null(annotate) & !is.null(ld_snp) & length(ld_snp$POS)> 0){
        main_plot <- main_plot %>%  add_annotation(plot_labels = ld_snp,annotate_with=annotate_with,angle=angle,label_size = label_size, nudge_y = nudge_y, nudge_x = nudge_x, 
                                                   label_color=label_color,scale=scale,segment.color=segment.color,segment.size=segment.size,
                                                   segment.linetype=segment.linetype,max.overlaps=max.overlaps)
        if(annot_with_vline){
          main_plot <- main_plot %>% add_vline(ld_snp$POS, vline_color=vline_color, vline_linetype = vline_linetype, vline_alpha=vline_alpha, vline_size=vline_size,scale=scale)
        }
      }
    }
  }
  main_plot <- main_plot %>% add_vline(vline, vline_color=vline_color, vline_linetype = vline_linetype, vline_alpha=vline_alpha, vline_size=vline_size,scale=scale)
  main_plot <- main_plot + scale_y_continuous(expand=c(.02,.02)) + scale_x_continuous(expand=c(.01,.01))

if(!is.null(sign_thresh)){
 main_plot <- main_plot %>% add_sign_thresh(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, using_ntop = using_ntop,sign_thresh_linetype = sign_thresh_linetype, sign_thresh_size = sign_thresh_size, scale=scale) %>%
    add_sign_thresh_labels(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, xmin=xmin,sign_thresh_label_size = sign_thresh_label_size,scale=scale)
}
    with_vline <- FALSE
  if(! is.null(rsids_with_vline)){
    rsids <- rsids_with_vline
    with_vline <- TRUE
  }
    rsids_df <- NULL
  if( ! is.null(rsids)){
    rsids_df <- get_rsids_from_df(dat,rsids)
    main_plot <-main_plot %>% add_rsids(rsids_df, rsids_color=rsids_color, nudge_x=nudge_x, nudge_y=nudge_y, label_size=label_size, angle=angle, label_color=label_color, scale=scale ,
                                        with_vline = with_vline, label_fontface=label_fontface,label_family=label_family,segment.size=segment.size,segment.color=segment.color,segment.linetype=segment.linetype)
    
  }

  if(is.null(ymin)){ ymin <- get_ymin(dat)* 0.99}
  if(is.null(ymax)){ ymax <- get_ymax(dat)* 1.04}

  main_plot <- main_plot %>% set_xmin_xmax_ymin_ymax(xmin,xmax,ymin,ymax)
  main_plot <- main_plot + labs(y=expression(-log[10](italic(P))))

  #get the overview plot
  overview_plot <- get_base_plot(df_full,color=color,show_legend = FALSE,overview_plot = TRUE)
  overview_plot <- rm_axis_text_and_ticks(overview_plot)
  overview_plot <- add_zoom_rectangle(overview_plot,df_full,xmin, xmax)
  overview_plot <- overview_plot %>% add_sign_thresh(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, using_ntop = using_ntop,sign_thresh_linetype = sign_thresh_linetype, sign_thresh_size = 0.2,scale=scale)
  if(! is.null(title)){
   if(show_overview){
    #add the title to the overview plot
     overview_plot <- overview_plot %>% add_title(title=title, title_text_size = title_text_size)
   }else{ # add the title to the mainplot
     main_plot <- main_plot %>% add_title(title=title, title_text_size = title_text_size)
   }
   }
  if(using_ntop){
    main_plot <- main_plot %>%  add_zero_hline()
    overview_plot <- overview_plot %>% add_zero_hline()
  }
  if(is.null(show_genes)){
    if(xmax-xmin < 1000001){
      show_genes <- FALSE
    }else{show_genes <- TRUE}
  }
  ngenes <- 0

  #get the genes in the region
  genes <- get_genes_in_region(chr=chr,xmin=xmin,xmax=xmax,protein_coding_only=protein_coding_only,show_genes=show_genes, build=build)
    ngenes <- length(genes$gene_symbol)
    if(ngenes>100){show_gene_names <- FALSE}
    if(ngenes < max_genes){
      if(show_genes){
      plot_args <- get_geneplot_coords(genes,xmin,xmax)
      } else{ #plot with exons
        plot_args <- get_exonplot_coords(genes,xmin,xmax)
      }
    }
    else{
      print("Aborting!")
      print(paste("There are [",ngenes, "] genes within this region, which is too many for plotting in this view", sep=""))
      print(paste("The maximum number of genes allowed for plotting is set to: ", max_genes, ". (alhough not recommended, this limit can be modified using the max_genes argument). ",sep=""))
      print("The number of genes can be reduced by either setting the protein_coding_only argument to TRUE (protein_coding_only=TRUE), or by reducing the gene_padding parameter (e.g. gene_padding=500)")
      stop("Aborting due to too many genes within the region!!")
    }
    if(is.null(verbose)){
      print(paste("Zoomed to region:  ",chr,":",xmin,"-",xmax,sep=""))
    }else if(verbose){
      print(paste("Zoomed to region:  ",chr,":",xmin,"-",xmax,sep=""))
    }
    if(! is.null(gene_label_size)){
      label_size<-gene_label_size
    }
    if(annot_with_vline){
      vline <- top_snps$POS
    }
    if(with_vline){
      vline <- rsids_df$POS
    }
    gene_plot <- get_gene_plot(plot_args,label_size=label_size,xmin = xmin, xmax = xmax, show_gene_names = show_gene_names,scale=scale, 
                               show_gene_legend=show_gene_legend,gene_color=gene_color,gene_label_fontface = gene_label_fontface,gene_label_family = gene_label_family)
  gene_plot <- gene_plot + scale_y_continuous(expand=c(.02,.02)) + scale_x_continuous(expand=c(.01,.01))
  gene_plot <- gene_plot + xlab(paste("Position on Chromosome ",gsub('chr','',chr), sep=""))
  gene_plot <- gene_plot %>% add_vline(vline, vline_color=vline_color, vline_linetype = vline_linetype, vline_alpha=vline_alpha, vline_size=vline_size,scale=scale)
  gene_plot <- gene_plot %>% set_legend_texts(legend_title_size=legend_title_size, legend_text_size=legend_text_size,scale=scale)
  gene_plot <- gene_plot %>% set_plot_text_sizes(axis_text_size=axis_text_size,axis_title_size = axis_title_size,legend_text_size=legend_text_size, legend_title_size=legend_title_size,scale=scale)
  plots <- list(main_plot=main_plot, overview_plot=overview_plot, gene_plot=gene_plot)
  if(extract_plots){
    return(plots)
  }else{
    plots %>% draw_plots_with_egg(show_overview = show_overview, unit_main=unit_main, unit_gene=unit_gene, unit_overview=unit_overview,unit_ratios=unit_ratios)
  }
  }
}


draw_plots_with_egg <- function(plots, show_overview=TRUE, unit_main=7, unit_gene=2, unit_overview=1.25,unit_ratios=NULL){
  if(! is.null(unit_ratios)){
    tmp <- unlist(stringr::str_split(unit_ratios, ":"))
    if(length(tmp)==3){
      unit_overview=as.numeric(tmp[1])
      unit_main=as.numeric(tmp[2])
      unit_gene=as.numeric(tmp[3])
    }else if(length(tmp)==2 & !show_overview){ #locuszoomplot
      unit_main=as.numeric(tmp[1])
      unit_gene=as.numeric(tmp[2])
    }
  }
  g2 <- ggplotGrob(plots$main_plot)
  g3 <- ggplotGrob(plots$gene_plot)
  fg2 <- egg::gtable_frame(g2, height=unit(unit_main, "null"),debug=FALSE)
  fg3 <- egg::gtable_frame(g3, height=unit(unit_gene, "null"),debug=FALSE)

  if(show_overview){
    g1 <- ggplotGrob(plots$overview_plot)
    fg1 <- egg::gtable_frame(g1, height=unit(unit_overview, "null"), debug=FALSE)
    fg12 <- egg::gtable_frame(gridExtra::gtable_rbind(fg1,fg2,fg3),debug=FALSE)
  }else{
    fg12 <- egg::gtable_frame(gridExtra::gtable_rbind(fg2,fg3),debug=FALSE)
  }
  grid::grid.newpage()
  grid::grid.draw(fg12)
}
