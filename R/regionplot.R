

#' Region plot
#'
#' @description
#'
#' @description
#'
#' \code{regionplot()} displays the association results for a smaler region within one chromosome
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase), and either gene, variant or region (represented by chr, xmin and xmax)
#'
#' All other input parameters are optional
#'
#' @param show_overview Show the overview plot (default show_overview=TRUE)
#' @param gene Zoom in on this gene (e.g. gene=FTO)
#' @param show_genes Show genes instead of exons (default show_genes=FALSE)
#' @param show_exons Show exons instead of genees (default show_exons=FALSE)
#' @param vline Parameter (optional) to add a vertical line to the plot at a specific chromosomal position, e.g \code{vline=204000066}. Multiple values can be provided in a vector, e.g  \code{vline=c(204000066,100500188)}
#' @param rsids A vector of rs ids to highlight on the plot, e.g. rsids=c("rs1234, rs45898")
#' @param rsids_color The color of the variants in variants_id (default color is red)
#' @param max_genes Only label the genes if they are fewer than max_genes (default values is 200).
#' @param show_gene_names Show the gene names even though they exceed the max_genes count
#' @param variant Zoom in on this variant. Can be either an rsid, or a dataframe (like returned from get_best_hit())
#' @param gene Zoom in on this gene (e.g. gene=FTO)
#' @param region genetic region, e.g. chr1:67038906-67359979
#' @param gene_padding The size of the region around the gene, if the gene argument was used (default: gene_padding=100000)
#' @inheritParams manhattan
#'
#' @return plots using egg (https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html)
#' @export
#'
#' @examples
#' \dontrun{
#' regionplot(df,gene="FTO")
#' }


regionplot <- function(df, ntop=3, annotate=NULL, xmin=0, size=2, shape=19, alpha=1,label_size=4, annotate_with="ID",
                    color=get_topr_colors(), axis_text_size=11,axis_title_size=12,title_text_size=13, show_genes=FALSE, show_overview=TRUE,
                       show_exons=FALSE, max_genes=200, sign_thresh=5e-09, sign_thresh_color="red", sign_thresh_label_size=3.5,
                       xmax=NULL,ymin=NULL,ymax=NULL,protein_coding_only=FALSE,region_size=1000000,gene_padding=100000,angle=0,legend_title_size=12,legend_text_size=12,
                       nudge_x=0.01,nudge_y=0.01, rsids=NULL, variant=NULL,rsids_color=NULL,legend_name="Data:",legend_position="right",
                       chr=NULL,vline=NULL,show_gene_names=NULL,legend_labels=NULL,gene=NULL, title=NULL, label_color=NULL,locuszoomplot=F,region=NULL){
  # three plots, overview_plot, main_plot and gene_plot
  #only include overview plot if df region is larger than the region between xmin and xmax
  if(is.null(gene) & is.null(variant) & is.null(region) & (is.null(chr)  & is.null(xmax))){
    stop("Regional argument is missing (either gene, variant, region or the 3 arguments (chr,xmin,xmax) has to be provided as input to this function")
  }
  else{
  dat <- dat_check(df) %>% set_log10p(ntop)  %>% set_size_shape_alpha(size, shape, alpha, locuszoomplot = locuszoomplot)
    if(! locuszoomplot){
      dat <- dat %>% set_color(color)
    }
  top_snps <- NULL
  using_ntop <- FALSE
  if(length(dat) > ntop){
    using_ntop <- TRUE
  }
  if(! is.null(gene)){
    gene_df <- get_gene(gene)
    if(dim(gene_df)[1]==0){  stop(paste("Could not find gene ",gene)) }
    else{
      chr <- gene_df$chrom
      xmin <- gene_df$gene_start-gene_padding
      xmax <- gene_df$gene_end+gene_padding
    }
  }else if(! is.null(variant)){
    if(is.data.frame(variant)){
      if("POS" %in% colnames(variant)){
        chr <- variant$CHROM
        xmin <- variant$POS-(region_size/4)
        xmax <- variant$POS+(region_size/4)
      }
      else{
        stop("There is no POS column in the variant dataframe")
      }
    }
    else{
      tmp <- dat[[1]] %>% filter(ID == variant) %>% arrange(P) %>% head(n=1)
      if(length(tmp$POS > 0 )){
        chr <- tmp$CHROM
        xmin <- tmp$POS-(region_size/4)
        xmax <- tmp$POS+(region_size/4)
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
    print(paste("CHR ",chr, " xmin: ",xmin, " xmax: ",xmax, sep=""))
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
  if(! is.null(annotate)){ top_snps <- get_annotation(dat, annotate=annotate, region_size=region_size,protein_coding_only = protein_coding_only) }
  # get the main plot
  main_plot <- get_base_plot(dat,color=color,legend_labels=legend_labels, legend_name = legend_name,legend_position = legend_position,locuszoomplot=locuszoomplot) %>% set_plot_text_sizes(axis_text_size=axis_text_size,axis_title_size = axis_title_size,legend_text_size=legend_text_size, legend_title_size=legend_title_size)

  main_plot <- main_plot+theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())
    if(! is.null(annotate)){
      main_plot <- main_plot %>%  add_annotation(plot_labels = top_snps,annotate_with=annotate_with,angle=angle,label_size = label_size, nudge_y = nudge_y, nudge_x = nudge_x, label_color=label_color)
    }
    if(locuszoomplot){
      #annotate the top variant
      if(is.null(rsids)){
      top_snps <- get_main_LD_snp(dat)
      if(!is.null(top_snps) & length(top_snps$POS)> 0){
        main_plot <- main_plot %>%  add_annotation(plot_labels = top_snps,annotate_with=annotate_with,angle=angle,label_size = label_size, nudge_y = nudge_y, nudge_x = nudge_x, label_color=label_color)
      }
      }
    }
  main_plot <- main_plot %>% add_vline(vline)
  main_plot <- main_plot + scale_y_continuous(expand=c(.02,.02)) + scale_x_continuous(expand=c(.01,.01))

if(!is.null(sign_thresh)){
 main_plot <- main_plot %>% add_sign_thresh(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, using_ntop = using_ntop) %>%
    add_sign_thresh_labels(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, xmin=xmin,sign_thresh_label_size = sign_thresh_label_size)
}

  if( ! is.null(rsids)){
   main_plot <-main_plot %>% add_rsids(dat, rsids, rsids_color=rsids_color, nudge_x=nudge_x, nudge_y=nudge_y, label_size=label_size, angle=angle)
  }

  if(is.null(ymin)){ ymin <- get_ymin(dat)* 1.04}
  if(is.null(ymax)){ ymax <- get_ymax(dat)* 1.04}

  main_plot <- main_plot %>% set_xmin_xmax_ymin_ymax(xmin,xmax,ymin,ymax)

  #get the overview plot
  overview_plot <- get_base_plot(df_full,color=color,show_legend = FALSE,overview_plot = TRUE)
  overview_plot <- rm_axis_text_and_ticks(overview_plot)
  overview_plot <- add_zoom_rectangle(overview_plot,df_full,xmin, xmax)
  overview_plot <- overview_plot %>% add_sign_thresh(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, size=0.2, using_ntop = using_ntop)
  if(! is.null(title)){
   if(show_overview){
    #add the title to the overview plot
     overview_plot <- overview_plot %>% add_title(title=title, title_text_size = title_text_size)
   }else{ # add the title to the mainplot
     main_plot <- main_plot + add_title(title=title, title_text_size = title_text_size)
   }
   }
  if(using_ntop){
    main_plot <- main_plot %>%  add_zero_hline()
    overview_plot <- overview_plot %>% add_zero_hline()
  }


 if((xmax-xmin < 1000001 || show_exons) & !show_genes){
    show_exons <- TRUE
    show_genes <- FALSE
  }else{
    show_genes <- TRUE
   show_exons <- FALSE
  }
  ngenes <- 0

  #get the genes in the region
    genes <- get_genes_in_region(chr=chr,xmin=xmin,xmax=xmax,protein_coding_only=protein_coding_only, show_exons=show_exons,show_genes=show_genes)
    ngenes <- length(genes$gene_symbol)
    if(ngenes>100){show_gene_names <- F}
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

  print(paste("Zoomed to region:  ",chr,":",xmin,"-",xmax,sep=""))
  gene_plot <- get_gene_plot(plot_args,label_size=label_size,xmin = xmin, xmax = xmax, show_gene_names = show_gene_names)
  gene_plot <- gene_plot + scale_y_continuous(expand=c(.02,.02)) + scale_x_continuous(expand=c(.01,.01))
  gene_plot <- gene_plot + xlab(paste("Position on Chromosome ",gsub('chr','',chr), sep=""))
  gene_plot <- gene_plot %>% add_vline(vline)
  gene_plot <- gene_plot %>% set_legend_texts(legend_title_size=legend_title_size, legend_text_size=legend_text_size)
  gene_plot <- gene_plot  + theme(axis.text=element_text(size=axis_text_size), axis.title = element_text(size=axis_title_size))
  plots <- list(main_plot=main_plot, overview_plot=overview_plot, gene_plot=gene_plot)
  plots %>% draw_plots(show_overview = show_overview)
  }
}


draw_plots <- function(plots, show_overview=TRUE){
  g2 <- ggplotGrob(plots$main_plot)
  g3 <- ggplotGrob(plots$gene_plot)
  fg2 <- egg::gtable_frame(g2, height=unit(7, "null"),debug=F)
  fg3 <- egg::gtable_frame(g3, height=unit(2, "null"),debug=F)

  if(show_overview){
    g1 <- ggplotGrob(plots$overview_plot)
    fg1 <- egg::gtable_frame(g1, height=unit(1.25, "null"), debug=F)
    fg12 <- egg::gtable_frame(gridExtra::gtable_rbind(fg1,fg2,fg3),debug=FALSE)
  }else{
    fg12 <- egg::gtable_frame(gridExtra::gtable_rbind(fg2,fg3),debug=FALSE)
  }
  grid::grid.newpage()
  grid::grid.draw(fg12)
}
