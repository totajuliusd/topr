#' Manhattan plot
#'
#' @description
#'
#' \code{manhattan()} displays the association results for the entire genome
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase)
#'
#' All other input parameters are optional
#'
#'
#' @param df Dataframe or a list of dataframes (required columns are \code{CHROM,POS,P}), in upper- or lowercase) of association results.
#' @param ntop Number of datasets (GWASes) to show on the top plot
#' @param chr The chromosome to plot (i.e. chr15), only required if the input dataframe contains results from more than one chromosome
#' @param title Plot title (optional
#' @param color Optional parameter setting the color of the plot points (default: \code{color="darkblue"})
#' @param size Optional parameter setting the size of the plot points (default: \code{size=1.2})
#' @param alpha A number or vector of numbers setting the transparancy of the plotted points
#' @param shape A number of vector of numers setting the shape of the plotted points
#' @param annotate Display annotation for variants with p-values below this threshold
#' @param annotate_with Annotate the variants with eiher Gene_Symbol or ID (default annotate_with="Gene_Symbol")
#' @param label_size Optional parameter to set the size of the plot labels (default: \code{label_size=3.5})
#' @param sign_thresh Optional parameter setting the threshold of the dashed red horizontal line representing the significance threshold (default: \code{sign_thresh=5.1e-9}). Multiple thresholds can be provided in a vector, e.g \code{sign_thresh=c(5.1e-9,1.0e-6)}). Set this parameter to NULL if you dont want this line to appear at all \code{sign_thresh=NULL}
#' @param sign_thresh_color set the color of the significance threshold line or lines
#' @param highlight_genes A vector of genes or genes to highlight the datapoints for on the plot
#' @param highlight_genes_ypos Display the genes at this position on the y-axis (default value is 1)
#' @param highlight_genes_color Colors for the hihglighted genes (default: green)
#' @param xmin,xmax Parameters setting the chromosomal range to display on the x-axis
#' @param ymin,ymax Optional parameters, min and max of the y-axis, (default values: \code{ymin=0, ymax=max(-log10(df$P))})
#' @param rect Rectangle to add to the plot
#' @param legend_labels Legend labels
#' @param legend_name Change the name of the legend (default: None)
#' @param legend_position Top,bottom,left or right
#' @param title_text_size Text size of the plot title (default: 13)
#' @param axis_text_size Text size of the x and y axes tick labels (default: 12)
#' @param axis_title_size Text size of the x and y title labels (default: 12)
#' @param legend_title_size Text size of the legend title
#' @param legend_text_size Text size of the legend text
#' @param protein_coding_only Set this parameter to TRUE to only use protein coding genes for annotation
#' @param region the size of the region used when annotating the top variant in a region (default value is 10000000 or 10 MB)
#' @param nudge_x  To vertically adjust the starting position of each gene label (this is a ggrepel parameter)
#' @param nudge_y  To horizontally adjust the starting position of each gene label (this is a ggrepel parameter)
#' @param angle The angle of the text label
#' @param gene Zoom in on this gene (e.g. gene=FTO)
#' @param gene_region The size of the region around the gene, if the gene argument was used (default: gene_region=100000)
#' @param sign_thresh_label_size  Set the text size of the label for the signficance thresholdds (default text sizze is 3.5)
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import dplyr
#' @import utils
#'
#' @examples
#' \dontrun{
#' data(gwas_CD)
#' manhattan(gwas_CD)
#' }
#'

manhattan <- function(df, ntop=3, title="",annotate=NULL, color=get_topr_colors(),
                   sign_thresh=5e-09,sign_thresh_color="red", sign_thresh_label_size=3.5, label_size=3, size=1,shape=19,alpha=1,highlight_genes_color="green",highlight_genes_ypos=1,
                   axis_text_size=11,axis_title_size=12, title_text_size=13,legend_title_size=12,legend_text_size=12, protein_coding_only=TRUE,angle=0,
                   legend_labels=NULL,gene=NULL,gene_region=100000, chr=NULL, annotate_with="Gene_Symbol",region=10000000,
                      legend_name=NULL,legend_position=NULL,rect=NULL, nudge_x=0.1,nudge_y=0.2,
                      xmin=NULL, xmax=NULL,ymin=NULL,ymax=NULL,highlight_genes=NULL){
    top_snps <- NULL
    genes_df <- NULL
    xaxis_label <- "Chromosome"
    dat <- dat_check(df) %>% set_size_shape_alpha(size, shape, alpha) %>% set_color(color) %>% set_log10p(ntop)
    using_ntop <- FALSE
    if(length(dat) > ntop){
      using_ntop <- TRUE
    }
    if(! is.null(gene)){
      gene_df <- get_gene(gene,chr)
      if(dim(gene_df)[1]==0){  stop(paste("Could not find gene ",gene)) }
      else{
        chr <- gene_df$chrom
        xmin <- gene_df$gene_start-gene_region
        xmax <- gene_df$gene_end+gene_region
      }
    }
    if(! is.null(chr)){  dat <- dat %>% filter_on_chr(chr)
      xaxis_label <- paste(xaxis_label, gsub("chr", "", chr), sep=" ")
      if(! is.null(xmin) & ! is.null(xmax)){
        dat <- dat %>% filter_on_xmin_xmax(xmin,xmax)
      }
    }
  if(is.null(ymin)){
    ymin <- get_ymin(dat)
    if(!is.null(highlight_genes)){ ymin <- ifelse(highlight_genes_ypos < ymin, highlight_genes_ypos, ymin) }
  }
  if(is.null(ymax)){ ymax <- get_ymax(dat) *1.04}
    # get the annotation
    if(! is.null(annotate)){ top_snps <- get_annotation(dat, region=region, annotate=annotate, protein_coding_only = protein_coding_only) }
    #get the genes
    if (! is.null(highlight_genes)){
      if(is.data.frame(highlight_genes)){
        genes_df <- highlight_genes
      }else{
        genes_df <- get_genes_by_Gene_Symbol(highlight_genes,chr)
      }
    }
    if(length(unique(dat[[1]]$CHROM))>1 & is.null(chr)){  #Manhattan plot
      incl_chrX <- include_chrX(dat)
      offsets <- get_chr_offsets(incl_chrX)
      shades <- get_shades(offsets,dat,ntop=ntop,include_chrX = incl_chrX,ymin=ymin,ymax=ymax)
      if(! is.null(annotate)){  top_snps <- top_snps %>%  get_pos_with_offset(offsets) }
      if (! is.null(highlight_genes)){
        genes_df$CHROM <- gsub("chr", "", genes_df$CHROM)
        genes_df <- genes_df  %>% get_pos_with_offset(offsets)
      }
      dat <- get_pos_with_offset4list(dat,offsets)
    }
    # Do the plotting
    main_plot <- get_base_plot(dat,color=color,legend_labels = legend_labels,legend_name=legend_name,legend_position = legend_position)
      if(! is.null(title)){
        main_plot <- main_plot %>% add_title(title=title, title_text_size = title_text_size)
    }
    if(is.null(chr)){
      ticks <- get_ticks(dat)
      main_plot <- main_plot %>% add_shades_and_ticks(shades,ticks)
    #  main_plot <- main_plot + scale_y_continuous(expand=c(.02,.02))
    }else{
      main_plot <- main_plot + scale_y_continuous(expand=c(.02,.02))  + scale_x_continuous(expand=c(.01,.01),labels = scales::comma)
    }

  main_plot <- set_axis_labels(main_plot,xaxis_label = xaxis_label)
  main_plot <- main_plot %>% set_plot_text_sizes(axis_text_size=axis_text_size,axis_title_size = axis_title_size, legend_text_size=legend_text_size, legend_title_size=legend_title_size)

  #add the significance threshold/s
  main_plot <- main_plot %>% add_sign_thresh(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, using_ntop = using_ntop) %>%
   add_sign_thresh_labels(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, xmin=xmin, sign_thresh_label_size = sign_thresh_label_size)
  if(using_ntop){
    main_plot <- main_plot %>%  add_zero_hline()
}
  main_plot <- main_plot %>%  add_annotation(plot_labels = top_snps,annotate_with=annotate_with,angle=angle,label_size = label_size, nudge_x=nudge_x, nudge_y=nudge_y)

  if (! is.null(highlight_genes)){
    main_plot <- add_genes2plot(main_plot, genes_df, highlight_genes_ypos=highlight_genes_ypos,highlight_genes_color=highlight_genes_color)
  }
  main_plot <- main_plot %>% set_ymin_ymax(ymin,ymax)
  main_plot <- change_axes(main_plot)
  return(main_plot)
  }

