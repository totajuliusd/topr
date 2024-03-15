#' Create a Manhattan plot
#'
#' @description
#'
#' \code{manhattan()} displays association results for the entire genome on a Manhattan plot.
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase)
#'
#' All other input parameters are optional
#'
#'
#' @param df Dataframe or a list of dataframes (required columns are \code{CHROM,POS,P}), in upper- or lowercase) of association results.
#' @param ntop An integer, number of datasets (GWAS results) to show on the top plot
#' @param chr A string or integer, the chromosome to plot (i.e. chr15), only required if the input dataframe contains results from more than one chromosome
#' @param title A string to set the plot title
#' @param color A string or a vector of strings, for setting the color of the datapoints on the plot
#' @param size A number or a vector of numbers, setting the size of the plot points (default: \code{size=1.2})
#' @param alpha A number or a vector of numbers setting the transparency of the plotted points
#' @param shape A number of a vector of numbers setting the shape of the plotted points
#' @param annotate A number (p-value). Display annotation for variants with p-values below this threshold
#' @param annotate_with A string. Annotate the variants with either Gene_Symbol or ID (default: "Gene_Symbol")
#' @param label_size An number to set the size of the plot labels (default: \code{label_size=3})
#' @param label_color A string or a vector of strings. To change the color of the gene or variant labels
#' @param sign_thresh A number or vector of numbers, setting the horizontal significance threshold (default: \code{sign_thresh=5e-8}). Set to NULL to hide the significance threshold.
#' @param sign_thresh_color A string or vector of strings to set the color/s of the significance threshold/s
#' @param highlight_genes A string or vector of strings, gene or genes to highlight at the bottom of the plot
#' @param highlight_genes_ypos An integer, controlling where on the y-axis the highlighted genes are placed (default value is 1)
#' @param highlight_genes_color A string, color for the highlighted genes (default: darkred)
#' @param xmin,xmax Integer, setting the chromosomal range to display on the x-axis
#' @param ymin,ymax Integer, min and max of the y-axis, (default values: \code{ymin=0, ymax=max(-log10(df$P))})
#' @param legend_labels A string or vector of strings representing legend labels for the input datasets
#' @param legend_name A string, use to change the name of the legend (default: None)
#' @param legend_position A string, top,bottom,left or right
#' @param title_text_size A number, size of the plot title (default: 13)
#' @param axis_text_size A number, size of the x and y axes tick labels (default: 12)
#' @param axis_title_size A number, size of the x and y title labels (default: 12)
#' @param legend_title_size A number, size of the legend title
#' @param legend_text_size A number, size of the legend text
#' @param legend_nrow An integer, sets the number of rows allowed for the legend labels
#' @param protein_coding_only A logical scalar, if TRUE, only protein coding genes are used for annotation
#' @param region_size An integer (default = 20000000) (or a string represented as 200kb or 2MB) indicating the window size for variant labeling. Increase this number for sparser annotation and decrease for denser annotation.
#' @param nudge_x  A number to vertically adjust the starting position of each gene label (this is a ggrepel parameter)
#' @param nudge_y  A number to horizontally adjust the starting position of each gene label (this is a ggrepel parameter)
#' @param angle A number, the angle of the text label
#' @param sign_thresh_label_size  A number setting the text size of the label for the significance thresholds (default text size is 3.5)
#' @param gene_label_size A number setting the size of the gene labels shown at the bottom of the plot
#' @param gene_label_angle A number setting the angle of the gene label shown at the bottom of the plot (default: 0)
#' @param scale A number, to change the size of the title and axes labels and ticks at the same time (default : 1)
#' @param show_legend A logical scalar, set to FALSE to hide the legend (default : TRUE)
#' @param sign_thresh_linetype A string, the line-type of the horizontal significance threshold (default : dashed)
#' @param sign_thresh_size A number, sets the size of the horizontal significance threshold line (default : 1)
#' @param rsids A string (rsid) or vector of strings to highlight on the plot, e.g. \code{rsids=c("rs1234, rs45898")}
#' @param rsids_color A string, the color of the variants in variants_id (default color is red)
#' @param rsids_with_vline  A string (rsid) or vector of strings to highlight on the plot with their rsids and vertical lines further highlighting their positions
#' @param annotate_with_vline A number (p-value). Display annotation and vertical lines for variants with p-values below this threshold
#' @param shades_color The color of the rectangles (shades) representing the different chromosomes on the Manhattan plot
#' @param shades_alpha The transparency (alpha) of the rectangles (shades)
#' @param segment.size 	line segment color (ggrepel argument)
#' @param segment.color line segment thickness (ggrepel argument)
#' @param segment.linetype line segment solid, dashed, etc.(ggrepel argument)
#' @param max.overlaps Exclude text labels that overlap too many things. Defaults to 10 (ggrepel argument)
#' @param label_fontface  A string or a vector of strings. Label font “plain”, “bold”, “italic”, “bold.italic” (ggrepel argument)
#' @param label_family A string or a vector of strings. Label font name (default ggrepel argument is "")
#' @param gene_label_fontface  Gene label font “plain”, “bold”, “italic”, “bold.italic” (ggrepel argument)
#' @param gene_label_family Gene label font name (default ggrepel argument is "")
#' @param build A number representing the genome build or a data frame. Set to 37 to change to build (GRCh37). The default is build 38 (GRCh38).
#' @param verbose A logical scalar (default: NULL). Set to FALSE to suppress printed messages 
#' @param label_alpha An number or vector of numbers to set the transparency of the plot labels (default: \code{label_alpha=1})
#' @param shades_line_alpha The transparency (alpha) of the lines around the rectangles (shades)
#' @param vline A number or vector of numbers to add a vertical line to the plot at a specific chromosomal position, e.g \code{vline="chr1:204000066"}. Multiple values can be provided in a vector, e.g  \code{vline=c("chr1:204000066","chr5:100500188")}
#' @param vline_color A string. The color of added vertical line/s (default: grey)
#' @param vline_linetype A string. The linetype of added vertical line/s (default : dashed)  
#' @param vline_alpha A number. The alpha of added vertical line/s (default : 1)  
#' @param vline_size A number.The size of added vertical line/s (default : 0.5)  
#' @param region A string representing a genetic region, e.g. chr1:67038906-67359979
#' @param theme_grey A logical scalar (default: FALSE). Use gray rectangles (instead of white to distinguish between chromosomes)
#' @param xaxis_label A string. The label for the x-axis (default: Chromosome)
#' @param use_shades A logical scalar (default: FALSE). Use shades/rectangles to distinguish between chromosomes
#' @param even_no_chr_lightness Lightness value for even numbered chromosomes. A number or vector of numbers between 0 and 1 (default: 0.8). If set to 0.5, the same color as shown for odd numbered chromosomes is displayed. A value below 0.5 will result in a darker color displayed for even numbered chromosomes, whereas a value above 0.5 results in a lighter color.
#' @param get_chr_lengths_from_data A logical scalar (default: TRUE). If set to FALSE, use the inbuilt chromosome lengths (from hg38), instead of chromosome lengths based on the max position for each chromosome in the input dataset/s.
#' @param log_trans_p A logical scalar (default: TRUE). By default the p-values in the input datasets are log transformed using -log10. Set this argument to FALSE if the p-values in the datasets have already been log transformed. 
#' @param chr_ticknames A vector containing the chromosome names displayed on the x-axis. If NULL, the following format is used: chr_ticknames <- c(1:16, '',18, '',20, '',22, 'X')
#' @param show_all_chrticks A logical scalar (default : FALSE). Set to TRUE to show all the chromosome names on the ticks on the x-axis
#' @param hide_chrticks_from_pos A number (default: 17). Hide every nth chromosome name on the x-axis FROM this position (chromosome number)
#' @param hide_chrticks_to_pos A number (default: NULL). Hide every nth chromosome name on the x-axis TO this position (chromosome number). When NULL this variable will be set to the number of numeric chromosomes in the input dataset.
#' @param hide_every_nth_chrtick A number (default: 2). Hide every nth chromosome tick on the x-axis (from the hide_chr_ticks_from_pos to the hide_chr_ticks_to_pos).
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import dplyr
#' @import utils
#'
#' @examples
#' \dontrun{
#' manhattan(CD_UKBB)
#' }

manhattan <- function(df, ntop=4, title="",annotate=NULL, color=get_topr_colors(),
                   sign_thresh=5e-08,sign_thresh_color="red", sign_thresh_label_size=3.5, label_size=3.5, size=0.8,shape=19,alpha=1,
                   highlight_genes_color="darkred",highlight_genes_ypos=1.5,axis_text_size=12,axis_title_size=14, title_text_size=15,
                   legend_title_size=13,legend_text_size=12, protein_coding_only=TRUE,angle=0,
                   legend_labels=NULL,chr=NULL, annotate_with="Gene_Symbol",region_size=20000000,legend_name=NULL,
                   legend_position="bottom", nudge_x=0.1,nudge_y=0.7,xmin=NULL, xmax=NULL,ymin=NULL,ymax=NULL,
                   highlight_genes=NULL,label_color=NULL,legend_nrow=NULL,gene_label_size=NULL,gene_label_angle=0,
                   scale=1,show_legend=TRUE,sign_thresh_linetype="dashed", sign_thresh_size=0.5,rsids=NULL, rsids_color=NULL,
                  rsids_with_vline=NULL,annotate_with_vline=NULL,shades_color=NULL,shades_alpha=0.5,segment.size=0.2, 
                  segment.color="black",segment.linetype="dashed",max.overlaps=10,label_fontface="plain",label_family="",
                  gene_label_fontface="plain",gene_label_family="",build=38,verbose=NULL,label_alpha=1,shades_line_alpha=1,vline=NULL,
                  vline_color="grey",vline_linetype="dashed", vline_alpha=1,vline_size=0.5,region=NULL, theme_grey=FALSE, xaxis_label="Chromosome",
                  use_shades=FALSE, even_no_chr_lightness=0.8, get_chr_lengths_from_data=TRUE, log_trans_p=TRUE,
                  chr_ticknames=NULL, show_all_chrticks=FALSE, hide_chrticks_from_pos=17, hide_chrticks_to_pos=NULL, hide_every_nth_chrtick=2){
    
    top_snps <- NULL
    genes_df <- NULL
    chr_map <- NULL
   
    if(theme_grey)
      use_shades=T
    dat <- dat_check(df, verbose=verbose, log_trans_p) 
    
    if(! is.null(chr)){
      dat <- dat %>% filter_on_chr(chr)
      xaxis_label <- paste(xaxis_label, gsub("chr", "", chr), sep=" ")
      if(! is.null(xmin) & ! is.null(xmax))
        dat <- dat %>% filter_on_xmin_xmax(xmin,xmax)
    }
    else{
      datl <- dat %>% convert_chrs_to_numeric(get_chr_lengths_from_data) #convert chromosomes to numeric so they can be numerically sorted
      dat <- datl$dat; chr_map <- datl$chr_map;
    }
    dat <- dat %>% set_size_shape_alpha(size, shape, alpha) %>% set_color(color,shades_alpha,use_shades,even_no_chr_lightness,chr) 
    if(log_trans_p) 
      dat <- dat %>% set_log10p(ntop)
    else{
      warning("Assuming p-values have already been log transformed since [log_trans_p] is set to FALSE and plotting the data as is (without log transforming the p-values)!  ")
      dat <- dat %>% add_log10p_wo_trans(ntop)
    }
    if(!use_shades & ! is.null(shades_color)){
      warning(paste0("Argument use_shades is set to FALSE by default. For the shades_color argument to have an effect, the use_shades argument has to be set to TRUE. Add the argument [use_shades=TRUE] and re-run."))
    }
   
    using_ntop <- FALSE
    if(length(unique(dat[[1]]$CHROM))==1){chr<-unique(dat[[1]]$CHROM)}
    
    if(! is.null(region)){
      tmp <- unlist(stringr::str_split(region, ":"))
      chr <- tmp[1]
      tmp_pos <- unlist(stringr::str_split(tmp[2], "-"))
      xmin <- as.numeric(tmp_pos[1])
      xmax <- as.numeric(tmp_pos[2])
    }
  if(length(dat) > ntop){
      using_ntop <- TRUE
  }
  annot_with_vline <- FALSE
  if(! is.null(annotate_with_vline)){
    annotate <- annotate_with_vline
    annot_with_vline <- TRUE
  }
   
  
  if(is.null(ymin)){
    ymin <- get_ymin(dat)
    if(!is.null(highlight_genes)){ 
      ymin <- ifelse(highlight_genes_ypos < ymin, highlight_genes_ypos, ymin) }
  }
  if(is.null(ymax)){ 
    ymax <- get_ymax(dat)
    ymax_addon <- ymax *0.04
    if(ymax_addon < 0.7) # to make sure there will be space for annotation
      ymax_addon <- 0.7
    ymax <- ymax+ymax_addon
  }

    # get the annotation
    if(! is.null(annotate)){
      top_snps <- get_annotation(dat, region_size = region_size, annotate=annotate, protein_coding_only = protein_coding_only,nudge_x=nudge_x,nudge_y=nudge_y,
                                 angle=angle,label_fontface=label_fontface,label_family=label_family, build=build, verbose = verbose, label_alpha=label_alpha, chr_map=chr_map)

    }
      #get the genes
    if (! is.null(highlight_genes)){
      if(is.data.frame(highlight_genes)){
        genes_df <- highlight_genes
      }else{
        genes_df <- get_genes_by_Gene_Symbol(highlight_genes,chr, build=build)
      }
    }
    offsets=NULL
    incl_chrX=T
    if(length(unique(dat[[1]]$CHROM))>1 & is.null(chr)){  #Manhattan plot
      incl_chrX <- include_chrX(dat)
      chr_lengths_and_offsets <- get_chr_lengths_and_offsets(dat, get_chr_lengths_from_data)
      offsets <- stats::setNames(chr_lengths_and_offsets$offset,chr_lengths_and_offsets$CHROM)
      if(! is.null(annotate)){  top_snps <- top_snps %>%  get_pos_with_offset(offsets) }
      if (! is.null(highlight_genes)){
        genes_df$CHROM <- gsub("chr", "", genes_df$CHROM)
        genes_df <- genes_df  %>% get_pos_with_offset(offsets)
      }
      dat <- get_pos_with_offset4list(dat,offsets)
    }
    
    main_plot <- get_base_plot(dat,color=color,legend_labels = legend_labels,legend_name=legend_name, legend_position = legend_position, 
                               legend_nrow = legend_nrow, show_legend = show_legend,scale=scale,verbose=verbose)
  
        if(! is.null(title)){
        main_plot <- main_plot %>% add_title(title=title, title_text_size = title_text_size,scale=scale)
     }
    if(is.null(chr)){

      ticks <- get_ticks(dat,chr_lengths_and_offsets,chr_ticknames,chr_map,show_all_chrticks, hide_chrticks_from_pos, hide_chrticks_to_pos, hide_every_nth_chrtick,get_chr_lengths_from_data)
      if(use_shades)
        shades <- get_shades(chr_lengths_and_offsets,dat,ntop=ntop,include_chrX = incl_chrX,ymin=ymin,ymax=ymax)
      
     main_plot <- main_plot %>% add_shades_and_ticks(shades,ticks,shades_color=shades_color,shades_alpha=shades_alpha,shades_line_alpha=shades_line_alpha, theme_grey=theme_grey, use_shades=use_shades)
    }else{
      main_plot <- main_plot + scale_y_continuous(expand=c(.02,.02))  + scale_x_continuous(expand=c(.01,.01),labels = scales::comma)
    }
 
  main_plot <- set_axis_labels(main_plot,xaxis_label = xaxis_label)
  main_plot <- main_plot %>% set_plot_text_sizes(axis_text_size=axis_text_size,axis_title_size = axis_title_size, 
                                                 legend_text_size=legend_text_size, legend_title_size=legend_title_size,scale=scale)


   #add the significance threshold/s
  if(!is.null(sign_thresh)){
  main_plot <- main_plot %>% add_sign_thresh(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, using_ntop = using_ntop, 
                                             sign_thresh_linetype = sign_thresh_linetype, sign_thresh_size = sign_thresh_size,scale=scale) %>%
   add_sign_thresh_labels(sign_thresh = sign_thresh, sign_thresh_color = sign_thresh_color, xmin=xmin, sign_thresh_label_size = sign_thresh_label_size,scale=scale)
  }
     if(using_ntop){
    main_plot <- main_plot %>%  add_zero_hline()
}
  if (! is.null(annotate)){
    main_plot <- main_plot %>%  add_annotation(plot_labels = top_snps,annotate_with=annotate_with,angle=angle,label_size = label_size, 
                                               label_color=label_color, nudge_x=nudge_x, nudge_y=nudge_y, scale=scale, 
                                               segment.size=segment.size,segment.color=segment.color,segment.linetype=segment.linetype, max.overlaps=max.overlaps)
  }

  
  if (! is.null(highlight_genes)){
    if(! is.null(gene_label_size)){
      label_size <- gene_label_size
    }
    main_plot <- add_genes2plot(main_plot, genes_df, highlight_genes_ypos=highlight_genes_ypos,highlight_genes_color=highlight_genes_color, 
                                label_size=label_size,gene_label_angle = gene_label_angle,scale=scale,gene_label_fontface=gene_label_fontface,gene_label_family=gene_label_family)
    
  }

  if(annot_with_vline){
    main_plot <- main_plot %>% add_vline(top_snps$POS, vline_color=vline_color, vline_linetype = vline_linetype, vline_alpha=vline_alpha, vline_size=vline_size,scale=scale)
  }
  with_vline <- FALSE
  if(! is.null(rsids_with_vline)){
    rsids <- rsids_with_vline
    with_vline <- TRUE
  }
  if(! is.null(rsids)){
    rsids_df <- get_rsids_from_df(dat,rsids)
    main_plot <-main_plot %>% add_rsids(rsids_df, rsids_color=rsids_color, nudge_x=nudge_x, nudge_y=nudge_y, label_size=label_size, angle=angle, label_color=label_color, scale=scale, with_vline = with_vline)
  }

  if(! is.null(vline)){
    v <- data.frame("x"=vline)
    v <- v %>% tidyr::separate("x", c("CHROM","POS"),":")
    v <- v %>% dplyr::mutate(CHROM=gsub('chr','',CHROM,ignore.case = T))
    v$CHROM <- as.numeric(v$CHROM)
    v$POS <- as.numeric(v$POS)
    if(is.null(offsets)){
      vlines <- v$POS
      }
    else{
      vlines_w_offsets <- v %>% get_pos_with_offset(offsets)
      vlines <- as.vector(vlines_w_offsets$POS)
    }
    main_plot <- main_plot %>% add_vline(vlines, vline_color=vline_color, vline_linetype = vline_linetype, vline_alpha=vline_alpha, vline_size=vline_size,scale=scale)
  }
   if(!is.null(ymax) & !is.null(ymin)){
     main_plot <- main_plot %>% set_ymin_ymax(ymin,ymax)
   }
   if(!is.null(xmin) && !is.null(xmax) & ! is.null(chr)){
    main_plot <- main_plot %>% set_xmin_xmax(xmin,xmax)
  }
  
  main_plot <- change_axes(main_plot)
  return(main_plot)
  }

