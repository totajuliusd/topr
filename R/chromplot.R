
#' Chromosome plot
#'
#' @description
#'
#' \code{chromplot()} displays the association results for an entire chromosome.
#' Required parameter is a dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase)
#'
#' All other input parameters are optional. A dataset including results from more than one
#' chromosome can be used as input if provided with the \code{chr parameter}.
#'
#' E.g. \code{chromplot(multi_chr_dataset,chr = "chr15")}
#'
#'
#' @param df Dataframe (required columns are \code{CHROM,POS,P}), in upper- or lowercase) to use for the plot. See examples on how to retrieve the dataset with GOR.
#' @param chr The chromosome to plot (i.e. chr15), only required if the input dataframe contains results from more than one chromosome
#' @param title Plot title (optional)
#' @param variants Optional parameter, a dataset (required columns are \code{CHROM,POS,P}) with a list of variants to display and label on the plot.
##' Variants will be labeled with Gene_Symbol if present in the dataframe, otherwise with ID.  If neither Gene_Symbol nor ID are in the dataset, the variant will be labelled with its position (CHROM and POS)
#' @param color Optional parameter setting the color of the plot points (default: \code{color="darkblue"})
#' @param size Optional parameter setting the size of the plot points (default: \code{size=1.2})
#' @param alpha A number or vector of numbers setting the transparancy of the plotted points
#' @param shape A number of vector of numers setting the shape of the plotted points
#' @param annotation_thresh Display annotation for variants with p-values below this threshold
#' @param label_size Optional parameter to set the size of the plot labels (default: \code{label_size=3.5})
#' @param sign_thresh Optional parameter setting the threshold of the dashed red horizontal line representing the significance threshold (default: \code{sign_thresh=5.1e-9}). Multiple thresholds can be provided in a vector, e.g \code{sign_thresh=c(5.1e-9,1.0e-6)}). Set this parameter to NULL if you dont want this line to appear at all \code{sign_thresh=NULL}
#' @param sign_thresh_color set the color of the significance threshold line or lines
#' @param variant_ids A vector of variant ids to highlight on the plot, e.g. variant_ids=c("rs1234, rs45898")
#' @param variant_ids_color The color of the variants in variants_id (default color is red)
#' @param highlight_genes A vector of genes or genes to highlight the datapoints for on the plot
#' @param highlight_genes_ypos Display the genes at this position on the y-axis (default value is 1)
#' @param highlight_genes_color Colors for the hihglighted genes (default: green)
#' @param xmin,xmax Parameters setting the chromosomal range to display on the x-axis
#' @param ymin,ymax Optional parameters, min and max of the y-axis, (default values: \code{ymin=0, ymax=max(-log10(df$P))})
#' @param vline Parameter (optional) to add a vertical line to the plot at a specific chromosomal position, e.g \code{vline=204000066}. Multiple values can be provided in a vector, e.g  \code{vline=c(204000066,100500188)}
#' @param rect Rectangle to add to the plot
#' @param rect_range Sets the size of the rectangle, default is 500 kb
#' @param variants_color One color or a vector of colors in case multiple variant dataframes are provided
#' @param variants_size A number or a vector of numbers (default: 1.2)
#' @param variants_alpha A number or a vector of numbers (default: 1)
#' @param variants_shape A number or a vector of numbers (default: 19)
#' @param legend_labels Legend labels
#' @param legend_name Change the name of the legend (default: None)
#' @param legend_position Top,bottom,left or right
#' @param title_text_size Text size of the plot title (default: 13)
#' @param axis_text_size Text size of the x and y axes tick labels (default: 12)
#' @param axis_title_size Text size of the x and y title labels (default: 12)
#' @param legend_title_size Text size of the legend title
#' @param legend_text_size Text size of the legend text
#' @param geneplot_label_size Size of the labels for the genes on the gene and exon plots
#' @param show_xaxis show the xaxis
#' @param protein_coding_only Set this parameter to TRUE to only use protein coding genes for annotation
#'
#' @return a chromosome plot (ggplot object)
#' @export
#'
#' @examples
#' \dontrun{
#' data(gwas_CD)
#' snps_CD <- get_best_snp_per_MB(gwas_CD, thresh = 1e-09, region = 10000000)
#' CHR="chr16"
#' chromplot(gwas_CD,chr=CHR, variants = snps_CD, annotation_thresh = 5e-09)
#' }


chromplot=function(dat, annotation_thresh = NULL, title = "", size = 1.2, shape = 19, alpha = 1,
                   color = c("darkblue","#E69F00","#00AFBB","#999999","#FC4E07","darkorange1"),
                   label_size=3.5,sign_thresh=5.1e-9, sign_thresh_color="red", legend_position="right",.checked = FALSE, show_xaxis = TRUE,
                   axis_text_size=11,title_text_size=12,axis_title_size=13,legend_title_size=12, legend_text_size = 12,
                   variant_ids_color=NULL,protein_coding_only=FALSE,
                   legend_labels = NULL, legend_name="Data",xmin=0, highlight_genes=NULL,highlight_genes_ypos=NULL,highlight_genes_color=NULL,
                   variants=NULL,xmax=NULL,ymin=NULL,ymax=NULL,chr=NULL,vline=NULL,variant_ids=NULL){



  if (! .checked) { # if input data has not already been checked
    #check and set data
    is_df_empty(dat, "main dat")

    if (is.data.frame(dat)) {
      dat <- list(dat)
    } else if (is.null(chr)) {
      chr <- get_chr_from_df(dat[[1]])
    }
    dat=dat_column_check_and_set(dat)
    dat <- dat_chr_check(dat) %>%
      filter_on_chr(chr) %>%
      filter_on_xmin_xmax(xmin, xmax) %>%
      set_size_shape_alpha(size, shape, alpha) %>%
      set_color(color)

    #check and set variants
    if(! is.null(annotation_thresh)){
      if (! is.null(variants)) {
        is_df_empty(variants, "variants")
        if(is.data.frame(variants)) variants=list(variants)
        variants=dat_column_check_and_set(variants)
        variants=dat_chr_check(variants)
      }
      else{
        variants=get_best_snp_per_MB(dat, thresh = annotation_thresh,protein_coding_only = protein_coding_only)
      }
      variants=filter_on_chr(variants,chr)
      variants=filter_on_xmin_xmax(variants,xmin,xmax)
      variant_color=color
      if(length(dat) == 1) variant_color="red"
      variants=set_size_shape_alpha(variants,size,shape,alpha)
      variants=set_color(variants,variant_color)
      variants=set_annotation_thresh(variants,annotation_thresh)
    }

  }
  df=dat[[1]]
  zoom_y=0
  if(is.null(ymax))
    ymax=max(-log10(df$P))*1.04
  else zoom_y=1

  if(is.null(ymin)) ymin=min(-log10(df$P))

  p1=ggplot(df)+geom_point(data=df, aes(x=POS, y=-log10(P),color=color),shape=df$shape, alpha=df$alpha, size=df$size) +theme_bw()+labs(shape="Max Impact")

  if(length(dat)>1){
    for(i in 2:length(dat)){
      df=dat[[i]]
      p1=p1+geom_point(data=df, aes(x=POS, y=-log10(P),color=color),shape=df$shape,size=df$size,alpha=df$alpha) # ,shape=df$shape, color=df$color, size=df$size)
    }
  }
  if(!is.null(vline)){
    for(i in 1:length(vline)){
      p1=p1+geom_vline(xintercept =vline[i] , colour="grey", linetype="dashed")
    }
  }

  if(!is.null(variants)){
    for(i in 1:length(variants)){
      p1=p1+geom_point(data=variants[[i]], aes(x=POS, y=-log10(P)),color=variants[[i]]$color, size=variants[[i]]$size, shape=variants[[i]]$shape)
    }
  }
  if(! is.null(annotation_thresh)){
    if(! is.null(variants)){
      for(i in 1:length(variants)){
        if(length(dat) == 1) {
          varcol="grey40"
        }else{
          varcol=variants[[i]]$color
        }
        if(show_xaxis) { #we are in whole chromosome view and want to display the Gene_Symbol for the variant
          p1=p1+suppressWarnings(ggrepel::geom_text_repel(data=variants[[i]], aes(x=POS,y=-log10(P),label=ifelse(P<annotation_thresh, Gene_Symbol,"")), color=varcol,size=label_size,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5))
        }else{ #in region view, display the variant ID instead of the gene_symbol
          p1=p1+suppressWarnings(ggrepel::geom_text_repel(data=variants[[i]], aes(x=POS,y=-log10(P),label=ifelse(P<annotation_thresh, ID,"")), color=varcol,size=label_size,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5))
        }
      }
    }
  }
  p1=color_genes(p1,dat, highlight_genes, highlight_genes_color,highlight_genes_ypos)
  p1=add_sign_thresh_to_plot(p1,df, sign_thresh, sign_thresh_color,xmin,xmax,ymin,ymax)

  if(! is.null(variant_ids)){
      for(i in 1:length(dat)){
        df=dat[[i]]
        variants2label=df %>% filter(ID %in% variant_ids) %>% distinct(ID, .keep_all = T)
        p1=p1+suppressWarnings(ggrepel::geom_text_repel(data=variants2label, aes(x=POS,y=-log10(P), label=ID), color="grey40",direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5))
        if(is.null(variant_ids_color)) variant_ids_color="red"
        p1=p1+geom_point(data=variants2label, aes(x=POS, y=-log10(P),color=variant_ids_color),shape=variants2label$shape,size=variants2label$size,alpha=variants2label$alpha)
      }
  }

  chr_label=gsub("chr", "", chr)

  p1=tidy_plot(p1, axis_text_size=axis_text_size,axis_title_size=axis_title_size, title_text_size=title_text_size,legend_title_size=legend_title_size,legend_text_size=legend_text_size)
  p1=p1+scale_x_continuous(expand=c(.01,.01),labels = scales::comma )

  p1=p1+labs(x=paste("Position on chr" ,chr_label, sep=""), y=expression(-log[10](italic(p))))

  if(show_xaxis){
    p1=p1+theme(panel.border=element_blank(),,axis.line=element_line(color="grey"))
  }else{
    p1=p1+theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())
  }



  if((! is.null(xmax) && (! is.null(xmin))))
    p1=p1 + coord_cartesian(xlim=c(xmin,xmax))

  if(zoom_y)  p1=p1+coord_cartesian(ylim=c(ymin,ymax))

    #add the legend
  if(! is.null(legend_name)) legend_name=legend_name
  legend_set2color=0
  if(is.null(legend_labels)) {
    legend_labels=color[1:length(dat)]
    legend_set2color=1
  }
  p1=p1+scale_color_identity(guide = "legend", name=legend_name,  breaks=color[1:length(dat)],labels=legend_labels)
  p1=p1+theme(legend.position =legend_position)

  #Dont include a legend if there is only one dataset and no legend label
  if((length(dat)<2) & legend_set2color){
    p1=p1+theme(legend.position ="none")
  }
  return(p1)
}
