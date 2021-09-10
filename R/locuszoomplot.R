
#' Locuszoom plot
#'
#' @description
#'
#' \code{locuszoom()} displays the association results for a smaler region within one chromosome
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P} in upper or lowercase)
#'
#' @inheritParams regionplot
#'
#' @return plots using egg (https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html)
#' @export
#'
#' @examples
#' \dontrun{
#' locuszoom(df)
#' }


locuszoom <- function(df, ntop=3, xmin=0, size=2, shape=19, alpha=1,label_size=4, annotate_with="ID",
                      color=NULL, axis_text_size=11,axis_title_size=12,title_text_size=13, show_genes=FALSE, show_overview=F,
                      show_exons=FALSE, max_genes=200, sign_thresh=5e-09, sign_thresh_color="red", sign_thresh_label_size=3.5,
                      xmax=NULL,ymin=NULL,ymax=NULL,protein_coding_only=FALSE,region_size=1000000,gene_padding=100000,angle=0,legend_title_size=12,legend_text_size=12,
                      nudge_x=0.01,nudge_y=0.01, rsids=NULL, variant=NULL,rsids_color="gray40",legend_name="Data:",legend_position="right",
                      chr=NULL,vline=NULL,show_gene_names=NULL,legend_labels=NULL,gene=NULL, title=NULL, label_color="gray40",region=NULL){
  dat <- dat_check(df)
  if(length(unique(dat[[1]]$CHROM)) > 1){
    stop("There are multiple chromosomes in the input dataset. The locuszoom plot only works for a small genetic region_size within one chromosome")
  }
  for(i in seq_along(dat)){
    df <- as.data.frame(dat[[i]])
    if(! "R2" %in% colnames(df)){
      stop("Aborting. Could not find the R2 column in the input data!! ")
    }
  }
  if(is.null(gene) & is.null(variant) & is.null(region) & (is.null(chr)  & is.null(xmax))){
    xmin <- min(dat[[1]]$POS)
    xmax <- max(dat[[1]]$POS)
    chr <- dat[[1]]$CHROM
     #get xmin and max from the input dataframe, since we expect the dataset provided to be of a smaller region
    for(i in 2:length(dat)){
      tmp_xmin <- min(dat[[1]]$POS)
      tmp_xmax <- max(dat[[1]]$POS)
      xmin <- ifelse(tmp_xmin < xmin, tmp_xmin, xmin)
      xmax <- ifelse(tmp_xmax < xmax, tmp_xmax, xmax)

    }
    xmin <- xmin - 100
    xmax <- xmax + 100
  }
  df <- set_lz_colors(dat)

  regionplot(df, ntop=ntop, annotate=NULL, xmin=xmin, size=size, shape=shape, alpha=alpha,label_size=label_size, annotate_with=annotate_with,
             color=color, axis_text_size=axis_text_size,axis_title_size=axis_title_size,title_text_size=title_text_size, show_genes=show_genes, show_overview=show_overview,
             show_exons=show_exons, max_genes=max_genes, sign_thresh=sign_thresh, sign_thresh_color=sign_thresh_color, sign_thresh_label_size=sign_thresh_label_size,
             xmax=xmax,ymin=ymin,ymax=ymax,protein_coding_only=protein_coding_only,region_size=region_size,gene_padding=gene_padding,angle=angle,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
             nudge_x=nudge_x,nudge_y=nudge_y, rsids=rsids, variant=variant,rsids_color=rsids_color,legend_name=legend_name,legend_position=legend_position,
             chr=chr,vline=vline,show_gene_names=show_gene_names,legend_labels=legend_labels,gene=gene, title=title, label_color=label_color,locuszoomplot=T,region=region)
}

set_lz_colors <- function(dat){
  for(i in seq_along(dat)){
    df <- dat[[i]]
    df$color <- "darkblue"
    df$color <- ifelse(df$R2 == 1, "purple", df$color)
    df$color <- ifelse(df$R2 != 1 & df$R2 > 0.8, "red", df$color)
    df$color <- ifelse(df$R2 < 0.8 & df$R2 > 0.6, "orange", df$color)
    df$color <- ifelse(df$R2 < 0.6 & df$R2 > 0.4, "green", df$color)
    df$color <- ifelse(df$R2 < 0.4 & df$R2 > 0.2, "turquoise", df$color)
    dat[[i]] <- df
  }
  return(dat)
}

