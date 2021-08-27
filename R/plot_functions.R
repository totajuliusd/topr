
change_axes <- function(p1){
  p1+theme(panel.border=element_blank(),axis.line=element_line(color="grey"))
}
rm_grids <- function(p1){
  p1 <- p1+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  return(p1)
}
rm_axis_text_and_ticks <- function(p1){
  #hide labels and ticks
  p1 <- p1+theme(axis.title = element_blank(),
                 axis.text=element_blank(),axis.ticks = element_blank())
  return(p1)
}

set_axis_labels <- function(p1,xaxis_label="Chromosome"){
  p1 <- p1 + labs(x=xaxis_label, y=expression(-log[10](italic(p))))
  return(p1)
}

set_plot_text_sizes <- function(p1, axis_text_size=12, axis_title_size=12, legend_title_size=12,legend_text_size=12){
  p1 <- p1 + theme(axis.text=element_text(size=axis_text_size),
              axis.title = element_text(size=axis_title_size),
              legend.title=element_text(size=legend_title_size),
              legend.text=element_text(size=legend_text_size))
  return(p1)
}
add_title <- function(p1, title="", title_text_size=14){
  p1 <- p1 + ggtitle(title)
  p1 <- p1 + theme(plot.title = element_text(size=title_text_size))
  return(p1)
}


set_legend_texts <- function(p1, legend_title_size=12,legend_text_size=12){
  p1 <- p1 + theme(legend.title=element_text(size=legend_title_size),
                   legend.text=element_text(size=legend_text_size))
  return(p1)
}

add_genes2plot <- function(p1,genes,highlight_genes_color="green",highlight_genes_ypos=1){
  p1 <- p1 + geom_point(data=genes, aes(x=POS, y=highlight_genes_ypos),color=highlight_genes_color, size=2, shape=15) +
    ggrepel::geom_text_repel(data=genes, aes(x=POS,y=highlight_genes_ypos,label=Gene_Symbol), color="black",size=3,direction="both",nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5)
  return(p1)
}

add_zero_hline <- function(p1){
  return(p1 + geom_hline(yintercept = 0,size=0.2, color="grey54"))
}

add_shades_and_ticks <- function(p1, shades, ticks){
  p1 <- p1 + geom_rect(data=shades, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1),color="#D3D3D3", alpha=0.1, size=0.1)+
     scale_x_continuous(breaks=ticks$pos, labels=ticks$names,expand=c(.01,.01)) + scale_y_continuous(expand=c(.02,.02))
  return(p1)
}

add_annotation <- function(p1,plot_labels=NULL, nudge_x=0.01, nudge_y=0.01, label_size=3.5, angle=0,annotate_with="Gene_Symbol"){
  if(! is.null(plot_labels)){
    p1 <- p1 + ggrepel::geom_text_repel(data=plot_labels, aes(x=POS, y=log10p, label=(plot_labels %>% dplyr::pull(annotate_with)) ),
                                    nudge_x=nudge_x,nudge_y=ifelse(plot_labels$log10p>0, nudge_y, -nudge_y),size=label_size,
                                    segment.size=0.2, color=plot_labels$color,segment.color = "black",
                                    max.iter=10000,direction="both",angle=angle)
  }
  return(p1)
}


add_zoom_rectangle <- function(p1,dat,xmin=NULL,xmax=NULL){
  #add the rectangle to the overview plot
  if(!is.null(xmax) && ! is.null(xmin)){
    ymin <- get_ymin(dat)
    ymax <- get_ymax(dat)
  }
  p1 <- p1+geom_rect(mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), size=0.2, alpha=0.3, color="red", fill=NA)
  return(p1)
}

add_vline <- function(p1, vline){
  if(!is.null(vline)){
    for(i in seq_along(vline)){
      p1 <- p1+geom_vline(xintercept = as.numeric(vline[i]) , colour="grey", linetype="dashed")
    }
  }
  return(p1)
}

add_annot <- function(p1, angle=0, annotate_with="Gene_Symbol",top_snps=NULL){
  if(! is.null(top_snps)){
    p1 <- p1 %>% add_annotation(top_snps,angle=angle, annotate_with=annotate_with)
  }
  return(p1)
}


add_sign_thresh <- function(p1, sign_thresh=1e-09,sign_thresh_color="red",size=0.5, using_ntop=FALSE){
  #add significance threholds to the plot
  if(is.vector(sign_thresh)){
    for(i in seq_along(sign_thresh)){
      color <- "red"
      if(!is.null(sign_thresh_color) & (length(sign_thresh) == length(sign_thresh_color))){
        color <- sign_thresh_color[i]
      }
      p1 <- p1 + geom_hline(yintercept = -log10(as.numeric(sign_thresh[[i]])),size=size, color=color, linetype="dashed")
      if(using_ntop){
        p1 <- p1 + geom_hline(yintercept = log10(as.numeric(sign_thresh[[i]])),size=size, color=color, linetype="dashed")
      }
    }
  }
  else{
    p1 <- p1+geom_hline(yintercept = -log10(as.numeric(sign_thresh[[i]])),size=size, colour="red", linetype="dashed")
    if(using_ntop){
      p1 <- p1 + geom_hline(yintercept = log10(as.numeric(sign_thresh[[i]])),size=size, color=color, linetype="dashed")
    }
  }
  return(p1)
}


add_sign_thresh_labels <- function(p1, sign_thresh=1e-09,sign_thresh_color="red",xmin=0,sign_thresh_label_size=3.5){
  #add significance threholds label to the plot
  if(is.null(xmin)){
    xmin <- 0
  }
  tmpdf <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("color", "label", "ypos","xpos"))))
  for(i in seq_along(sign_thresh)){
    color <- "red"
    if(!is.null(sign_thresh_color) & (length(sign_thresh) == length(sign_thresh_color))){
      color <- sign_thresh_color[i]
    }
    tmpdf <- rbind(tmpdf, data.frame(color=color, label=sign_thresh[[i]],ypos=-log10(sign_thresh[[i]])* 1.02,xpos=xmin))

  }
  p1 <- p1+ggrepel::geom_text_repel(data=tmpdf, aes(x=xpos, y=ypos, label=label,color=color), nudge_y=0.02,size=sign_thresh_label_size,
                                 max.iter=10000,direction="y")


  return(p1)
    #p1+geom_text(data=tmpdf, aes(x=xlabelpos,y=ypos,label=label,color=color),size=sign_thresh_label_size))
}

add_rsids <- function(p1,dat,rsids,rsids_color="gray40"){
  rsids_df=data.frame(POS=numeric(),log10p=numeric,ID=character())
  for(i in seq_along(dat)){
    ids_found <- dat[[i]] %>% filter(ID %in% rsids)
    rsids_df <- rbind(rsids_df, ids_found)
  }
  p1 <- p1+ggrepel::geom_text_repel(data=rsids_df, aes(x=POS, y=log10p, label=ID,color=rsids_color), nudge_y=0.02,max.iter=10000,direction="both")
  return(p1)
}

set_ymin_ymax <- function(p1,ymin,ymax){
  return(p1+coord_cartesian(ylim=c(ymin,ymax)))
}
set_xmin_xmax <- function(p1,xmin,xmax){
  return(p1+coord_cartesian(xlim=c(xmin,xmax)))
}

set_xmin_xmax_ymin_ymax <- function(p1,xmin,xmax,ymin,ymax){
  return(p1+coord_cartesian(xlim=c(xmin,xmax),ylim=c(ymin,ymax),))
}