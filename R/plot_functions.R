
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

set_plot_text_sizes <- function(p1, axis_text_size=12, axis_title_size=12, legend_title_size=12,legend_text_size=11,scale=1){
  p1 <- p1 + theme(axis.text=element_text(size=axis_text_size*scale),
              axis.title = element_text(size=axis_title_size*scale),
              legend.title=element_text(size=legend_title_size*scale),
              legend.text=element_text(size=legend_text_size*scale))
  return(p1)
}
add_title <- function(p1, title="", title_text_size=14,scale=1){
  p1 <- p1 + ggtitle(title)
  p1 <- p1 + theme(plot.title = element_text(size=title_text_size*scale))
  return(p1)
}


set_legend_texts <- function(p1, legend_title_size=13,legend_text_size=12,scale=1){
  p1 <- p1 + theme(legend.title=element_text(size=legend_title_size*scale),
                   legend.text=element_text(size=legend_text_size*scale))
  return(p1)
}

add_genes2plot <- function(p1,genes,highlight_genes_color="green",highlight_genes_ypos=1, label_size=3,gene_label_angle=0,scale=scale,gene_label_fontface="plain",gene_label_family=""){
  p1 <- p1 + geom_point(data=genes, aes(x=POS, y=highlight_genes_ypos),color=highlight_genes_color, size=2*scale, shape=15) +
    ggrepel::geom_text_repel(data=genes, aes(x=POS,y=highlight_genes_ypos,label=Gene_Symbol), fontface=gene_label_fontface,family=gene_label_family,color="black",size=label_size*scale,direction="both",angle=gene_label_angle,nudge_x = 0.01,nudge_y = 0.01,segment.size=0.2,segment.alpha =.5)
  return(p1)
}

add_zero_hline <- function(p1){
  return(p1 + geom_hline(yintercept = 0,size=0.2, color="grey54"))
}

add_shades_and_ticks <- function(p1, shades, ticks,shades_color=NULL,shades_alpha=0.5, shades_line_alpha=1,theme_grey=F){
   if(theme_grey){
     shades_color="#D3D3D3"
     shades_alpha=0.1
     p1 <- p1 + geom_rect(data=shades, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1),
                         color=alpha(shades_color,shades_line_alpha),  alpha=shades_alpha, size=0.1)+
      scale_x_continuous(breaks=ticks$pos, labels=ticks$names,expand=c(.01,.01)) + scale_y_continuous(expand=c(.02,.02))
  }
  else{
    if(is.null(shades_color))
       shades_color="white"
    p1 <- p1 + geom_rect(data=shades, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1),
                         color=alpha(shades_color,shades_line_alpha),
                         fill=alpha(shades_color, shades_alpha), size=0.1)+
   scale_x_continuous(breaks=ticks$pos, labels=ticks$names,expand=c(.01,.01)) + scale_y_continuous(expand=c(.02,.02))
  }
   return(p1)
}

add_annotation <- function(p1,plot_labels=NULL, nudge_x=0.01, nudge_y=0.01, label_size=3.5, angle=0,annotate_with="Gene_Symbol", label_color=NULL,scale=1, 
                           segment.size=0.2,segment.color="black",segment.linetype="solid",max.overlaps=10){
  if(! is.null(label_color)){
     if(is.vector(label_color) & length(label_color) > 1){
      label_color <- label_color[1]; 
      print("label color can only be assigned one color, so using the first color from the vector!  The arguments label_alpha,label_font and label_family take vectors as input and can be used to distinguish between labels.")
    }
     if(nrow(plot_labels) >0 ) plot_labels$color <- label_color
    }
   if(! is.null(plot_labels)){
    if(is.null(label_color) & length(unique(plot_labels$color)) == 1){plot_labels$color="black"} 
        p1 <- p1 + ggrepel::geom_text_repel(data=plot_labels, 
                                        aes(x=POS, y=log10p, label=(plot_labels %>% dplyr::pull(annotate_with)) ),
                                        nudge_x=plot_labels$nudge_x,nudge_y=ifelse(plot_labels$log10p>0, plot_labels$nudge_y, -plot_labels$nudge_y),
                                        segment.size=segment.size,size=label_size*scale, fontface=plot_labels$fontface, family=plot_labels$family, alpha=plot_labels$alpha,
                                       color=plot_labels$color,
                                        segment.color = segment.color,
                                        segment.linetype=segment.linetype,max.iter=10000,direction="both",angle=plot_labels$angle, max.overlaps = max.overlaps)
    
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

add_vline <- function(p1, vline, vline_color="grey",vline_linetype="dashed", vline_alpha=1, vline_size=1, scale=1){
  if(!is.null(vline)){
    for(i in seq_along(vline)){
      p1 <- p1+geom_vline(xintercept = as.numeric(vline[i]) , colour=vline_color, linetype=vline_linetype, alpha=vline_alpha, size=vline_size*scale)
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


add_sign_thresh <- function(p1, sign_thresh=1e-09,sign_thresh_color="red", using_ntop=FALSE, sign_thresh_linetype="dashed", sign_thresh_size=0.5, scale=1){
  #add significance threholds to the plot
  if(is.vector(sign_thresh)){
    for(i in seq_along(sign_thresh)){
      color <- "red"
      if(!is.null(sign_thresh_color) & (length(sign_thresh) == length(sign_thresh_color))){
        color <- sign_thresh_color[i]
      }
      p1 <- p1 + geom_hline(yintercept = -log10(as.numeric(sign_thresh[[i]])),size=sign_thresh_size*scale, color=color, linetype=sign_thresh_linetype)
      if(using_ntop){
        p1 <- p1 + geom_hline(yintercept = log10(as.numeric(sign_thresh[[i]])),size=sign_thresh_size*scale, color=color, linetype=sign_thresh_linetype)
      }
    }
  }
  else{
    p1 <- p1+geom_hline(yintercept = -log10(as.numeric(sign_thresh[[i]])),size=sign_thresh_size*scale, colour="red", linetype=sign_thresh_linetype)
    if(using_ntop){
      p1 <- p1 + geom_hline(yintercept = log10(as.numeric(sign_thresh[[i]])),size=sign_thresh_size*scale, color=color, linetype=sign_thresh_linetype)
    }
  }
  return(p1)
}


add_sign_thresh_labels <- function(p1, sign_thresh=1e-09,sign_thresh_color="red",xmin=0,sign_thresh_label_size=3.5,scale=1){
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
  p1 <- p1+ggrepel::geom_text_repel(data=tmpdf, aes(x=xpos, y=ypos, label=label,color=color), nudge_y=0.02,size=sign_thresh_label_size*scale,
                                 max.iter=10000,direction="y")


  return(p1)
}

get_rsids_from_df <- function(dat, rsids){
  rsids_df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("POS", "log10p", "ID","color"))))
  if(!is.vector(rsids)){
    rsids <- c(rsids)
  }
  for(i in seq_along(dat)){
    ids_found <- dat[[i]] %>% dplyr::filter(ID %in% rsids) %>% dplyr::select("POS","log10p","ID","P","color")
    rsids_df <- base::rbind(rsids_df, ids_found)
  }
  if(length(rsids_df$POS) == 0){
    print("Could not find any of the requested rsids")
  }
  else{
    rsids_df <- rsids_df %>% dplyr::arrange(P) %>% distinct(ID,color,.keep_all=T) #%>% dplyr::distinct("ID", "color", .keep_all=T)
  }
  if(length(rsids_df$POS) < length(rsids)){
    print("Could not find all the requested rsids!! ")
  }
  return(rsids_df)
}


add_rsids <- function(p1,rsids_df, rsids_color=NULL, nudge_x=0.01, nudge_y=0.01, label_size=3.5, angle=0,label_color=NULL, scale=1, with_vline=F, label_fontface="plain",label_family="",segment.size=0.2,segment.color="black",segment.linetype="solid"){
  if(!is.null(label_color)){
    rsids_df$color <- label_color
  }
  if(!is.null(rsids_color)){
    rsids_df$color <- rsids_color
  }
  p1 <- p1+ggrepel::geom_text_repel(data=rsids_df, aes(x=POS, y=log10p, label=ID,color=color), nudge_x=nudge_x, nudge_y=nudge_y, size=label_size*scale, 
                                    angle=angle,max.iter=10000,direction="both",fontface=label_fontface, family=label_family)
  if(with_vline){
    p1 <- p1 %>% add_vline(rsids_df$POS)
  }
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