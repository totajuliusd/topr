

get_base_plot <- function(dat, color=get_topr_colors(), show_legend=TRUE, legend_labels=NULL,scale=1,
                          legend_name=NULL,overview_plot=FALSE,legend_position="right",locuszoomplot=FALSE, legend_nrow=NULL,verbose=NULL){
  p1 <- ggplot()+theme_bw() #+geom_point(data=dat[[1]]$gwas, aes(dat[[1]]$gwas$pos_adj, dat[[1]]$gwas$log10p),color=colors[1] alpha=0.7,size=1)+theme_bw()
  for(i in seq_along(dat)){
    if(overview_plot){
      dat[[i]]$size <- 0.4
    }
     p1 <- p1+geom_point(data=dat[[i]], aes(x=POS, y=log10p, color=color), alpha=dat[[i]]$alpha, size=dat[[i]]$size*scale,shape=dat[[i]]$shape)
  }
  #add the legend for the different datasets
  if(! is.null(legend_name)) legend_name <- legend_name
  if(is.null(legend_labels) & show_legend){
    if(length(dat) > 1){
      legend_labels <- color[seq_along(dat)]
      if(is.null(verbose)){
        print("Use the legend_labels argument to change the legend labels from color names to meaningful labels! ")
      }else if(verbose){
        print("Use the legend_labels argument to change the legend labels from color names to meaningful labels! ")
      }
    }
  }
  if(locuszoomplot & show_legend){
    colors <-c("darkblue","turquoise","green","orange","red")
    p1 <- p1+scale_color_identity(guide = "legend", name="R2",  breaks=colors, labels=c("R2 < 0.2", "0.2 < R2 < 0.4", "0.4 < R2 < 0.6","0.6 < R2 <0.8", "0.8 < R2"))
  }
  else if(length(dat) == 1 & is.null(legend_labels)) { #no need to use a legend when there is only one phenotype on display
    p1 <- p1+scale_color_identity(breaks=color[seq_along(dat)])
  }
  else if(show_legend){
    p1 <- p1+scale_color_identity(guide = "legend", name=legend_name,  color[seq_along(dat)], labels=legend_labels)+
      theme(legend.position = legend_position)
  }
  if(! is.null(legend_nrow)){
    p1 <- p1 + guides(color=guide_legend(nrow=legend_nrow,byrow=TRUE))
  }
  p1 <- rm_grids(p1)
  return(p1)
}


get_gene_plot <- function(df,label_size=3.5,xmin=xmin,xmax=xmax, show_gene_names=NULL,scale=1,show_gene_legend=TRUE, gene_color=NULL,gene_label_fontface="plain",gene_label_family=""){
   if(is.null(gene_color)){
       color="#00AFBB"
  }else{
    color=gene_color
    show_gene_legend=FALSE
  }
  p1 <- ggplot(data=df$genes)+theme_bw()
  if(is.null(show_gene_names)){  #if the user wants to view the gene names even though the region is large, he/she can do so by setting show_gene_names to TRUE
    if(xmax-xmin > 10000000){
      show_gene_names <- FALSE
    } else{ show_gene_names <- TRUE}
  }
  if(nrow(df$genes)>0){
    if(show_gene_names){
      p1 <- p1+geom_text(data=df$genes, aes(x=gene_start,y=as.numeric(y),label=gene_symbol,fontface=gene_label_fontface,family=gene_label_family),hjust=0, vjust=0, size=label_size*scale)
    }
    if(length(df) == 4){ # show genes as lines and exons as rectangles
      p1 <- p1+geom_segment(data=df$gene_coords, aes(x=x1,y=y1,xend=x2,yend=y2))
    }
    if(show_gene_names){
       if(show_gene_legend){
        p1 <- p1+geom_rect(data=df$shades, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1,fill=biotype))+labs(fill="Biotype")
      }
      else{
         p1 <- p1+geom_rect(data=df$shades, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1), color=color, fill=color)
       }
    } #make the genes a bit thicker, since there doesnt have to be room for the text
    else{
      if(show_gene_legend){
        p1 <- p1+geom_rect(data=df$shades, mapping=aes(ymax=y2+1,xmin=x1, xmax=x2, ymin=y1+1,fill=biotype))+labs(fill="Biotype")
      }
      else{
        p1 <- p1+geom_rect(data=df$shades, mapping=aes(ymax=y2+1,xmin=x1, xmax=x2, ymin=y1+1), color=color, fill=color)
      }
    }
  }
  p1 <- p1+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
  p1 <- p1+coord_cartesian(xlim=c(xmin, xmax), ylim=c(0.5,df$level_count+1))
  p1 <- rm_grids(p1)
  return(p1)
}


