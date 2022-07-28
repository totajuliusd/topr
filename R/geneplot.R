
get_geneplot_coords <- function(genes,xmin=xmin, xmax=xmax){
  level_count <- 0
  #first check whether there are any genes/exons in the region
  if(nrow(genes)>0){
    min_space_between_genes <-as.numeric((xmax-xmin)*0.1)
    x1 <- vector(); x2 <- vector(); y1 <- vector(); y2 <- vector();biotype <- vector()
    genes$y <- NA
    current_level <- 1
    gene_struct <-  data.frame(gene_end=integer(),level=integer(),gene=character())
    chr <- gsub("chr","", unique(genes$chrom))
    for(i in seq_along(genes$gene_symbol)){
      gene_start <- genes$gene_start[i]; gene_end <- genes$gene_end[i]; gene_name <- genes$gene_symbol[i]
      x1 <- append(gene_start,x1)
      x2 <- append(gene_end, x2)
      biotype <- append(genes$biotype[i], biotype)
      add_level_to_gene_struct <- 1
      if(length(gene_struct$level)>0){
      for(j in seq_along(gene_struct$level)){
        #if(gene_start > (gene_struct$gene_end[j]+min_space_between_genes)){
        if(gene_start > (gene_struct$gene_end[j]+min_space_between_genes)){
          if(add_level_to_gene_struct){
            current_level <- gene_struct$level[j]
            #reset the gene end at this level
            gene_struct$gene_end[j] <- gene_end
            add_level_to_gene_struct <- 0
          }
        }
      }
    }
    if(add_level_to_gene_struct){
      level_count <- level_count+1
      gene_struct <- rbind(gene_struct, data.frame(level=level_count, gene_end=gene_end, gene=gene_name))
      current_level <- level_count

    }
    y1 <- append(current_level-0.20,y1)
    y2 <- append(current_level+0.20,y2)
    genes$y[i] <- current_level+0.25
  }
  }
  shades <- NULL
  if(nrow(genes)>0){
    shades <- data.frame(x1=x1, x2=x2, y1=y1,y2=y2,biotype=biotype)
  }
  return(list(shades=shades,genes=genes, level_count=level_count))
}

plot_exons <- function(shades=shades,genes=genes,genes_coords=genes_coords, level_count=level_count) { #Do the plotting

  p1 <- ggplot(data =  df)+theme_bw()
  if(nrow(df)>0){
    p1 <- p1+geom_text(data=genes, aes(x=gene_start,y=as.numeric(y),label=gene_symbol),hjust=0, vjust=0, size=label_size)
    p1 <- p1+geom_segment(data=genes_coords, aes(x=x1,y=y1,xend=x2,yend=y2))
    p1 <- p1+geom_rect(data=shades, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1,fill=biotype))+labs(fill="Biotype")

  }
  p1 <- p1+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
  p1 <- p1+xlab(paste("Position on chr",chr, sep=""))
  #p1=p1+scale_y_discrete(breaks=NULL)

  p1 <- p1+suppressMessages(coord_cartesian(xlim=c(xmin, xmax), ylim=c(1,level_count+1)))
  if(!is.null(vline)){
    p1 <- p1+geom_vline(xintercept =vline , colour="grey", linetype="dashed")
  }
  return(p1)
}
