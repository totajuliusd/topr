
get_exonplot_coords <- function(df, xmin, xmax){
  
  level_count <- 0
  #first check whether there are any genes/exons in the region
  if(nrow(df)>0){
    min_space_between_genes <- as.numeric((xmax-xmin)*0.1)
    x1 <- vector(); x2 <- vector();  y1 <- vector(); y2 <- vector();  df$y <-NA
    x1_exon <- vector(); x2_exon <- vector(); y1_exon <- vector(); y2_exon <- vector(); biotype <- vector();
    current_level <- 1
    gene_struct <- data.frame(gene_end=integer(),level=integer(),gene=character())
    for(i in seq_along(df$gene_symbol)){
      gene_start <- df$gene_start[i]
      gene_end <- df$gene_end[i]
      gene_name <- df$gene_symbol[i]
      x1 <- append(gene_start,x1)
      x2 <- append(gene_end, x2)
      add_level_to_gene_struct <- 1
      if(length(gene_struct$level)>0){
        for(j in seq_along(gene_struct$level)){
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
      y1 <- append(current_level,y1)
      y2 <- append(current_level,y2)
      
      for(k in seq_along(unlist(strsplit(df$exon_chromstart[i], ",")))){
        x1_exon <- as.numeric(append(unlist(strsplit(df$exon_chromstart[i], ","))[k], x1_exon))
        x2_exon <- as.numeric(append(unlist(strsplit(df$exon_chromend[i], ","))[k], x2_exon))
        y1_exon <- append(current_level+0.3, y1_exon)
        y2_exon <- append(current_level-0.3, y2_exon)
        biotype <- append(df$biotype[i],biotype)
      }
      df$y[i] <- current_level+0.35
    }
    gene_coords <- data.frame(x1=x1, x2=x2, y1=y1,y2=y2)
    shades <- data.frame(x1=x1_exon, x2=x2_exon,y1=y1_exon,y2=y2_exon,biotype=biotype)
  } else {
    gene_coords <- data.frame()
    shades <- data.frame()
  }
  return(list(shades=shades,genes=df,gene_coords=gene_coords,level_count=level_count))
}

