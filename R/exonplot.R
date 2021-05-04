
exonplot=function(df, xmin, xmax, label_size=3.5,vline=NULL){
  min_space_between_genes=as.numeric((xmax-xmin)*0.1)
  x1=vector()
  x2=vector()
  y1=vector()
  y2=vector()
  df$y=NA
  x1_exon=vector()
  x2_exon=vector()
  y1_exon=vector()
  y2_exon=vector()
  biotype=vector()

  current_level=1
  level_count=0
  gene_struct = data.frame(gene_end=integer(),level=integer(),gene=character())
  chr=gsub("chr","", unique(df$chrom))

  for(i in 1:length(df$gene_symbol)){
    gene_start=df$gene_start[i]
    gene_end=df$gene_end[i]
    gene_name=df$gene_symbol[i]
    x1=append(gene_start,x1)
    x2=append(gene_end, x2)
    add_level_to_gene_struct=1
    if(length(gene_struct$level)>0){
      for(j in 1:length(gene_struct$level)){
        if(gene_start > (gene_struct$gene_end[j]+min_space_between_genes)){
          if(add_level_to_gene_struct){
            current_level=gene_struct$level[j]
            #reset the gene end at this level
            gene_struct$gene_end[j]=gene_end
            add_level_to_gene_struct=0
          }
        }
      }
    }
    if(add_level_to_gene_struct){
      level_count=level_count+1
      gene_struct=rbind(gene_struct, data.frame(level=level_count, gene_end=gene_end, gene=gene_name))
      current_level=level_count

    }
    y1=append(current_level,y1)
    y2=append(current_level,y2)
    # print(paste("gene ",gene_start,"-",gene_end," ",gene_name))
    #  print(pase("exon", ))
    for(k in 1:length(unlist(strsplit(df$exon_chromstart[i], ",")))){
      x1_exon=as.numeric(append(unlist(strsplit(df$exon_chromstart[i], ","))[k], x1_exon))
      x2_exon=as.numeric(append(unlist(strsplit(df$exon_chromend[i], ","))[k], x2_exon))
      y1_exon=append(current_level+0.3, y1_exon)
      y2_exon=append(current_level-0.3, y2_exon)
      biotype=append(df$biotype[i],biotype)
    }
    df$y[i]=current_level+0.35
  }
  genes=data.frame(x1=x1, x2=x2, y1=y1,y2=y2)
  exons=data.frame(x1=x1_exon, x2=x2_exon,y1=y1_exon,y2=y2_exon)
  p1=ggplot(data =  df)+theme_bw()



  p1=p1+geom_text(data=df, aes(x=gene_start,y=as.numeric(y),label=gene_symbol),hjust=0, vjust=0, size=label_size)
  p1=p1+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
  p1=p1+xlab(paste("Position on chr",chr, sep=""))
  p1=p1+scale_y_discrete(breaks=NULL)

  p1=p1+geom_segment(data=genes, aes(x=x1,y=y1,xend=x2,yend=y2))
  p1=p1+geom_rect(data=exons, mapping=aes(ymax=y2,xmin=x1, xmax=x2, ymin=y1,fill=biotype))+labs(fill="Biotype")
  p1=p1+suppressMessages(coord_cartesian(xlim=c(xmin, xmax), ylim=c(1,level_count+0.2)))
  if(!is.null(vline)){
    p1=p1+geom_vline(xintercept =vline , colour="grey", linetype="dashed")
  }
  return(p1)
}
