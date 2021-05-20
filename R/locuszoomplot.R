#' Locuszoom plot
#'
#' @description
#'
#' \code{locuszoom()} displays the association results for a smaler region within one chromosome
#' Required parameter one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P,R2} in upper or lowercase), xmin and xmax
#'
#' All other input parameters are optional
#'
#' @inheritParams regionplot
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' locuszoomt(df,xmin=xmin,xmax=xmax)
#' }

#guides(colour = guide_colourbar(order = 1),
#http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software

locuszoom=function(df,ids2label,xmin=0,xmax=max(df$POS),label_size=3.5, gene_label_size=3.5, genes=NULL,vline=NULL,sign_thresh=NULL){

  lz_plot=locuszoom_plot2(df, ids2label,xmin=xmin,xmax=xmax, label_size=label_size, sign_thresh=sign_thresh,vline=vline)

  if("exon_chromstart" %in% colnames(genes))
    gene_plot=exonplot(genes, xmin, as.numeric(xmax),gene_label_size=gene_label_size,vline=vline)
  else
    gene_plot=geneplot(genes, xmin, as.numeric(xmax),gene_label_size=gene_label_size,vline=vline)
  g1=ggplotGrob(lz_plot)
  g2=ggplotGrob(gene_plot)
  fg1=gtable_frame(g1, height=unit(7, "null"), debug=F)
  fg2=gtable_frame(g2, height=unit(2, "null"),debug=F)
  fg12=gtable_frame(gtable_rbind(fg1,fg2),debug=FALSE)
  grid.newpage()
  grid.draw(fg12)
}



locuszoom_plot2=function(df, ids2label,xmin=0, xmax=max(df$POS),label_size=3.5, vline=NULL,sign_thresh=NULL,axis_text_size=10,axis_title_size=12,title_text_size=13,legend_title_size=12,legend_text_size=12){
  colors=c("darkblue","turquoise","green","orange","red")
  # 0.2, 0.4, 0.6, 0.8
  if(!is.null(xmin)){
    df=df %>% filter(POS>xmin) %>% filter(POS<xmax)
  }
  #retrieve the IDs fro
  zoom2=df %>%  filter(ID==ids2label[1])%>% arrange(P)%>% head(1)

  df$color="darkblue"
  df$color=ifelse(df$R2 > 0.8, "red", df$color)
  df$color=ifelse(df$R2 < 0.8 & df$R2 > 0.6, "orange", df$color)
  df$color=ifelse(df$R2 < 0.6 & df$R2 > 0.4, "green", df$color)
  df$color=ifelse(df$R2 < 0.4 & df$R2 > 0.2, "turqoise", df$color)

  colorMap=c("0.8 < r2 < 1"="red", "0.6 < r2 < 0.8"="orange","0.4 < r2 < 0.6"="green","0.2 < r2 < 0.4"="turquoise","0 < r2 < 0.2"="darkblue")

  p1=ggplot(data = df)+theme_bw()

  p1=p1+geom_point(data=df, aes(POS, -log10(P),color=color), size=ifelse(df$POS==zoom2$POS,3,2), shape=ifelse(df$POS==zoom2$POS,18,19))

  if(!is.null(vline))
    p1=p1+geom_vline(xintercept =vline , colour="grey", linetype="dashed")

  zoom_tmp=df[FALSE,]
  for(i in ids2label){
    zoom2 = df %>%  filter(ID==i) %>% arrange(P) %>% head(1)
    p1=p1+ggrepel::geom_text_repel(data=unique(zoom2), aes(x=POS,y=-log10(P),label=ID), color="grey40",size=label_size,direction="both",nudge_x = 0.5,nudge_y = 0.5,segment.size=0.2,segment.alpha =.5)
    zoom_tmp=rbind(zoom_tmp, zoom2) %>% distinct(ID,POS, .keep_all = T)
  }
  p1=p1+ggrepel::geom_text_repel(data=unique(zoom_tmp), aes(x=POS,y=-log10(P),label=ID), color="grey40",size=label_size,direction="both",nudge_x = 0.5,nudge_y = 0.5,segment.size=0.2,segment.alpha =.5)

  #if(!is.null(sign_thresh))
  #  p1=p1+geom_hline(yintercept = -log10(sign_thresh), colour="red", linetype="dashed")

  p1=tidy_plot(p1, axis_text_size=axis_text_size,axis_title_size=axis_title_size, title_text_size=title_text_size,legend_title_size=legend_title_size,legend_text_size=legend_text_size)

  p1=p1+theme(axis.text.y=element_text(size=10), title = element_text(size=11))
  p1=p1+xlim(xmin, xmax)+ ylim(min(-log10(df$P)), max(-log10(df$P)*1.1))

  #Add the legend
  p1=p1+scale_color_identity(guide = "legend", name="R2",  breaks=colors, labels=c("R2 < 0.2", "0.2< R2 < 0.4", "0.4 < R2 < 0.6","0.6 < R2 <0.8", "0.8 < R2"))
  p1=p1+theme(legend.position = "right")

  return(p1)
}


locuszoom_plot=function(df, ids2label,xmin=0, xmax=max(df$POS),label_size=3.5, vline=NULL,sign_thresh=NULL){
  df$color="darkblue"
  df$color=ifelse(df$R2 > 0.2, "turquoise", df$color)
  df$color=ifelse(df$R2 > 0.4, "green", df$color)
  df$color=ifelse(df$R2 > 0.6, "orange", df$color)
  df$color=ifelse(df$R2 > 0.8, "red", df$color)

  if(!is.null(xmin)){
    df=df %>% filter(POS>xmin) %>% filter(POS<xmax)
  }

  zoom2=df %>%  filter(ID==ids2label[1])%>% arrange(P)%>% head(1)


  colorMap=c("0.8 < r2 < 1"="red", "0.6 < r2 < 0.8"="orange","0.4 < r2 < 0.6"="green","0.2 < r2 < 0.4"="turquoise","0 < r2 < 0.2"="darkblue")

  p1=ggplot(data = df)+theme_bw()

  #remove space between axis and plots
  #remove space between axis and plots
  p1=p1+scale_y_continuous(expand=c(.02,.02))
  p1=p1+scale_x_continuous(expand=c(.01,.01),labels = scales::comma )

  if(!is.null(vline)){
    p1=p1+geom_vline(xintercept =vline , colour="grey", linetype="dashed")
  }
  # p1=p1+geom_point(aes(POS, -log10(P),color=colorMap), size=ifelse(df$POS==zoom2$POS,3,2), color=ifelse(df$POS==zoom2$POS, "purple", df$color), shape=ifelse(df$POS==zoom2$POS,18,19))
  tmp=df %>% filter(R2>0.8)
  p1=p1+geom_point(data=tmp, aes(POS, -log10(P),colour="0.8 < r2 < 1"), size=ifelse(tmp$POS==zoom2$POS,3,2), shape=ifelse(tmp$POS==zoom2$POS,18,19))
  tmp=df %>% filter(R2<0.8 & R2>0.6)
  if(length(tmp$POS)> 0){ p1=p1+geom_point(data=tmp, aes(POS, -log10(P),colour="0.6 < r2 < 0.8"), size=ifelse(tmp$POS==zoom2$POS,3,2), shape=ifelse(tmp$POS==zoom2$POS,18,19))}

  tmp=df %>% filter(R2<0.6 & R2>0.4)
  if(length(tmp$POS)> 0){ p1=p1+geom_point(data=tmp, aes(POS, -log10(P),colour="0.4 < r2 < 0.6"), size=ifelse(tmp$POS==zoom2$POS,3,2), shape=ifelse(tmp$POS==zoom2$POS,18,19))}

  tmp=df %>% filter(R2<0.4 & R2>0.2)
  if(length(tmp$POS)> 0){ p1=p1+geom_point(data=tmp, aes(POS, -log10(P),colour="0.2 < r2 < 0.4"), size=ifelse(tmp$POS==zoom2$POS,3,2), shape=ifelse(tmp$POS==zoom2$POS,18,19))}

  tmp=df %>% filter(R2<0.2)
  if(length(tmp$POS)> 0){ p1=p1+geom_point(data=tmp, aes(POS, -log10(P),colour="0 < r2 < 0.2"), size=ifelse(tmp$POS==zoom2$POS,3,2), shape=ifelse(tmp$POS==zoom2$POS,18,19))}


  zoom_tmp=df[FALSE,]
  for(i in ids2label){
    zoom2 = df %>%  filter(ID==i) %>% arrange(P) %>% head(1)
    p1=p1+ggrepel::geom_text_repel(data=unique(zoom2), aes(x=POS,y=-log10(P),label=ID), color="grey40",size=label_size,direction="both",nudge_x = 0.5,nudge_y = 0.5,segment.size=0.2,segment.alpha =.5)

    if("p" %in% colnames(zoom2)){
      zoom2$P=zoom2$p
    }
    zoom_tmp=rbind(zoom_tmp, zoom2) %>% distinct(ID,POS, .keep_all = T)
  }
  p1=p1+ggrepel::geom_text_repel(data=unique(zoom_tmp), aes(x=POS,y=-log10(P),label=ID), color="grey40",size=label_size,direction="both",nudge_x = 0.5,nudge_y = 0.5,segment.size=0.2,segment.alpha =.5)

  if(!is.null(sign_thresh))
    p1=p1+geom_hline(yintercept = -log10(sign_thresh), colour="red", linetype="dashed")



  p1=p1+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  #ßp1=p1+theme(panel.border=element_blank(),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(color="grey"))
  p1=suppressMessages(p1)+theme(axis.text.y=element_text(size=10),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(), title = element_text(size=11))

  p1=suppressMessages(p1+xlim(xmin, xmax)+ ylim(min(-log10(df$P)), max(-log10(df$P)*1.1)))
  p1=p1+scale_color_manual(name="R2", values=colorMap)
  return(p1)
}
