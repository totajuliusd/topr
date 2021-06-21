
mk_overview_plot=function(dat,xmin=0,xmax=0,sign_thresh=1e-09, sign_thresh_color="red",chr=NULL){
  if(is.null(chr))
    chr=get_chr_from_df(dat[[1]])
  df=dat[[1]]
  p1=ggplot(df)+geom_point(data=df, aes(x=POS, y=-log10(P)), color=df$color, size=0.3, shape=df$shape)+theme_bw()

  if(length(dat)>1){
    for(i in 2:length(dat)){
      df=dat[[i]]
      p1=p1+geom_point(data=df, aes(x=POS, y=-log10(P)), shape=df$shape,size=0.3,alpha=df$alpha,color=df$color) # ,shape=df$shape, color=df$color, size=df$size)
    }
  }

  #add the rectangle to the plot
  if(!is.null(xmax) && ! is.null(xmin)){
    #add the overview rectangle
    ymin=min(-log10(dat[[1]]$P))
    ymax=max(-log10(dat[[1]]$P))
    if(length(dat)>1){
      for(i in 2:length(dat)){
        if( min(-log10(dat[[i]]$P)) < ymin)
          ymin=min(-log10(dat[[i]]$P))
        if( max(-log10(dat[[i]]$P)) > ymax)
          ymax=max(-log10(dat[[i]]$P))
      }
    }
    p1=p1+geom_rect(mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), size=0.2, alpha=0.3, color="red", fill=NA)
  }

  p1=add_sign_thresh_to_plot(p1,df, sign_thresh, sign_thresh_color,min(df$POS),max(df$POS),ymin,ymax, overviewplot = TRUE)

  #remove space between axis and plots
  p1=tidy_plot(p1)
  #dont label the axis or add any ticks to the overviewplot
  p1=p1+theme(axis.title = element_blank(),
              axis.text=element_blank(),axis.ticks = element_blank())
  return(p1)

}