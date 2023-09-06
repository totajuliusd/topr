filter_on_chr <- function(dat,chr){
  chr <- gsub("chr", "", chr)
  for(i in seq_along(dat)){
    dat[[i]] <- dat[[i]] %>% dplyr::filter(CHROM==chr)
  }
  return(dat)
}

filter_on_xmin_xmax <- function(dat,xmin=NULL,xmax=NULL){
  for(i in seq_along(dat)){
    if(! is.null(xmin)) {
      if(xmin < 0){
        xmin <- 0
      }
      dat[[i]] <- dat[[i]] %>% dplyr::filter(POS >= xmin)
    }
    if(! is.null(xmax)) dat[[i]] <-  dat[[i]] %>% dplyr::filter(POS <= xmax)
  }
  return(dat)
}


include_chrX <- function(dat){
  for(i in seq_along(dat)){
    if("23" %in% unique(dat[[i]]$CHROM)){
      return(TRUE)
    }
  }
  return(FALSE)
}

is_df_empty <- function(df, type){
  if(is.null(type)) type <- "main data "
  if(is.data.frame(df)) df <- list(df)
  for(i in seq_along(df)){
    datrows <- nrow(df[[i]])
    if(datrows == 0){
      stop(paste("One or more of the dataframes provided as input is empty!!! (", type," n=",i,")",sep=""))
    }
  }
}

rename_value <- function(x, value) {
  if (length(x) == 0L) {
    character()
  } else {
    value
  }
}



#' @importFrom grDevices col2rgb

lightness = function(col, light=0.5) {
  
  if (length(col)==1) {
    if (length(light)>1) {
      tmp = col
      col = light
      col[] = tmp
    }
  } else {
    if (length(light)==1) {
      light = rep(light,length(col))
    } else {
      if (length(col)!=length(light)) stop('col and light must have the same length (or be of length 1).')
    }
  }
  
  for (i in seq(length(col))) {
    rgb = col2rgb(col[i],alpha=T)/255
    if (light[i]<0.5) {
      rgb[1:3] = rgb[1:3]*light[i]*2
    } else {
      rgb[1:3] = 1-(1-rgb[1:3])*(1-light[i])*2
    }
    if (rgb[4]==1) {
      col[i] = rgb(rgb[1],rgb[2],rgb[3])
    } else {
      col[i] = rgb(rgb[1],rgb[2],rgb[3],rgb[4])
    }
  }
  
  return(col)
  
}

