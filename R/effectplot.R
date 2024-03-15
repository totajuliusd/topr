#R

#' Flip to the positive allele for dataset 1
#'
#' @description
#'
#' \code{flip_to_positive_allele_for_dat1()}
#'
#' @param df A dataframe that is in the snpset format (like returned by the get_snpset() function)
#'
#' @return The input dataframe after flipping to the positive effect allele in dataframe 1
#' @export
#'
#' @examples
#' \dontrun{
#' CD_UKBB_index_snps <- get_lead_snps(CD_UKBB)
#' snpset <- get_snpset(CD_UKBB_index_snps, CD_FINNGEN)
#' flip_to_positive_allele_for_dat1(snpset$matched)
#' }
#'

flip_to_positive_allele_for_dat1 <- function(df){
  if(! is.null(df)){
    pos_alleles <- df %>% dplyr::filter(E1 >= 0)
    neg_alleles <- df %>% dplyr::filter(E1 < 0)
    #do an allele and effect swap on negative alleles
    neg_alleles$E1 <- neg_alleles$E1 * (-1)
    neg_alleles$E2 <- neg_alleles$E2 * (-1)
    neg_alleles$ALT1tmp <-  neg_alleles$ALT1
    neg_alleles$ALT1 <- neg_alleles$REF1
    neg_alleles$REF1 <- neg_alleles$ALT1tmp
    neg_alleles$ALT2tmp <-  neg_alleles$ALT2
    neg_alleles$ALT2 <- neg_alleles$REF2
    neg_alleles$REF2 <- neg_alleles$ALT2tmp
    neg_alleles <- neg_alleles %>% dplyr::select(-ALT1tmp, -ALT2tmp)
    return(rbind(pos_alleles, neg_alleles))
  }
  else{
    warning("Input dataframe is empty")
  }
}

#' Match the variants in the snpset by their alleles
#'
#' @description
#'
#' \code{match_by_alleles()}
#'
#' @param df A dataframe that is in the snpset format (like returned by the \code{\link{get_snpset}} function)
#' @param verbose A logical scalar (default: FALSE). Assign to TRUE to get information on which alleles are matched and which are not.
#' @param show_full_output A logical scalar (default:FALSE). Assign to TRUE to show the full output from this function
#'
#' @return The input dataframe containing only those variants with matched alleles in the snpset
#' @export
#'
#' @examples
#' \dontrun{
#' CD_UKBB_lead_snps <- get_lead_snps(CD_UKBB)
#' snpset <- get_snpset(CD_UKBB_lead_snps, CD_FINNGEN)
#' match_by_alleles(snpset$found)
#' }
#'

match_by_alleles <- function(df, verbose=NULL, show_full_output=FALSE){
  if(! is.null(df)){
    df$ID_tmp <- paste(df$CHROM, df$POS, sep="_")
    matched_snps <- df %>% dplyr::filter(REF1 == REF2 & ALT1 == ALT2)
    #check whether ref and alt are reversed
    if(length(matched_snps$POS) != length(df$POS)){
      swapped_allele <-  df %>% dplyr::filter(REF1 == ALT2 & REF2 == ALT1)
      swapped_allele$E2 <- swapped_allele$E2*(-1)
      swapped_allele$REF2tmp <- swapped_allele$REF2
      swapped_allele$REF2 <- swapped_allele$ALT2
      swapped_allele$ALT2 <- swapped_allele$REF2tmp
      swapped_allele <- swapped_allele %>% dplyr::select(-REF2tmp)
      matched_snps <- rbind(matched_snps, swapped_allele)
    }
    #check whether all variants were matched
    not_matched <- df %>% dplyr::filter(! ID_tmp %in% matched_snps$ID_tmp)
    if(! is.null(verbose)){
      if(verbose){
        if(length(not_matched$POS)>0){
          print(paste("Could not match ",length(not_matched$POS), " snps on their alleles."))
          if(length(not_matched$POS)>10 & !show_full_output){
            print("Note! Showing truncated output for snps not matched. For the full output rerun with show_full_output=TRUE")
            print(not_matched[c(1:10),])
            print("...")
          }else{
            print(not_matched)
          }
          print(paste(length(matched_snps$POS), " SNPs have matching REF and ALT alleles. ", sep=""))
        }else{
          print("All SNPs found were matched by their REF and ALT alleles")
        }
      }
    }
    return(list(matched=matched_snps %>% dplyr::select(-ID_tmp), nomatch=not_matched))
  }else{
    warning("The input dataframe is empty")
  }
}

#' Get variants that overlap between two datasets
#'
#' @description
#'
#' \code{match_by_pos()}
#'
#' @param df1 A dataframe of variants, has to contain CHROM and POS
#' @param df2 A dataframe of variants, has to contain CHROM and POS
#' @param verbose A logical scalar (default: FALSE). Assign to TRUE to get information on which alleles are matched and which are not.
#' @param show_full_output A logical scalar (default:FALSE). Assign to TRUE to show the full output from this function
#'
#' @return A list containing two dataframes, one of overlapping snps and the other snps not found in the second input dataset 
#' @export
#'
#' @examples
#' \dontrun{
#' CD_UKBB_index_snps <- get_lead_snps(CD_UKBB)
#' match_by_pos(CD_UKBB_index_snps, CD_FINNGEN)
#' }
#'
match_by_pos <- function(df1, df2, verbose=NULL, show_full_output=FALSE){
  dat2 <- dat_check(df2, verbose=verbose)
  dat1 <- dat_check(df1, verbose=verbose)
  df1 <- dat1[[1]]
  df1_orig <- df1
  df2 <- dat2[[1]]
  snpset <- NULL
  not_found <-NULL
  df2$ID_tmp <- paste(df2$CHROM, df2$POS, sep="_")
  df1$ID_tmp <- paste(df1$CHROM, df1$POS, sep="_")
  if((! "BETA" %in% colnames(df2)) & ("OR" %in% colnames(df2))){
    df2$BETA <- log(df2$OR)
  }
  if((! "BETA" %in% colnames(df1)) & ("OR" %in% colnames(df1))){
    df1$BETA <- log(df1$OR)
  }
  if((! "REF" %in% colnames(df1)) & (! "REF" %in% colnames(df2))){
    warning("The input datasets have to include REF and ALT columns to be able to use this function")
  }else{
    df2 <- df2 %>% dplyr::select(CHROM,POS,REF,ALT,P,BETA) %>% dplyr::rename(P2=P,E2=BETA,ALT2=ALT,REF2=REF)
    if("Gene_Symbol" %in% colnames(df1)){
      df1 <- df1 %>% dplyr::select(CHROM,POS,REF,ALT,P,BETA,Gene_Symbol,ID_tmp) %>% dplyr::rename(P1=P,E1=BETA,ALT1=ALT,REF1=REF)
    }
    else{
      df1 <- df1 %>% dplyr::select(CHROM,POS,REF,ALT,P,BETA, ID_tmp ) %>% dplyr::rename(P1=P,E1=BETA,ALT1=ALT,REF1=REF)
    }
    snpset <- df1 %>% dplyr::inner_join(df2,by=c("CHROM","POS")) %>% dplyr::arrange(CHROM,POS,P1) %>% dplyr::distinct(CHROM,POS, .keep_all = TRUE)
    if("ID" %in% colnames(df1_orig)){
      #add the ID column back to the dataframe
      snpset <- snpset %>% dplyr::inner_join(df1_orig %>% dplyr::select(CHROM,POS,ID), by=c("CHROM","POS")) %>% dplyr::distinct(CHROM,POS, .keep_all = TRUE)
    }
    else{
      print("ID is not in colnames df1")
      print(colnames(df1))
    }
    not_found <- df1 %>% dplyr::filter(! ID_tmp %in% snpset$ID_tmp)
    if(! is.null(verbose)){
      if(verbose){
        print("Overlapping SNPs: ")
        print(paste("There are a total of ",length(df1$POS), " SNPs in the first dataset. Thereof [",length(snpset$POS), "] were are also the second dataset (dat2) and [", length(not_found$POS), "] are not.", sep=""))
        print(paste("SNPs FOUND in dat2: ",length(snpset$POS), sep="" ))
        print(snpset)
        if(length(not_found$POS)>0){
          print(paste("SNPs NOT FOUND in dat2: " ,length(not_found$POS), sep=""))
          if(length(not_found$POS)>10 & !show_full_output){
            print("Note! Showing truncated output for snps not found. For the full output rerun with show_full_output=TRUE")
            print(not_found[c(1:10),])
            print("...")
          }else{
            print(not_found)
          }
        }
      }
    }
    snpset <-snpset %>% dplyr::select(-ID_tmp)
  }
  return(list(found=snpset, not_found=not_found))
}

#' Create a dataframe that can be used as input for making effect plots
#'
#' @description
#'
#' \code{get_snpset()}
#'
#' @param df1 The dataframe to extract the top snps from (with p-value below thresh)
#' @param df2 The dataframe in which to search for overlapping SNPs from dataframe1
#' @param verbose Logical, (default: FALSE). Assign to TRUE to get information on which alleles are matched and which are not.
#' @param show_full_output A logical scalar (default:FALSE). Assign to TRUE to show the full output from this function
#' @param build A string, genome build, choose between builds 37 (GRCh37) and 38 (GRCh38) (default is 38)
#' @inheritParams get_lead_snps
#' 
#' 
#' @return Dataframe of overlapping snps (snpset)
#' @export
#'
#' @examples
#' \dontrun{
#' CD_UKBB_index_snps <-get_lead_snps(CD_UKBB)
#' get_snpset(CD_UKBB_index_snps, CD_FINNGEN)
#' }
#'

get_snpset <- function(df1, df2, thresh=5e-08, protein_coding_only=TRUE, region_size=1000000, verbose=NULL, show_full_output=FALSE, build=38){
  snpset <- NULL
  
  if(is.list(df1) & length(df1)==2){
    if(is.data.frame(df1[[1]]) & is.data.frame(df1[[2]])){
     df2=df1[[2]]
     df1=df1[[1]]
    }
  }
  
  if(! is.numeric(region_size)){
    region_size <- convert_region_size(region_size)
  } 
  if((! "REF" %in% colnames(df1)) & (! "REF" %in% colnames(df2))){
    warning("The input datasets have to include REF and ALT columns to be able to use this function")
  }else{
    snpset_list <- df1 %>% get_lead_snps(thresh = thresh, region_size=region_size, keep_chr=FALSE,verbose=verbose) %>% match_by_pos(df2, verbose=verbose, show_full_output=show_full_output) 
    snps_matched <- snpset_list$found %>% match_by_alleles(verbose=verbose, show_full_output = show_full_output)
    snps_matched_flipped <- snps_matched$matched %>% flip_to_positive_allele_for_dat1()
    if(! "Gene_Symbol" %in% colnames(df1)  ||  ! "gene_symbol" %in% colnames(df1)){
      snps_matched_flipped <- snps_matched_flipped %>% annotate_with_nearest_gene(protein_coding_only = protein_coding_only, build=build)
    }
  }
  return(list(matched=snps_matched_flipped, snp_not_found_in_df2=snpset_list$not_found, snp_found_different_alleles_in_df2=snps_matched$nomatch ))
}

#' Show the code/functions used to get a snpset
#'
#' @description
#'
#' \code{get_snpset_code()}
#'
#' @return Dataframe containing the top hit
#' @export
#'
#' @examples
#' \dontrun{
#' get_snpset_code()
#' }
#'
get_snpset_code <-function(){
  print("snpset_list <- df1 %>% get_lead_snps(thresh = 5e-08,region_size=1000000) %>% match_by_pos(df2)")
  print("snps_matched <- snpset_list$found %>% match_by_alleles()")
  print("snps_matched_flipped %>% snps_matched$matched %>% flip_to_positive_allele_for_dat1()")
}


#' Create a plot comparing variant effects in two datasets
#'
#' @description
#'
#' \code{effectplot()}
#'
#' @param df The input dataframe (snpset) containing one row per variant and P values (P1 and P2) and effects (E1 and E2) from two datasets/phenotypes OR a list containing two datasets.
#' @param pheno_x A string representing the name of the phenotype whose effect is plotted on the x axis
#' @param pheno_y A string representing the name of the phenotype whose effect is plotted on the y axis
#' @param annotate_with  A string, The name of the column that contains the label for the datapoints (default value is Gene_Symbol)
#' @param thresh A number. Threshold cutoff, datapoints with P2 below this threshold are shown as filled circles whereas datapoints with P2 above this threshold are shown as open circles
#' @param ci_thresh A number.Show the confidence intervals if the P-value is below this threshold
#' @param gene_label_thresh Deprecated: A number, label datapoints with P2 below this threshold
#' @param color A string, default value is the first of the topr colors
#' @param subtitle_text_size A number setting the text size of the subtitle (default: 11)
#' @param snpset_thresh A number representing the threshold used to create the snpset used for plotting (Only applicable if the input dataframe is a list containing two datasets)
#' @param snpset_region_size A number representing the region size to use when creating the snpset used for plotting (Only applicable if the input dataframe is a list containing two datasets)
#' @param annotate A number, label datapoints with p-value below below this number (in the second df) by their nearest gene
#' @inheritParams manhattan
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' effectplot(list(CD_UKBB, CD_FINNGEN))
#' }
#'

effectplot <- function(df,pheno_x="x_pheno", pheno_y="y_pheno", annotate_with="Gene_Symbol", thresh=5e-08, ci_thresh=1,gene_label_thresh = 5e-08, 
                       color=get_topr_colors()[1],scale=1, build=38,label_fontface="italic",label_family="",nudge_y=0.001,nudge_x=0.001,size=2,
                       segment.size=0.2,segment.linetype="solid",segment.color="transparent",angle=0,title=NULL,
                       axis_text_size=10,axis_title_size=12, title_text_size=13,subtitle_text_size=11,gene_label_size=3.2,snpset_thresh=5e-08,
                       snpset_region_size=1000000,max.overlaps=10, annotate=0,label_color=NULL){
  if (!missing(gene_label_thresh)) deprecated_argument_msg(gene_label_thresh)
   dat <- df
  if(is.list(dat) & length(dat)==2){
    if(is.data.frame(dat[[1]]) & is.data.frame(dat[[2]])){
       snpset <- get_snpset(df1=dat[[1]],df2=dat[[2]],thresh=snpset_thresh, region_size = snpset_region_size, build=build)
       dat <- snpset$matched
    }
  }
  if(! is.null(dat)){
    if(min(dat$E2) < 0){
      ymin <- abs(min(dat$E2))* (-1.1)
    }else{ymin=0}
    
    ymax <- max(dat$E2) 
    xmax <- max(dat$E1)
    xymax <- ifelse(ymax > xmax, ymax, xmax)
    xymax <- xymax *1.1
    xmin <- 0
    
    dat <- dat%>% dplyr::rename(y_P=P2,y_Effect=E2, x_P=P1, x_Effect=E1)
    P1 <- dat %>% dplyr::filter(y_P < thresh)
    P2 <- dat %>% dplyr::filter(thresh < y_P & y_P < 0.05)
    P3 <- dat %>% dplyr::filter(y_P > 0.05)
    dat$y_P <- signif(dat$y_P,2)
    dat$sigma <- NA

    for (i in seq_len(nrow(dat))){dat$sigma[i] <- -abs(dat$y_Effect[i])/stats::qnorm(dat$y_P[i]/2)}
    dat$C1 <- ifelse(dat$y_P<ci_thresh,dat$y_Effect-1.96*dat$sigma,dat$y_Effect)
    dat$C2 <- ifelse(dat$y_P<ci_thresh,dat$y_Effect+1.96*dat$sigma,dat$y_Effect)
    dat$sigma_x <- NA
    #error bars for phenotype x
    for (i in seq_len(nrow(dat))){dat$sigma_x[i] <- -abs(dat$x_Effect[i])/stats::qnorm(dat$x_P[i]/2)}
    dat$C1_x <- ifelse(dat$x_P<ci_thresh,dat$x_Effect-1.96*dat$sigma_x,dat$x_Effect)
    dat$C2_x <- ifelse(dat$x_P<ci_thresh,dat$x_Effect+1.96*dat$sigma_x,dat$x_Effect)
    if(annotate_with %in% colnames(dat)){
      dat$label <- ifelse(dat$y_P < annotate, dat %>% pull(annotate_with), "")
    } else if("ID" %in% colnames(dat) & annotate_with =="ID"){
      dat$label <- ifelse(dat$y_P < annotate, dat$ID, "")
    } else{
      #annotate the data
      dat_tmp <- dat %>% select(CHROM,POS)
      dat_annot <- annotate_with_nearest_gene(dat_tmp, build=build)
      dat <- merge(dat, dat_annot%>% select(CHROM,POS,Gene_Symbol), by=c("CHROM","POS"))
      dat$label <- ifelse(dat$y_P < annotate, dat$Gene_Symbol, "")
      }
    p1 <- ggplot(data=dat, aes(y=y_Effect,x=x_Effect)) #+geom_point(data=dat[[1]]$gwas, aes(dat[[1]]$gwas$pos_adj, dat[[1]]$gwas$log10p),color=colors[1] alpha=0.7,size=1)+theme_bw()
    
    
    p1 <- p1+geom_errorbar(data=dat,mapping=aes(ymin=C1,ymax=C2),color="grey",alpha=0.4)
    #horizontal errorbar
    p1 <- p1+geom_errorbarh(data=dat,mapping=aes(xmin=C1_x,xmax=C2_x),color="grey",alpha=0.4)
    
    p1 <- p1+geom_point(data=P1, aes(x=x_Effect, y=y_Effect),color=color,shape=19,stroke=1.5,alpha=1,size=size)
    p1 <- p1+geom_point(data=P2, aes(x=x_Effect, y=y_Effect),color=color,shape=1,stroke=1.5, alpha=0.5,size=size)
    p1 <- p1+geom_point(data=P3, aes(x=x_Effect, y=y_Effect),color="#CCCCCC",shape=19,stroke=1.5,alpha=0.1,size=size)
    
    # Remove panel borders and grid lines
    p1 <- p1+theme(legend.background=element_blank(),legend.key=element_rect(fill="white"))
    p1 <- p1+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank())
    p1 <- p1+geom_vline(xintercept=0, colour="grey")+geom_hline(yintercept=0, colour="grey")  
    
    
    p1 <- p1+theme(plot.subtitle=element_text(size=subtitle_text_size*size),plot.title = element_text(size=title_text_size*scale), 
                   axis.title= element_text(size= axis_title_size*scale), axis.text.x = element_text(size=axis_text_size*scale))
    sign_col="black"
    nonsign_col="grey"
    if(! is.null(label_color)){
      sign_col <- label_color; nonsign_col <-label_color
    }
    p1 <- p1+ggrepel::geom_text_repel(aes(label=label),segment.color = segment.color,fontface=label_fontface,family=label_family,angle=angle,
                                      nudge_y=nudge_y,nudge_x=nudge_x,size=gene_label_size*scale,show.legend = FALSE,segment.size=segment.size,segment.linetype=segment.linetype,
                                      colour=ifelse(dat$y_P<thresh,sign_col,nonsign_col),max.overlaps=max.overlaps)
    
    p1 <-p1+geom_abline(slope=1,intercept=0,color="grey44", size=0.2)
    
    if(min(dat$y_Effect) < 0){
      p1 <- p1+geom_abline(slope=-1,intercept=0,color="grey44", size=0.2)
    }
    p1 <- p1+coord_cartesian(xlim=c(xmin,xymax), ylim=c(ymin,xymax))
    if(is.null(title)){
      title <- paste(pheno_x, " vs ", pheno_y, " (N=",length(dat$x_P),")", sep="")
    }
    subtitle <- paste("Filled circle: P < ",thresh, " ; Open circle: ",thresh, " < P < 0.05 ; Grey circle: P > 0.05", sep="")
    p1 <- p1+labs(title=title, y=paste(pheno_y," (Effect)"), x=paste(pheno_x,"(Effect)"), caption=subtitle)
    #set the origin / interception of the y-axis and x-axis to zero
    p1 <- p1 + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    return(p1)
  }
}
