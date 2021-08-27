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
#' flip_to_positive_allele_for_dat1(df)
#' }
#'

flip_to_positive_allele_for_dat1 <- function(df){
  pos_alleles <- df %>% filter(E1 >= 0)
  neg_alleles <- df %>% filter(E1 < 0)
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

#' Match the variants in the snpset by their alleles
#'
#' @description
#'
#' \code{matched_alleles()}
#'
#' @param df A dataframe that is in the snpset format (like returned by the get_snpset() function)
#'
#' @return The input dataframe containing only those variants whith matched alleles in the snpset
#' @export
#'
#' @examples
#' \dontrun{
#' match_alleles(df)
#' }
#'

match_alleles <- function(df){
  df$ID_tmp <- paste(df$CHROM, df$POS, sep="_")
  matched_snps <- df %>% dplyr::filter(REF1 == REF2 & ALT1 == ALT2)
  #check whether ref and alt are reversed
  if(length(matched_snps$POS) != length(df$POS)){
    swapped_allele <-  snpset1 %>% dplyr::filter(REF1 == ALT2 & REF2 == ALT1)
    swapped_allele$E2 <- swapped_allele$E2*(-1)
    swapped_allele$REF2tmp <- swapped_allele$REF2
    swapped_allele$REF2 <- swapped_allele$ALT2
    swapped_allele$ALT2 <- swapped_allele$REF2tmp
    swapped_allele <- swapped_allele %>% dplyr::select(-REF2tmp)
    matched_snps <- rbind(matched_snps, swapped_allele)
  }
  #check whether all variants were mached
  not_matched <- df %>% dplyr::filter(! ID_tmp %in% matched_snps$ID_tmp)
  if(length(not_matched$POS)>0){
    print(paste("Could not match ",length(not_matched$POS), " snps on their alleles."))
    print(not_matched)
    print(paste(length(matched_snps$POS), " SNPs have matching REF and ALT alleles. ", sep=""))
  }else{
    print("All SNPs were matched by their REF and ALT alleles")
  }
  return(matched_snps %>% dplyr::select(-ID_tmp))
}

#' Get variants that overlap between two datasets
#'
#' @description
#'
#' \code{get_overlapping_snps_by_pos()}
#'
#' @param df1 A dataframe of variants, has to contain CHROM and POS
#' @param df2 A dataframe of variants, has to contain CHROM and POS
#'
#' @return The input dataframe containing only those variants whith matched alleles in the snpset
#' @export
#'
#' @examples
#' \dontrun{
#' get_overlapping_snps_by_pos(dat1, dat2)
#' }
#'
get_overlapping_snps_by_pos <- function(df1, df2){
  dat2 <- dat_check(df2)
  dat1 <- dat_check(df1)
  df1 <- dat1[[1]]
  df1_orig <- df1
  df2 <- dat2[[1]]
  df2$ID_tmp <- paste(df2$CHROM, df2$POS, sep="_")
  df1$ID_tmp <- paste(df1$CHROM, df1$POS, sep="_")
  if((! "BETA" %in% colnames(df2)) & ("OR" %in% colnames(df2))){
    df2$BETA <- log(df2$OR)
  }
  if((! "BETA" %in% colnames(df1)) & ("OR" %in% colnames(df1))){
    df1$BETA <- log(df1$OR)
  }
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
  print("Overlapping SNPs: ")
  not_found <- df1 %>% dplyr::filter(! ID_tmp %in% snpset$ID_tmp)
  print(paste("There are a total of ",length(df1$POS), " SNPs in the first dataset. Thereof [",length(snpset$POS), "] were are also the second dataset (dat2) and [", length(not_found$POS), "] are not.", sep=""))
  print(paste("SNPs FOUND in dat2: ",length(snpset$POS), sep="" ))
  print(snpset)
  if(length(not_found$POS)>0){
    print(paste("SNPs NOT FOUND in dat2: " ,length(not_found$POS), sep=""))
    print(not_found)
  }
  return(snpset %>% dplyr::select(-ID_tmp))
}

#' Create a dataframe that can be used as input for making effect plots
#'
#' @description
#'
#' \code{create_snpset()}
#' All other input parameters are optional
#'
#' @param df1 The dataframe to extract the top snps from (with p-value below thresh)
#' @param df2 The dataframe in which to search for overlapping SNPs from dataframe1
#' @param thresh The p-value threshold used for extracting the top snps from dataset 1
#' @param region The size of the interval which to extract the top snps from
#' @param protein_coding_only Set this variable to TRUE to only use protein_coding genes for the annotation

#' @return Dataframe containing the top hit
#' @export
#'
#' @examples
#' \dontrun{
#' create_snpset(CD_UKBB,CD_FINNGEN, thresh=1e-09)
#' }
#'

create_snpset <- function(df1, df2, thresh=1e-08,protein_coding_only=TRUE, region=1000000){
  snpset <- df1 %>% get_best_snp_per_MB(thresh = thresh,region=region) %>% get_overlapping_snps_by_pos(df2) %>% match_alleles() %>% flip_to_positive_allele_for_dat1()
  if(! "Gene_Symbol" %in% colnames(df1)  ||  ! "gene_symbol" %in% colnames(df1)){
     snpset <- snpset %>% annotate_with_nearest_gene(protein_coding_only = protein_coding_only)
  }
  return(snpset)
}

#' Show the code/functions used to create a snpset
#'
#' @description
#'
#' \code{create_snpset_code()}
#'
#' @return Dataframe containing the top hit
#' @export
#'
#' @examples
#' \dontrun{
#' create_snpset_code()
#' }
#'
create_snpset_code <-function(){
  print("snpset <- df1 %>% get_best_snp_per_MB(thresh = 1e-08,region=1000000) %>% get_overlapping_snps_by_pos(df2) %>% match_alleles() %>% flip_to_positive_allele_for_dat1() %>% annotate_with_nearest_gene(protein_coding_only = TRUE)")
}


#' Create a plot comparing effects within two datasets
#'
#' @description
#'
#' \code{effect_plot()}
#'
#' @param dat The input dataframe (snpset) containing one row per variant and P values (P1 and P2) and effects (E1 and E2) from two datasets/phenotypes
#' @param pheno_x The name of the phenotype whose effect is plotted on the x axis
#' @param pheno_y The name of the phenotype whose effect is plotted on the y axis
#' @param annotate_with  The name of the column that contains the label for the datapoints (default value is Gene_Symbol)
#' @param thresh Threshold cutoff, datapoints with P2 below this threshold are shown as filled circles whereas datapoints with P2 above this threshold are shown as open circles
#' @param ci_thresh Show the confidence intervals if the P-value is below this threshold
#' @param gene_label_thresh Label datapoints with P2 below this threshold
#' @param color default value is the first of the topr colors
#' @export
#'
#' @examples
#' \dontrun{
#' effect_plot(dat)
#' }
#'

effect_plot <- function(dat,pheno_x="x_pheno", pheno_y="y_pheno", annotate_with="Gene_Symbol", thresh=5E-07, ci_thresh=1,gene_label_thresh = 1e-04, color=get_topr_colors()[1]){
  ymax <- max(dat$E2)*1.1
  ymin <- abs(min(dat$E2))* (-1.1)
  xmax <- max(dat$E1)*1.1
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
    dat$label <- ifelse(dat$y_P < gene_label_thresh, dat %>% pull(annotate_with), "")
  } else if("ID" %in% colnames(dat)){
    print("Could not find Gene_Symbol in the input data, so using ID to label the datapoints instead")
     dat$label <- ifelse(dat$y_P < gene_label_thresh, dat$ID, "")
  } else{
    dat$label <- ifelse(dat$y_P < gene_label_thresh, paste(dat$CHROM, dat$POS, sep="_"), "")
    print("Could not find Gene_Symbol nor ID in the input data, so using CHROM and POS to label the datapoints instead")
  }

  p1 <- ggplot(data=dat, aes(y=y_Effect,x=x_Effect)) #+geom_point(data=dat[[1]]$gwas, aes(dat[[1]]$gwas$pos_adj, dat[[1]]$gwas$log10p),color=colors[1] alpha=0.7,size=1)+theme_bw()
  p1 <- p1+geom_errorbar(data=dat,mapping=aes(ymin=C1,ymax=C2),width=0.002,color="grey",alpha=0.4)
  #horizontal errorbar
  p1 <- p1+geom_errorbarh(data=dat,mapping=aes(xmin=C1_x,xmax=C2_x),color="grey",alpha=0.4)
  p1 <- p1+geom_point(data=P1, aes(x=x_Effect, y=y_Effect),color=color,shape=19,stroke=1.5,alpha=1,size=2)
  p1 <- p1+geom_point(data=P2, aes(x=x_Effect, y=y_Effect),color=color,shape=1,stroke=1.5, alpha=0.5,size=2)
  p1 <- p1+geom_point(data=P3, aes(x=x_Effect, y=y_Effect),color="#CCCCCC",shape=19,stroke=1.5,alpha=0.1,size=2)
  p1 <- p1+theme(legend.background=element_blank(),legend.key=element_rect(fill="white"))
  p1 <- p1+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_blank())
  p1 <- p1+theme(plot.subtitle=element_text(size=12),
              axis.title.x= element_text(size=12), axis.title.y = element_text(size=12),
              legend.text = element_text(size=11))
  p1 <- p1+geom_vline(xintercept=0, colour="grey")+geom_hline(yintercept=0, colour="grey")    #+guides(colour=FALSE)
  p1 <- p1+geom_vline(xintercept=0, colour="grey")+guides(colour=FALSE)+geom_hline(yintercept=0, colour="grey")
  p1 <- p1+ggrepel::geom_text_repel(aes(label=label),segment.color = "transparent",fontface="italic",nudge_y=0.001,nudge_x=0.001,size=3.7,show.legend = FALSE,colour=ifelse(dat$y_P<gene_label_thresh,"black","grey"))
  #p1=p1+geom_smooth(method=lm, aes(weight=W),se=FALSE, colour="black", linetype="longdash", size=0.2, alpha=0.4,fullrange = TRUE)
  p1 <-p1+geom_abline(slope=1,intercept=0,color="grey44", size=0.2)

  if(min(dat$y_Effect) < 0){
    p1 <- p1+geom_abline(slope=-1,intercept=0,color="grey44", size=0.2)
  }
  p1 <- p1+coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  maintitle <- paste(pheno_x, " vs ", pheno_y, " (N=",length(dat$x_P),")", sep="");
  subtitle <- paste("Filled circle: P < ",thresh, " ; Open circle: ",thresh, " < P < 0.05 ; Grey circle: P > 0.05", sep="")
  p1 <- p1+labs(title=maintitle, y=paste(pheno_y," (Effect)"), x=paste(pheno_x,"(Effect)"), caption=subtitle)
               #fill=paste("Pval for pheno y", sep="")) #caption=subTitle
  return(p1)
}

