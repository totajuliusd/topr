#' UKBB Crohns disease (ICD 10 code K50)
#'
#' Dataset retrieved from the UK biobank consisting of 2,799 crohn´s cases (K50) and 484,515 controls. The dataset has been filtered on variants with P <1e-03.
#'
#' @format A data frame  with 26,824 rows and 10 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{REF}{Reference allele}
#'   \item{ALT}{Alternative allele}
#'   \item{SE}{Standard Error}
#'   \item{P}{P-value from Plink run, additive model, reggresion model GLM_FIRTH}
#'   \item{OR}{Odds Ratio}
#'   \item{Max_Impact}{Max variant impact (HIGH,MODERATE,LOW or LOWEST)}
#'   \item{cohort_AF}{variant allele frequency in the UKBB cohort)}
#' }
#' @source Crohn's UKBB ICD10 code K50 snps with pvalues below p1e-03
"CD_UKBB"

#' UKBB Ulcerative colitis (ICD 10 code K51)
#'
#' Dataset retrieved from the UK biobank including of 5,452 UC cases (K51) and 481,862 controls. The dataset has been filtered on variants with P <1e-03.
#'
#' @format A data frame  with 57,383 rows and 10 variables
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{REF}{Reference allele}
#'   \item{ALT}{Alternative allele}
#'   \item{SE}{Standard Error}
#'   \item{P}{P-value from Plink run, additive model, reggresion model GLM_FIRTH}
#'   \item{OR}{Odds Ratio}
#'   \item{Max_Impact}{Max variant impact (HIGH,MODERATE,LOW or LOWEST)}
#'   \item{cohort_AF}{variant allele frequency in the UKBB cohort)}
#' }
#' @source UC_K51_UKBB_p1e-03_5452_cases_481862_ctrls
"UC_UKBB"

#' Finngen v5 Crohn‘s disease (CHRONSMALL)
#'
#' Dataset retrieved from the Finngen database (version 5) including 968 crohn´s cases and 210,300 controls. The dataset has been filtered on variants with P <1e-03.
#'
#' @format A data frame  with 29,926rows and 9 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{REF}{Reference allele}
#'   \item{ALT}{Alternative allele}
#'   \item{SE}{Standard Error}
#'   \item{P}{P-value from Plink run, additive model, reggresion model GLM_FIRTH}
#'   \item{beta}{Variant effect}
#'   \item{cohort_AF}{variant allele frequency in the UKBB cohort)}
#' }
#' @source CHRONSMALL_968_cases_210300_controls_p_1e-03
"CD_FINNGEN"

#' Ensembl genes build 38
#'
#' @format A data frame  with 60,564 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Chromosome on build version 38 (GRCh38/hg38)}
#'   \item{gene_start}{genetic position of gene start on build version 38}
#'   \item{gene_end}{genetic position of gene end on build version 38 }
#'   \item{gene_symbol}{The name of the gene}
#'   \item{biotype}{the biotype of the gene}
#'   }
"ENSGENES"


#' Ensembl exons build 38
#'
#' chrY and chrM were excluded from the exon dataset
#'
#' @format A data frame  with 60,006 rows and 7 variables:
#' \describe{
#'   \item{chrom}{Chromosome on build version 38 (GRCh38/hg38)}
#'   \item{gene_start}{genetic position of gene start on build version 38}
#'   \item{gene_end}{genetic position of gene end on build version 38 }
#'   \item{gene_symbol}{The name of the gene}
#'   \item{exon_chromstart}{genetic positions of exon start}
#'   \item{exon_chromend}{genetic position of exon end}
#' }
"ENSEXONS"
