#' UKBB Crohns disease (ICD 10 code K50)
#'
#' Dataset retrieved from the UK biobank consisting of 2,799 crohn´s cases (K50) and 484,515 controls. The dataset has been filtered on variants with P <1e-03.
#'
#' @format A data frame  with 21,717 rows and 8 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{REF}{the reference allele}
#'   \item{ALT}{the alternative allele}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{P}{P-value from Plink run, additive model, regression model GLM_FIRTH}
#'   \item{OR}{Odds Ratio}
#'   \item{AF}{Allele frequency}
#' }
#' @source Crohn's UKBB ICD10 code K50, only including variants with P<1e-03
"CD_UKBB"

#' UKBB Ulcerative colitis (ICD 10 code K51)
#'
#' Dataset retrieved from the UK biobank including of 5,452 UC cases (K51) and 481,862 controls. The dataset has been filtered on variants with P<1e-03.
#'
#' @format A data frame  with 45,012 rows and 8 variables
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'  \item{REF}{the reference allele}
#'   \item{ALT}{the alternative allele}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{P}{P-value from Plink run, additive model, regression model GLM_FIRTH}
#'   \item{OR}{Odds Ratio}
#'   \item{AF}{Allele frequency}
#' }
#' @source Ulcerative Colitis UKBB ICD10 code K51, only including variants with P<1e-03
"UC_UKBB"

#' Finngen r7 Crohn‘s disease (K11_CROHNS)
#'
#' Dataset retrieved from the Finngen database (version 7) including 3147 crohn´s cases (K50) and 296,100 controls. The dataset has been filtered on variants with P <1e-03.
#' FinnGen data are publicly available and were downloaded from https://finngen.fi.
#'
#' @format A data frame  with 32,303 rows and 8 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{REF}{the reference allele}
#'   \item{ALT}{the alternative allele}
#'   \item{P}{P-value from Plink run, additive model, regression model GLM_FIRTH}
#'   \item{BETA}{Variant effect}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{AF}{Allele frequency}
#' }
#' @source  Crohn's K50 (K11_CROHNS), only including variants with P<1e-03
"CD_FINNGEN"



#' Example dataset including the R2 column for the locuszoom plot function
#'
#' The dataset is a subset of CD_UKBB and only includes variants above and near the IL23R gene on chromosome 1
#'
#' @format A data frame  with 329 rows and 5 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{P}{P-value from Plink run, additive model, regression model GLM_FIRTH}
#'   \item{R2}{variant correlation (r^2)}
#' }
#' @source A subset of the CD_UKBB dataset
"R2_CD_UKBB"
