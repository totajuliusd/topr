#' UKBB Crohns disease (ICD 10 code K50)
#'
#' Dataset retrieved from the UK biobank consisting of 2,799 crohn´s cases (K50) and 484,515 controls. The dataset has been filtered on variants with P <1e-03.
#'
#' @format A data frame  with 26,824 rows and 10 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{P}{P-value from Plink run, additive model, reggresion model GLM_FIRTH}
#'   \item{OR}{Odds Ratio}
#' }
#' @source Crohn's UKBB ICD10 code K50, only including variants with P<1e-03
"CD_UKBB"

#' UKBB Ulcerative colitis (ICD 10 code K51)
#'
#' Dataset retrieved from the UK biobank including of 5,452 UC cases (K51) and 481,862 controls. The dataset has been filtered on variants with P<1e-03.
#'
#' @format A data frame  with 57,383 rows and 10 variables
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{P}{P-value from Plink run, additive model, reggresion model GLM_FIRTH}
#'   \item{OR}{Odds Ratio}
#' }
#' @source Ulcerative Colitis UKBB ICD10 code K51, only including variants with P<1e-03
"UC_UKBB"

#' Finngen v5 Crohn‘s disease (CHRONSMALL)
#'
#' Dataset retrieved from the Finngen database (version 5) including 968 crohn´s cases and 210,100 controls. The dataset has been filtered on variants with P <1e-03.
#'
#' @format A data frame  with 29,926rows and 9 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{P}{P-value from Plink run, additive model, reggresion model GLM_FIRTH}
#'   \item{beta}{Variant effect}
#' }
#' @source  Crohn's small intestines (CHRONSMALL), only including variants with P<1e-03
"CD_FINNGEN"

#' Ensembl genes build HG38.104-5-2
#'
#' https://www.ensembl.info/2021/05/05/ensembl-104-has-been-released/
#'
#' genes on chrY and chrM were excluded
#'
#' @format A data frame  with 40,122 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Chromosome on build version 38 (GRCh38/hg38)}
#'   \item{gene_start}{genetic position of gene start on build version 38}
#'   \item{gene_end}{genetic position of gene end on build version 38 }
#'   \item{gene_symbol}{The name of the gene}
#'   \item{biotype}{the biotype of the gene}
#'   }
"ENSGENES"


#' Ensembl exons build HG38-104-5-2
#'
#' https://www.ensembl.info/2021/05/05/ensembl-104-has-been-released/
#'
#' exons on chrY and chrM were excluded from the exon dataset
#'
#' @format A data frame  with 40,122 rows and 7 variables:
#' \describe{
#'   \item{chrom}{Chromosome on build version 38 (GRCh38/hg38)}
#'   \item{gene_start}{genetic position of gene start on build version 38}
#'   \item{gene_end}{genetic position of gene end on build version 38 }
#'   \item{gene_symbol}{The name of the gene}
#'   \item{exon_chromstart}{genetic positions of exon start}
#'   \item{exon_chromend}{genetic position of exon end}
#' }
"ENSEXONS"

#' Example dataset including the R2 column for the locuszoom plot function
#'
#' The dataset is a subset of CD_UKBB and only includes variants above and near the IL23R gene on chromosome 1
#'
#' @format A data frame  with 26,824 rows and 10 variables:
#' \describe{
#'   \item{CHROM}{Chromosome, written as for example chr1 or 1}
#'   \item{POS}{genetic position of the variant}
#'   \item{ID}{Variant identifier, e.g. rsid }
#'   \item{P}{P-value from Plink run, additive model, reggresion model GLM_FIRTH}
#'   \item{R2}{variant correlation (r^2)}
#' }
#' @source A subset of the CD_UKBB dataset
"R2_CD_UKBB"
