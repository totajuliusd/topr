## ----include=FALSE------------------------------------------------------------
library(topr)
library(dplyr)
library(knitr)

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(CD_UKBB)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(CD_UKBB,
#            annotate = 5e-09,
#            title = "Crohn's disease")

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_annotate.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  genes = c("IL23R","NOTCH4","NOD2","JAK2","TTC33")
#  manhattan(CD_UKBB,
#            annotate = 5e-09,
#            title = "Crohn's disease",
#            highlight_genes = genes)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_genes.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(CD_UKBB,
#            annotate = 5e-09,
#            chr = "chr1")

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_chr1.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(list(UC_UKBB, CD_UKBB),
#            legend_labels = c("UC UKBB", "CD UKBB"))

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_multi.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(list(CD_UKBB, CD_FINNGEN, UC_UKBB),
#            legend_labels = c("CD UKBB", "CD FinnGen","UC UKBB"),
#            annotate = c(5e-9,5e-12,1e-15),
#            region_size = "3Mb",
#            ntop = 2,
#            highlight_genes = genes,
#            highlight_genes_ypos = -0.5 ,
#            title = "Inflammatory Bowel Disease")

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_multi_ntop.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(list(CD_UKBB, CD_FINNGEN, UC_UKBB),
#            legend_labels = c("CD UKBB", "CD FinnGen","UC UKBB"),
#            annotate = c(5e-9,5e-12,1e-15),
#            region_size = "3Mb",
#            ntop = 2,
#            highlight_genes = genes,
#            highlight_genes_ypos = -1.5,
#            title = "Inflammatory Bowel Disease",
#            ymax = 65,
#            ymin = -55,
#            nudge_y = 12,
#            angle = 90)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_multi_ntop_tidy.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(list(CD_UKBB, CD_FINNGEN, UC_UKBB),
#            legend_labels = c("CD UKBB", "CD FinnGen","UC UKBB"),
#            annotate = c(5e-9,5e-12,1e-15),
#            region_size = "3Mb",
#            ntop = 2,
#            highlight_genes = genes,
#            highlight_genes_ypos = -1.5 ,
#            title = "Inflammatory Bowel Disease",
#            ymax = 65,
#            ymin = -55,
#            nudge_y = 12,
#            angle = 90,
#            theme_grey = T)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_multi_ntop_tidy_grey.jpg')

## ----eval=FALSE---------------------------------------------------------------
#  manhattan(list(CD_UKBB, CD_FINNGEN, UC_UKBB),
#            legend_labels = c("CD UKBB", "CD FinnGen","UC UKBB"),
#            annotate_with_vline  = c(5e-9,5e-100,1e-100),
#            region_size="3Mb",
#            ntop=1,
#            highlight_genes = genes,
#            highlight_genes_ypos = -0.5,
#            title = "Inflammatory Bowel Disease",
#            nudge_y = 6,
#            angle = 90,
#            alpha = c(1, 1, 0.7),
#            size = c(1, 1.5, 0.9),
#            shape = c(19,6,19),
#            ymax=33)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/manhattan_multi_ntop_tidy_vline.jpg')

