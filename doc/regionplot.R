## ---- include=FALSE-----------------------------------------------------------
    library(topr)
    library(dplyr)

## ----regionplot, eval=FALSE---------------------------------------------------
#  regionplot(CD_UKBB,
#             gene="IL23R")

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/regionplot.jpg')

## ---- eval=FALSE--------------------------------------------------------------
#  regionplot(CD_UKBB,
#             gene="IL23R",
#             annotate=5e-9)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/regionplot_annotate.jpg')

## ----regionplot_annotate, eval=FALSE------------------------------------------
#  regionplot(CD_UKBB,
#             gene = "IL23R",
#             annotate_with_vline = 5e-9,
#             region_size = 100000)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/regionplot_annotate_vline.jpg')

## ----regionplot_multi, eval=FALSE---------------------------------------------
#  regionplot(list(UC_UKBB,CD_UKBB,CD_FINNGEN),
#             gene = "IL23R",
#             annotate_with_vline = 5e-6,
#             legend_labels = c("UC UKBB","CD UKBB","CD FINNGEN"))

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/regionplot_multi.jpg')

## ----locuszoom, eval=FALSE----------------------------------------------------
#  locuszoom(R2_CD_UKBB)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/locuszoom.jpg')

## ----locuszoom_annotate, eval=FALSE-------------------------------------------
#  locuszoom(R2_CD_UKBB,
#            annotate_with_vline = 1e-9,
#            region_size = 100000)

## ----echo=FALSE, out.width='100%'---------------------------------------------
knitr::include_graphics('figures/locuszoom_annotate.jpg')

