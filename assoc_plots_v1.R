## Plots
#
library(gorr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library("egg")
library("grid")

#chr_lengths_file="~/Workflows/lib/build38_chr_lengths.tsv"

chr_lengths_file="/Users/thorhildur/IdeaProjects/notebooks-workflows/lib/build38_chr_lengths.tsv"

#colors

#and plot using egg
#https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

source("/Users/thorhildur/IdeaProjects/topR/R/chromplot.R")
source("/Users/thorhildur/IdeaProjects/topR/R/setters_and_getters.R")
source("/Users/thorhildur/IdeaProjects/topR/R/geneplot.R")
source("/Users/thorhildur/IdeaProjects/topR/R/exonplot.R")
source("/Users/thorhildur/IdeaProjects/topR/R/help_functions.R")
source("/Users/thorhildur/IdeaProjects/topR/R/locuszoomplot.R")
source("/Users/thorhildur/IdeaProjects/topR/R/regionplot.R")
source("/Users/thorhildur/IdeaProjects/topR/R/dat_checks.R")
source("/Users/thorhildur/IdeaProjects/topR/R/miscellaneous.R")
source("/Users/thorhildur/IdeaProjects/topR/R/overviewplot.R")
source("/Users/thorhildur/IdeaProjects/topR/R/manhattan.R")


#darkblue, orange(E69F00), turqoise(00AFBB), grey (999999), red ("#FC4E07")
color=c("darkblue","#E69F00","#00AFBB","#999999","#FC4E07","darkorange1")

manhattan

