# Title     : TODO
# Objective : TODO
# Created by: thorhildur
# Created on: 27/04/2021

#' QQ Plot
#'
#' @param df dataframe with required data
#' @param scale plot elements scale, default: 1
#' @param n_variants number of total variants used in the study
#' @param fontfamily font
#'
#' @return ggplot
#' @export
#' @import ggplot2
#' @import showtext
#'
#'
#'
# Title     : TODO
# Objective : TODO
# Created by: thorhildur
# Created on: 27/04/2021

#' QQ Plot
#'
#' @param df dataframe with required data
#' @param scale plot elements scale, default: 1
#' @param n_variants number of total variants used in the study
#' @param fontfamily font
#'
#' @return ggplot
#' @export
#' @import ggplot2
#' @import showtext
#'

qqplot <- function(df, scale = 1, n_variants = 0, fontfamily = "", breaks = 15) {
  bonf_thres <- if (n_variants > 0) -log10(0.05 / n_variants) else 0

  if (!"exp" %in% colnames(df)) {
    log_line("Expected p-values column (exp) not found in data, creating it")
    if (n_variants == 0) {
      stop("Missing parameter for total variants, which is required for computing expected p-values")
    }

    df <- df %>%
      dplyr::arrange(pval) %>%
      dplyr::mutate(rank = rank(pval, ties.method = "first"),
                    exp = rank/n_variants)
  }

  if (class(df$pval) == "character") {
    df$pval <- readr::parse_number(df$pval)
  }
  if (class(df$exp) == "character") {
    df$exp <- readr::parse_number(df$exp)
  }

  zero_pvals <- dplyr::filter(df, pval <= .Machine$double.xmin )
  caption <- ""
  if (nrow(zero_pvals) > 0) {
    exclusion_list <- with(zero_pvals, paste(chrom, pos, ref, alt, sep = ":", collapse = ", "))
    exclusion_list <- stringr::str_trunc(exclusion_list, width = 100)
    caption <- stringr::str_glue("{nrow(zero_pvals)} zero-value p-value{if(nrow(zero_pvals) > 1) 's' else ''} removed from plot: {exclusion_list}")
  }

  df <- df %>%
    dplyr::filter(pval > .Machine$double.xmin) %>%
    dplyr::mutate(pval = -log10(pval),
                  exp = -log10(exp))

  min_theoretical <- min(df$exp)

  df=df %>% arrange(cohort_AF)
  num_groups=4
  df_by_AF=df %>% group_by((row_number()-1) %/% (n()/num_groups)) %>%nest %>% pull(data)


  ggplot2::ggplot(df_by_AF[[1]],aes(exp,pval)) +
    # geom_hline(aes(yintercept = val), color = "firebrick3", lty = 2, data = tibble(val = bonf_thres)) +
    ggplot2::geom_hline(aes(yintercept = val), color = "#606060", lty = 3, size = .5, data = dplyr::tibble(val = min_theoretical)) +
    ggplot2::geom_vline(aes(xintercept = val), color = "#606060", lty = 3, size = .5, data = dplyr::tibble(val = min_theoretical)) +
    #annotate(geom = "text", label = "Bonferroni threshold",
    #         x = min(df$theoretical), y = bonf_thres, family = fontname, vjust = -.5, size = 4*scale,
    #        hjust = 0, color = "firebrick") +
    ggplot2::geom_smooth(aes(y = exp), method = "lm", se = F, color = "#808080", size = .5, fullrange = T) +
    ggplot2::geom_point(size = scale,color="red")+

    #  p1=p1+geom_point(data=variants[[i]], aes(x=pos_adj, y=-log10(P),color=color), size=variants[[i]]$size, shape=variants[[i]]$shape)
    #  ggplot2::geom_point(df_by_AF[[1]], mapping=aes(x=exp,y=pval, size = scale,color="red"))+

    ggplot2::geom_point(df_by_AF[[2]],aes(x=exp,y=pval,color="blue"),size=scale) +
    #   ggplot2::geom_point(df_by_AF[[3]],mapping=aes(x=exp,y=pval,size = scale,color="green")) +
    #  ggplot2::geom_point(df_by_AF[[4]],mapping=aes(x=exp,y=pval, size = scale,color="oragne")) +

    #ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(breaks), labels = log_labs, limits = c(0,NA) ) +
    #ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(breaks), labels = log_labs, limits = c(0,NA) ) +
    #ggplot2::scale_color_brewer(palette = "red") +
    ggplot2::labs(title = "", x = "Theoretical", y = "Observed", color = "Test", caption = caption) +
    ggplot2::theme_minimal(18*scale, base_family = fontfamily) +
    ggplot2::theme(legend.position = "right", legend.direction = "vertical")
}


qqplot_original <- function(df, scale = 1, n_variants = 0, fontfamily = "", breaks = 15) {
  bonf_thres <- if (n_variants > 0) -log10(0.05 / n_variants) else 0

  if (!"exp" %in% colnames(df)) {
    log_line("Expected p-values column (exp) not found in data, creating it")
    if (n_variants == 0) {
      stop("Missing parameter for total variants, which is required for computing expected p-values")
    }

    df <- df %>%
      dplyr::arrange(pval) %>%
      dplyr::mutate(rank = rank(pval, ties.method = "first"),
                    exp = rank/n_variants)
  }

  if (class(df$pval) == "character") {
    df$pval <- readr::parse_number(df$pval)
  }
  if (class(df$exp) == "character") {
    df$exp <- readr::parse_number(df$exp)
  }

  zero_pvals <- dplyr::filter(df, pval <= .Machine$double.xmin )
  caption <- ""
  if (nrow(zero_pvals) > 0) {
    exclusion_list <- with(zero_pvals, paste(chrom, pos, ref, alt, sep = ":", collapse = ", "))
    exclusion_list <- stringr::str_trunc(exclusion_list, width = 100)
    caption <- stringr::str_glue("{nrow(zero_pvals)} zero-value p-value{if(nrow(zero_pvals) > 1) 's' else ''} removed from plot: {exclusion_list}")
  }

  df <- df %>%
    dplyr::filter(pval > .Machine$double.xmin) %>%
    dplyr::mutate(pval = -log10(pval),
                  exp = -log10(exp))

  min_theoretical <- min(df$exp)

  ggplot2::ggplot(df, aes(exp, pval)) +
    # geom_hline(aes(yintercept = val), color = "firebrick3", lty = 2, data = tibble(val = bonf_thres)) +
    ggplot2::geom_hline(aes(yintercept = val), color = "#606060", lty = 3, size = .5, data = dplyr::tibble(val = min_theoretical)) +
    ggplot2::geom_vline(aes(xintercept = val), color = "#606060", lty = 3, size = .5, data = dplyr::tibble(val = min_theoretical)) +
    #annotate(geom = "text", label = "Bonferroni threshold",
    #         x = min(df$theoretical), y = bonf_thres, family = fontname, vjust = -.5, size = 4*scale,
    #        hjust = 0, color = "firebrick") +
    ggplot2::geom_smooth(aes(y = exp), method = "lm", se = F, color = "#808080", size = .5, fullrange = T) +
    ggplot2::geom_point(size = scale) +
    #ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(breaks), labels = log_labs, limits = c(0,NA) ) +
    #ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(breaks), labels = log_labs, limits = c(0,NA) ) +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::labs(title = "", x = "Theoretical", y = "Observed", color = "Test", caption = caption) +
    ggplot2::theme_minimal(18*scale, base_family = fontfamily) +
    ggplot2::theme(legend.position = "right", legend.direction = "vertical")
}