#' Create a QQ plot
#'
#' @description
#'
#' \code{qqtopr()} displays QQ plots for association data.
#' Required parameter is at least one dataset (dataframe) containing the association data (with columns \code{CHROM,POS,P}
#'
#' @param dat Dataframe or a list of dataframes (required columns is \code{P})) of association results.
#' @param scale An integer, plot elements scale, default: 1
#' @param n_variants An integer, total number of variants used in the study
#' @param color A string or vector of strings setting the color/s for the input datset/s
#' @param breaks A number setting the breaks for the axes
#' @inheritParams manhattan
#'
#' @return ggplot
#' @export
#' 
qqtopr <- function(dat, scale = 1, n_variants = 0, breaks = 15, title=NULL, color=get_topr_colors(),size=1,
                   legend_name="",legend_position="right", legend_labels=NULL,
                   axis_text_size=11,axis_title_size=12, title_text_size=13,legend_title_size=12,legend_text_size=12) {
   #,
    nvars <- n_variants
    if(is.data.frame(dat)){dat <- list(dat)}
    if(is.null(legend_labels)){
        if(length(dat) > 1){
            legend_labels <- color[seq_along(dat)]
            print("Use the legend_labels argument to change the legend labels from color names to meaningful labels! ")
        }
    }

    for(i in seq_along(dat)){
        df <- dat[[i]]
        n_variants <- nvars[i]
        #df$bonf_thres <- if (n_variants > 0) -log10(0.05 / n_variants) else 0
         if (!"exp" %in% colnames(df)) {
            print("The Expected P column (exp) was not found in the data, creating it")
            if (nvars == 0) {
                n_variants <- length(df$P)
            }
         df <- df %>%
              dplyr::arrange(P) %>%
            dplyr::mutate(rank = rank(P, ties.method = "first"),
                        exp = rank/n_variants)
        }
        if (class(df$P) == "character") {
            df$P <- readr::parse_number(df$P)
        }
        if (class(df$exp) == "character") {
            df$exp <- readr::parse_number(df$exp)
        }

        zero_pvals <- dplyr::filter(df, P <= .Machine$double.xmin )
        caption <- ""
        if (nrow(zero_pvals) > 0) {
            exclusion_list <- with(zero_pvals, paste(chrom, pos, ref, alt, sep = ":", collapse = ", "))
            exclusion_list <- stringr::str_trunc(exclusion_list, width = 100)
            caption <- stringr::str_glue("{nrow(zero_pvals)} zero-value p-value{if(nrow(zero_pvals) > 1) 's' else ''} removed from plot: {exclusion_list}")
        }

        df <- df %>%
        dplyr::filter(P > .Machine$double.xmin) %>%
        dplyr::mutate(P = -log10(P),
                    exp = -log10(exp))

        df$min_theoretical <- min(df$exp)
    dat[[i]] <- df
    }

    df <- dat[[1]]
     p1 <- ggplot2::ggplot(df, aes(exp, P)) +
      ggplot2::geom_hline(aes(yintercept = val), color = "#606060", lty = 3, size = .5, data = dplyr::tibble(val = df$min_theoretical)) +
      ggplot2::geom_vline(aes(xintercept = val), color = "#606060", lty = 3, size = .5, data = dplyr::tibble(val = df$min_theoretical)) +
      ggplot2::geom_smooth(aes(y = exp), method = "lm", se = F, color = "#808080", size = .5, fullrange = T) +
      ggplot2::geom_point(size = size, aes(color=color[1])) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(breaks), limits = c(0,NA) ) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(breaks), limits = c(0,NA) ) +
      ggplot2::labs(title = title, x = "Theoretical", y = "Observed", color = "Test", caption = caption) +
       ggplot2::theme_bw()
       #ggplot2::theme_minimal(18*scale, base_family = fontfamily)
    #+
     # ggplot2::theme(legend.position = "right", legend.direction = "vertical")

    if(length(dat)>1){
        for(i in 2:length(dat)){
            df <- dat[[i]]
            p1 <- p1 + ggplot2::geom_point(data=df, aes(exp,P, color=color[i]), size = size)
              #ggplot2::geom_smooth(aes(y = exp), method = "lm", se = F, color = color[i], size = .5, fullrange = T)
        }
    }
    p1 <- p1 %>% set_plot_text_sizes(axis_text_size=axis_text_size,axis_title_size = axis_title_size, legend_text_size=legend_text_size, legend_title_size=legend_title_size,scale=scale)
   p1 <- p1 %>% change_axes()
    if(! is.null(title)){
        p1 <- p1 %>% add_title(title=title, title_text_size = title_text_size,scale=scale)
    }

    #add the legend
    if(length(dat) == 1) {
        p1 <- p1+ggplot2::scale_color_identity(breaks=color[seq_along(dat)])
    }
    else{
        p1 <- p1+ggplot2::scale_color_identity(guide = "legend", name=legend_name, breaks=color[seq_along(dat)], labels=legend_labels)+
         theme(legend.position = legend_position)
    }

  return(p1)
}

