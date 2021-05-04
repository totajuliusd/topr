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
qq_plot <- function(df, scale = 1, n_variants = 0, fontfamily = "", breaks = 15) {
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
        # annotate(geom = "text", label = "Bonferroni threshold",
        #          x = min(df$theoretical), y = bonf_thres, family = fontname, vjust = -.5, size = 4*scale,
        #          hjust = 0, color = "firebrick") +
        ggplot2::geom_smooth(aes(y = exp), method = "lm", se = F, color = "#808080", size = .5, fullrange = T) +
        ggplot2::geom_point(size = scale) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(breaks), labels = log_labs, limits = c(0,NA) ) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(breaks), labels = log_labs, limits = c(0,NA) ) +
        ggplot2::scale_color_brewer(palette = "Dark2") +
        ggplot2::labs(title = "", x = "Theoretical", y = "Observed", color = "Test", caption = caption) +
        ggplot2::theme_minimal(18*scale, base_family = fontfamily) +
        ggplot2::theme(legend.position = "right", legend.direction = "vertical")
}


#' @export
#' @importFrom purrr %||%
run_qq_plot <- function(opts) {
    # get a [w,h] vector out of the resolution string, defaulting to "800x600"
    resolution <-
        opts$resolution %||% "800x600" %>%
        stringr::str_split("x", simplify = F) %>%
        unlist() %>%
        as.integer()

    scale <- readr::parse_double(opts$scale %||% "1")

    #threshold <- as.integer(opts$threshold %||% 0)
    n_variants <- as.integer(opts$nvariants %||% 0)
    #bonf_threshold <- if (n_variants > 0) 0.05 / n_variants else 0
    output_filename <- opts$outputfile

    required_columns <- c("chrom", "pos", "ref", "alt", "pval", "exp")

    columns <-
        opts[required_columns] %>%
        purrr::imap(~ .x %||% .y)

    if (!is.null(opts$color)) {
        columns <- c(columns, color = opts$color)
    }

    log_bullet("Loading input file (TSV)")
    if (!fs::file_exists(opts$inputfile))
        stop(stringr::str_glue("Input file {opts$inputfile} not found"))

    log_line("Column mapping: ", paste(names(columns), columns, sep = crayon::green("->"), collapse = ", "))
    inputfile_data <-
        opts$inputfile %>%
        data.table::fread(header = T, sep = "\t")

    # warn and exit if input data has 1 rows or less, rather than failing with stop
    # since this doesn't necessarily need to be an error in input
    if (nrow(inputfile_data) <= 1) {
        warning("The input data has 1 rows or less, nothing to do here")
        # create a placeholder plot file so sequence miner can show why it didn't generate
        generate_empty_plot(output_filename, width = resolution[1], height = resolution[2])
        return()
    }

    # See if any required columns are missing from the dataframe we just loaded
    if (!all(columns %in% colnames(inputfile_data))) {
        missing_cols <- columns[which(!columns %in% colnames(inputfile_data))]

        # If the only column missing is exp (expected p-value), no worries at this point (we will generate it later)
        if (length(missing_cols) == 1 && missing_cols == "exp") {
            columns$exp <- NULL
        } else {
            stop("The following required columns were not found in the data: ",
                 paste(missing_cols, collapse = ", "), call. = F)
        }
    }

    inputfile_data <- inputfile_data[,as.character(columns), with = F]
    inputfile_data <- purrr::set_names(inputfile_data, names(columns))

    log_line(scales::comma(nrow(inputfile_data)), " rows loaded")


    log_bullet("Generating Q-Q plot")
    log_line("File:", output_filename)
    png(output_filename, width = resolution[1], height = resolution[2])
    print(qq_plot(inputfile_data, scale = scale, n_variants = n_variants))
    dev.off()

    log_bullet("Done!")


}