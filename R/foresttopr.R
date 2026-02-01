
##foresttopr( dat = list(dat_main,CD_FINNGEN,UC_UKBB), legend_labels = c("UKBB","FinnGen","UC"),xbreaks=c(0.7,0.9,1.5,2,3),color=c("skyblue4","skyblue3","red"),
##legend_position = "right",legend_nrow = 3,points_dist = 0.6,size=2, key_col = "Gene_Symbol",band_border_linewidth=0.000001)
# =========================
# Forest plot pipeline with:
# - effect_type = c("OR","beta") (case-insensitive), default OR
# - mixed inputs: OR-only converted to beta (log OR) when effect_type="beta";
#                beta-only converted to OR (exp beta) when effect_type="OR"
# - CI preference: use LCL/UCL if present; else SE if present; else derive SE from p
# - correct null line + axis for each effect_type:
#     OR   -> log10 scale, null = 1, labeled breaks
#     beta -> linear scale, null = 0, pretty breaks
# =========================

# ---- helper: find first matching column name (case-insensitive) ----
find_col_ci <- function(df, candidates) {
  nms <- names(df)
  idx <- match(toupper(candidates), toupper(nms), nomatch = 0L)
  idx <- idx[idx != 0L]
  if (length(idx) == 0) return(NULL)
  nms[idx[1]]
}

# ---- CI helper: prefer LCL/UCL > SE > P-derived SE ----
compute_ci_prefer <- function(est,
                              p = NA_real_,
                              se = NA_real_,
                              lcl = NA_real_,
                              ucl = NA_real_,
                              effect_type = c("OR", "beta")) {
  
  effect_type <- tolower(effect_type[1])
  stopifnot(effect_type %in% c("or", "beta"))
  
  # 0) Use explicit CI if present
  if (is.finite(lcl) && is.finite(ucl)) {
    return(c(lcl, ucl))
  }
  
  # 1) Use SE if present
  if (is.finite(se) && se > 0 && is.finite(est)) {
    if (effect_type == "or") {
      if (est <= 0) return(c(NA_real_, NA_real_))
      # IMPORTANT: se is assumed to be SE(log(OR)) if provided from logistic regression
      return(c(exp(log(est) - 1.96 * se),
               exp(log(est) + 1.96 * se)))
    } else {
      return(c(est - 1.96 * se,
               est + 1.96 * se))
    }
  }
  
  # 2) Fallback: derive SE from p-value
  if (!is.finite(est) || is.na(p) || p <= 0 || p > 1) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(p / 2, lower.tail = FALSE)
  if (!is.finite(z) || z <= 0) return(c(NA_real_, NA_real_))
  
  if (effect_type == "or") {
    if (est <= 0) return(c(NA_real_, NA_real_))
    se2 <- abs(log(est)) / z
    c(exp(log(est) - 1.96 * se2),
      exp(log(est) + 1.96 * se2))
  } else {
    se2 <- abs(est) / z
    c(est - 1.96 * se2,
      est + 1.96 * se2)
  }
}

# ---- standardize: effect + p + optional se/lcl/ucl into *_std fields ----
standardize_effects <- function(df, effect_type = c("OR", "beta")) {
  effect_type <- tolower(effect_type[1])
  stopifnot(effect_type %in% c("or", "beta"))
  
  # effect columns (case-insensitive)
  or_col   <- find_col_ci(df, c("OR"))
  beta_col <- find_col_ci(df, c("BETA"))
  eff_col  <- find_col_ci(df, c("EFFECT"))
  
  # p column (case-insensitive)
  p_col <- find_col_ci(df, c("P", "PVAL", "PVALUE", "P_VALUE", "P.VALUE"))
  
  # SE columns (case-insensitive, common variants)
  se_col <- find_col_ci(df, c("SE", "STDERR", "STD_ERR", "SE_BETA", "BETA_SE", "SEBETA"))
  
  # CI columns (case-insensitive, common variants)
  lcl_col <- find_col_ci(df, c("LCL", "LCI", "LOWER", "LOWER_CI", "CI_LOWER", "L95", "LOW95"))
  ucl_col <- find_col_ci(df, c("UCL", "UCI", "UPPER", "UPPER_CI", "CI_UPPER", "U95", "UP95"))
  
  # ---- standardize estimate to requested scale ----
  if (effect_type == "or") {
    if (!is.null(or_col)) {
      df$est_std <- suppressWarnings(as.numeric(df[[or_col]]))
    } else if (!is.null(beta_col)) {
      df$est_std <- exp(suppressWarnings(as.numeric(df[[beta_col]])))
    } else if (!is.null(eff_col)) {
      df$est_std <- exp(suppressWarnings(as.numeric(df[[eff_col]])))
    } else {
      stop("No OR/BETA/EFFECT column found to derive OR.")
    }
  } else { # beta
    if (!is.null(beta_col)) {
      df$est_std <- suppressWarnings(as.numeric(df[[beta_col]]))
    } else if (!is.null(eff_col)) {
      df$est_std <- suppressWarnings(as.numeric(df[[eff_col]]))
    } else if (!is.null(or_col)) {
      df$est_std <- log(suppressWarnings(as.numeric(df[[or_col]])))
    } else {
      stop("No BETA/EFFECT/OR column found to derive beta.")
    }
  }
  
  # ---- standardize uncertainty inputs (optional) ----
  df$p_std   <- if (!is.null(p_col))  suppressWarnings(as.numeric(df[[p_col]])) else NA_real_
  df$se_std  <- if (!is.null(se_col)) suppressWarnings(as.numeric(df[[se_col]])) else NA_real_
  df$lcl_std <- if (!is.null(lcl_col)) suppressWarnings(as.numeric(df[[lcl_col]])) else NA_real_
  df$ucl_std <- if (!is.null(ucl_col)) suppressWarnings(as.numeric(df[[ucl_col]])) else NA_real_
  
  df
}

# ---- detect gene column name ----
detect_gene_col <- function(df, gene_col = NULL) {
  if (!is.null(gene_col) && gene_col %in% names(df)) {
    return(gene_col)
  }
  candidates <- c("Gene_Symbol", "gene", "GENE", "Gene")
  cand <- candidates[candidates %in% names(df)]
  if (length(cand) == 0) {
    stop(
      "Could not find a key/gene column (tried: Gene_Symbol, gene, GENE, Gene). ",
      "Provide `key_col` explicitly."
    )
  }
  cand[1]
}


# =========================
# Matching functions
# =========================
#' Match association results across datasets by a key column
#'
#' @description
#' `get_matching_genes()` aligns rows from multiple association result tables
#' using a shared key column (e.g. gene or feature identifier). Effect estimates
#' and confidence intervals are standardized across datasets, while row labels
#' are taken exclusively from the reference dataset (the first element of `dfs`).
#'
#' This function is typically used internally by \code{\link{foresttopr}},
#' but may be useful on its own when preparing matched effect tables for
#' visualization or downstream analysis.
#'
#' @param dfs A list of data frames containing association results. Each data
#'   frame must contain a key column and effect size information.
#'
#' @param labels A character vector of dataset labels of the same length as
#'   `dfs`. These labels identify the source of each matched effect estimate.
#'
#' @param gene_col Character scalar specifying the column name used to match
#'   rows across datasets (e.g. gene identifier). If `NULL`, a suitable column
#'   is inferred from common gene identifier names.
#'
#' @param label_col Optional character scalar specifying the column name in the
#'   reference dataset (the first element of `dfs`) to use as a human-readable
#'   row label. If `NULL`, the matching key (`gene_col`) is used for labeling.
#'
#' @param effect_type Character scalar specifying the effect scale to use.
#'   Either `"OR"` (odds ratio) or `"beta"` (regression coefficient).
#'   Matching is case-insensitive. Effect estimates are converted between
#'   scales as needed.
#'
#' @details
#' Rows are matched across datasets using the key column specified by
#' `gene_col`. The set of keys present in the reference dataset defines the
#' universe of rows retained. For each dataset, confidence intervals are
#' derived preferentially from explicit bounds, standard errors, or p-values,
#' depending on availability.
#'
#' The returned table contains one row per matched key per dataset.
#'
#' @return A data frame containing matched effect estimates with the following
#'   columns:
#'   \describe{
#'     \item{key}{Matching key used to align rows across datasets.}
#'     \item{label}{Row label used for display purposes.}
#'     \item{set}{Dataset identifier corresponding to `labels`.}
#'     \item{or}{Effect estimate on the requested scale.}
#'     \item{p}{P-value associated with the effect estimate.}
#'     \item{lcl}{Lower confidence interval bound.}
#'     \item{ucl}{Upper confidence interval bound.}
#'   }
#'
#' @examples
#' matched <- get_matching_genes(
#'   dfs = list(CD_tmp, CD_FINNGEN),
#'   labels = c("CD_UKBB", "CD_FINNGEN"),
#'   gene_col = "ID",
#'   label_col = "Gene_Symbol",
#'   effect_type = "beta"
#' )
#'
#' @seealso \code{\link{foresttopr}}
#'
#' @keywords internal

# === match on gene (key) only ===
get_matching_genes <- function(
    dfs,
    labels,
    gene_col = NULL,     # this is the MATCH key column
    label_col = NULL,    # this is the DISPLAY label column (from ref only)
    effect_type = c("OR", "beta")
) {
  effect_type <- tolower(effect_type[1])
  stopifnot(effect_type %in% c("or", "beta"))
  stopifnot(length(dfs) == length(labels))
  
  # standardize effect columns first
  dfs_std <- lapply(dfs, standardize_effects, effect_type = effect_type)
  
  # ----- reference dataset (first one) -----
  key_ref <- detect_gene_col(dfs_std[[1]], gene_col)
  
  # label column comes from reference only
  if (is.null(label_col)) {
    label_ref <- key_ref
  } else {
    if (!label_col %in% names(dfs_std[[1]])) {
      stop("label_col='", label_col, "' not found in reference dataset (first element of dat).")
    }
    label_ref <- label_col
  }
  
  ref <- dfs_std[[1]] %>%
    dplyr::mutate(
      key   = .data[[key_ref]],
      label = .data[[label_ref]]
    )
  
  ref_out <- ref %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ci  = list(compute_ci_prefer(
        est = est_std, p = p_std, se = se_std, lcl = lcl_std, ucl = ucl_std,
        effect_type = effect_type
      )),
      lcl = ci[1],
      ucl = ci[2],
      set = labels[1]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      key,
      label,
      set,
      or  = est_std,
      p   = p_std,
      lcl,
      ucl
    )
  
  # keys + labels we want to keep for alignment
  ref_keys <- ref %>%
    dplyr::distinct(key, label)
  
  # ----- other datasets -----
  others_out <- dplyr::bind_rows(
    purrr::map2(dfs_std[-1], labels[-1], function(df, lab) {
      key_other <- detect_gene_col(df, gene_col)
      
      df2 <- df %>%
        dplyr::mutate(
          key  = .data[[key_other]],
          or   = est_std,
          p    = p_std,
          se   = se_std,
          lcl0 = lcl_std,
          ucl0 = ucl_std
        ) %>%
        dplyr::select(key, or, p, se, lcl0, ucl0)
      
      ref_keys %>%
        dplyr::left_join(df2, by = "key") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          ci  = list(if (!is.na(or)) compute_ci_prefer(
            est = or, p = p, se = se, lcl = lcl0, ucl = ucl0,
            effect_type = effect_type
          ) else c(NA_real_, NA_real_)),
          lcl = ci[1],
          ucl = ci[2]
        ) %>%
        dplyr::ungroup() %>%
        dplyr::transmute(
          key,
          label,       # keep ref label
          set = lab,
          or,
          p,
          lcl,
          ucl
        ) %>%
        dplyr::filter(!is.na(or))  # keep only matches
    })
  )
  
  dplyr::bind_rows(ref_out, others_out) %>%
    dplyr::mutate(p = suppressWarnings(as.numeric(p)))
}


# === match on SNP key (CHROM:POS) with allele flip handling ===

#' Match association results across datasets by variant key
#'
#' @description
#' `get_matching_snps()` aligns rows from multiple variant-level association
#' result tables using a genomic variant key constructed from chromosome and
#' position (e.g. `"CHR:POS"`). Effect estimates and confidence intervals are
#' standardized across datasets, while row labels are taken from the reference
#' dataset (the first element of `dfs`).
#'
#' The function supports allele-aware matching and automatically corrects
#' effect directions when reference and alternate alleles are flipped between
#' datasets.
#'
#' This function is typically used internally by \code{\link{foresttopr}} when
#' variant-level columns are detected, but may also be called directly to
#' prepare matched variant-level effect tables.
#'
#' @param dfs A list of data frames containing variant-level association
#'   results. Each data frame must contain chromosome and position columns
#'   (`CHROM`, `POS`) and effect size information.
#'
#' @param labels A character vector of dataset labels of the same length as
#'   `dfs`. These labels identify the source of each matched effect estimate.
#'
#' @param gene_col Character scalar specifying a column name in the reference
#'   dataset used for labeling rows (e.g. gene or variant identifier).
#'   Defaults to `"ID"`. This column is not used for matching.
#'
#' @param label_col Optional character scalar specifying an alternative column
#'   name in the reference dataset to use as a human-readable row label.
#'   If `NULL`, `gene_col` is used for labeling.
#'
#' @param effect_type Character scalar specifying the effect scale to use.
#'   Either `"OR"` (odds ratio) or `"beta"` (regression coefficient).
#'   Matching is case-insensitive. Effect estimates are converted between
#'   scales as needed.
#'
#' @details
#' Variants are matched across datasets using a key constructed from
#' chromosome and position. Reference and alternate alleles are compared
#' between datasets, and effect estimates are automatically flipped when
#' allele orientation differs.
#'
#' Confidence intervals are derived preferentially from explicit bounds,
#' standard errors, or p-values, depending on availability.
#'
#' The returned table contains one row per matched variant per dataset.
#'
#' @return A data frame containing matched variant-level effect estimates with
#'   the following columns:
#'   \describe{
#'     \item{key}{Variant key constructed from chromosome and position.}
#'     \item{label}{Row label used for display purposes.}
#'     \item{set}{Dataset identifier corresponding to `labels`.}
#'     \item{or}{Effect estimate on the requested scale.}
#'     \item{p}{P-value associated with the effect estimate.}
#'     \item{lcl}{Lower confidence interval bound.}
#'     \item{ucl}{Upper confidence interval bound.}
#'   }
#'
#' @examples
#' matched <- get_matching_snps(
#'   dfs = list(CD_tmp, CD_FINNGEN),
#'   labels = c("CD_UKBB", "CD_FINNGEN"),
#'   gene_col = "ID",
#'   label_col = "Gene_Symbol",
#'   effect_type = "beta"
#' )
#'
#' @seealso \code{\link{foresttopr}}, \code{\link{get_matching_genes}}
#'
#' @keywords internal

get_matching_snps <- function(
    dfs,
    labels,
    gene_col = "ID",     # used only for label if you want; you can keep as-is
    label_col = NULL,    # NEW: label from ref dataset
    effect_type = c("OR", "beta")
) {
  effect_type <- tolower(effect_type[1])
  stopifnot(effect_type %in% c("or", "beta"))
  stopifnot(length(dfs) == length(labels))
  
  dfs_std <- lapply(dfs, standardize_effects, effect_type = effect_type)
  
  # reference
  ref <- dfs_std[[1]] %>%
    dplyr::mutate(
      key = paste0(CHROM, ":", POS),
      label = if (!is.null(label_col) && label_col %in% names(.)) .data[[label_col]] else .data[[gene_col]],
      or = est_std,
      p  = p_std,
      se = se_std,
      lcl0 = lcl_std,
      ucl0 = ucl_std
    )
  
  ref_out <- ref %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ci  = list(compute_ci_prefer(or, p, se, lcl0, ucl0, effect_type = effect_type)),
      lcl = ci[1],
      ucl = ci[2],
      set = labels[1]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(key, label, set, or, p, lcl, ucl)
  
  others_out <- dplyr::bind_rows(
    purrr::map2(dfs_std[-1], labels[-1], function(df, lab) {
      
      df2 <- df %>%
        dplyr::mutate(
          key  = paste0(CHROM, ":", POS),
          or   = est_std,
          p    = p_std,
          se   = se_std,
          lcl0 = lcl_std,
          ucl0 = ucl_std
        ) %>%
        dplyr::select(key, REF, ALT, or, p, se, lcl0, ucl0)
      
      joined <- ref %>%
        dplyr::select(key, label, REF_ref = REF, ALT_ref = ALT) %>%
        dplyr::left_join(df2, by = "key") %>%
        dplyr::mutate(
          flip = dplyr::case_when(
            !is.na(REF) & REF == REF_ref & ALT == ALT_ref ~ FALSE,
            !is.na(REF) & REF == ALT_ref & ALT == REF_ref ~ TRUE,
            TRUE ~ NA
          ),
          or_adj = dplyr::case_when(
            is.na(flip) ~ NA_real_,
            !flip       ~ or,
            flip & effect_type == "or"   ~ 1 / or,
            flip & effect_type == "beta" ~ -or,
            TRUE ~ NA_real_
          )
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          ci  = list(if (!is.na(or_adj)) compute_ci_prefer(
            est = or_adj, p = p, se = se, lcl = lcl0, ucl = ucl0,
            effect_type = effect_type
          ) else c(NA_real_, NA_real_)),
          lcl = ci[1],
          ucl = ci[2]
        ) %>%
        dplyr::ungroup() %>%
        dplyr::transmute(key, label, set = lab, or = or_adj, p, lcl, ucl) %>%
        dplyr::filter(!is.na(or))
      
      joined
    })
  )
  
  dplyr::bind_rows(ref_out, others_out) %>%
    dplyr::mutate(p = suppressWarnings(as.numeric(p)))
}


# =========================
# Forest plot function
# =========================

#' Create a forest plot from one or more association result tables
#'
#' @description
#' `foresttopr()` creates a forest plot visualizing effect estimates and
#' confidence intervals across one or more datasets. The function supports
#' odds ratios (OR) and regression coefficients (beta), allows matching rows
#' across datasets by a key column, and optionally displays human-readable
#' labels from a separate annotation column.
#'
#' Effect estimates are automatically standardized across input datasets,
#' and confidence intervals are derived preferentially from explicit bounds,
#' standard errors, or p-values when necessary.
#'
#' @param dat A data frame or a list of data frames containing association
#'   results. Each data frame must contain an effect estimate column (e.g.
#'   OR or BETA) and a p-value column. If a single data frame is provided,
#'   it is internally wrapped into a list.
#'
#' @param legend_labels A character vector of labels corresponding to each
#'   dataset in `dat`. These labels are used in the plot legend. Defaults to
#'   `"Set1"`, `"Set2"`, etc.
#'
#' @param colors A character vector of colors to use for each dataset.
#'   If `NULL`, a default color palette is used.
#'
#' @param key_col Character scalar giving the column name used to match rows
#'   across datasets (e.g. gene identifier or variant ID). Defaults to `"gene"`.
#'
#' @param label_col Optional character scalar giving the column name in the
#'   reference dataset (the first element of `dat`) to use for labeling rows
#'   on the y-axis. If `NULL`, `key_col` is used for labeling.
#'
#' @param effect_type Character scalar specifying the effect scale to plot.
#'   Either `"OR"` (odds ratio; default) or `"beta"` (regression coefficient).
#'   Matching is case-insensitive. When required, effect estimates are
#'   automatically converted between scales.
#'
#' @param xlim Numeric length-2 vector giving x-axis limits. If `NULL`,
#'   limits are computed automatically from the data.
#'
#' @param xbreaks Numeric vector or function specifying x-axis breaks.
#'   If `NULL`, reasonable defaults are chosen based on `effect_type`.
#'
#' @param xlabel Character scalar giving the x-axis label. If `NULL`,
#'   a default label is chosen based on `effect_type`.
#'
#' @param size Numeric scalar or vector controlling point sizes for each dataset.
#'
#' @param shape Integer scalar or vector specifying point shapes for each dataset.
#'
#' @param alpha Numeric scalar or vector specifying point transparency.
#'
#' @param points_dist Numeric scalar controlling horizontal separation of
#'   points from different datasets within the same row.
#'
#' @param band_color Background color for alternating row bands.
#'
#' @param band_border_color Color for row band borders.
#'
#' @param band_border_linewidth Numeric scalar giving the line width for
#'   row band borders.
#'
#' @param sign_thresh Optional numeric scalar specifying a p-value threshold
#'   for highlighting statistically significant points via shape encoding.
#'
#' @param ylabel_order Optional character vector specifying the order of
#'   rows on the y-axis. If `NULL`, rows are ordered as they appear in the
#'   reference dataset.
#'
#' @param scale Numeric scalar used to globally scale text and point sizes.
#'
#' @param title Optional character scalar giving the plot title.
#'
#' @param title_text_size Numeric scalar controlling title text size.
#'
#' @param axis_text_size Numeric scalar controlling axis text size.
#'
#' @param axis_title_size Numeric scalar controlling axis title text size.
#'
#' @param show_shape_legend Logical; whether to display the shape legend.
#'
#' @param show_color_legend Logical; whether to display the color legend.
#'
#' @param legend_position Character string specifying legend position.
#'   One of `"right"`, `"top"`, or `"bottom"`.
#'
#' @param legend_nrow Optional integer specifying the number of rows in the legend.
#'
#' @param legend_name Optional character scalar giving the legend title.
#'
#' @param legend_title_size Numeric scalar controlling legend title text size.
#'
#' @param legend_text_size Numeric scalar controlling legend text size.
#'
#' @param match_on_gene Logical; if `FALSE` and vdariant-level columns
#'   (e.g. REF/ALT) are detected, matching is performed at the variant level.
#'   Otherwise, matching is performed using `key_col`.
#'
#' @return A `ggplot2` object representing the forest plot.
#'
#' @examples
#' foresttopr(
#'   dat = list(CD_UKBB %>% arrange(P) %>% head(n=10) %>% annotate_with_nearest_gene(), CD_FINNGEN),
#'   key_col = "ID",
#'   label_col = "Gene_Symbol",
#'   legend_labels = c("CD_UKBB", "CD_FINNGEN"),
#'   effect_type = "beta"
#' )
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#'
#' @export


foresttopr <- function(
    dat = NULL,
    legend_labels = NULL,
    colors = NULL,
    key_col = "ID",
    label_col = NULL,            # NEW
    effect_type = c("OR", "beta"),
    xlim = NULL,
    xbreaks = NULL,
    xlabel = NULL,
    size = 2.5,
    shape = 16,
    alpha = 1,
    points_dist = 0.6,
    band_color = "grey96",
    band_border_color = "grey96",
    band_border_linewidth = 0.01,
    sign_thresh = NULL,
    ylabel_order = NULL,
    scale = 1,
    title = NULL,
    title_text_size = 14,
    axis_text_size = 12,
    axis_title_size = 13,
    show_shape_legend = TRUE,
    show_color_legend = TRUE,
    legend_position = "right",
    legend_nrow = NULL,
    legend_name = NULL,
    legend_title_size = axis_text_size * 0.95,
    legend_text_size  = axis_text_size * 0.85,
    match_on_gene = FALSE
) {
  effect_type <- tolower(effect_type[1])
  stopifnot(effect_type %in% c("or", "beta"))
  
  # ---- normalize inputs ----
  dat <- if (is.data.frame(dat)) list(dat) else if (is.list(dat)) dat else list(dat)
  if (is.null(legend_labels)) legend_labels <- paste0("Set", seq_along(dat))
  stopifnot(length(dat) == length(legend_labels))
  n_sets <- length(dat)
  
  if (is.null(colors)) {
    colors <- topr::get_topr_colors()[seq_len(n_sets)]
  } else {
    colors <- colors[seq_len(n_sets)]
  }
  stopifnot(length(colors) == n_sets)
  
  legend_position <- match.arg(legend_position, c("right", "top", "bottom"))
  if (is.null(legend_nrow)) legend_nrow <- n_sets
  legend_nrow <- as.integer(legend_nrow)
  
  if (is.null(xlabel)) {
    xlabel <- if (effect_type == "or") "Odds ratio (95% CI)" else "Effect size (β, 95% CI)"
  }
  
  recycle_to_n <- function(x, n, nm) {
    if (length(x) == 1) rep(x, n) else if (length(x) == n) x else
      stop(sprintf("`%s` must have length 1 or %d.", nm, n))
  }
  size  <- recycle_to_n(size * scale, n_sets, "size")
  shape <- recycle_to_n(shape,        n_sets, "shape")
  alpha <- recycle_to_n(alpha,        n_sets, "alpha")
  
  cohort_cols   <- setNames(colors, legend_labels)
  cohort_shapes <- setNames(shape,  legend_labels)
  cohort_sizes  <- setNames(size,   legend_labels)
  cohort_alpha  <- setNames(alpha,  legend_labels)
  
  # ---- matching strategy ----
  has_ref_alt <- function(df) all(c("ref", "alt") %in% tolower(names(df)))
  all_have_ref_alt <- all(vapply(dat, has_ref_alt, logical(1)))
  
  df <- if (!match_on_gene && all_have_ref_alt) {
    get_matching_snps(dfs = dat, labels = legend_labels, gene_col = key_col,
                      label_col = label_col, effect_type = effect_type)
  } else {
    get_matching_genes(dfs = dat, labels = legend_labels, gene_col = key_col,
                       label_col = label_col, effect_type = effect_type)
  }
  
  df <- dplyr::mutate(df, set = factor(set, levels = legend_labels))
  df$p <- suppressWarnings(as.numeric(df$p))
  
  # ---- choose y-axis label variable ----
  y_var <- if ("label" %in% names(df)) "label" else if ("gene" %in% names(df)) "gene" else stop("No label/gene column found.")
  df[[y_var]] <- as.character(df[[y_var]])
  
  # ---- significance shape override ----
  if (!is.null(sign_thresh)) {
    stopifnot(is.numeric(sign_thresh), length(sign_thresh) == 1, is.finite(sign_thresh),
              sign_thresh > 0, sign_thresh < 1)
    df$.sig <- ifelse(!is.na(df$p) & df$p <= sign_thresh, "sig", "ns")
    sig_shapes <- c(sig = 16, ns = 1)
  }
  
  # ---- x limits & breaks ----
  if (effect_type == "or") {
    if (is.null(xlim)) {
      xmin <- min(df$lcl, na.rm = TRUE) * 0.95
      xmax <- max(df$ucl, na.rm = TRUE) * 1.05
      xlim <- c(min(0.7, xmin), xmax)
    }
    if (xlim[1] <= 0) xlim[1] <- 1e-5
    if (is.null(xbreaks)) xbreaks <- function(lims) sort(unique(c(1, scales::log_breaks()(lims))))
    break_vals <- if (is.function(xbreaks)) xbreaks(xlim) else xbreaks
    break_vals <- sort(unique(break_vals[is.finite(break_vals) & break_vals > 0]))
  } else {
    if (is.null(xlim)) {
      xmin <- min(df$lcl, na.rm = TRUE)
      xmax <- max(df$ucl, na.rm = TRUE)
      pad  <- 0.05 * (xmax - xmin)
      if (!is.finite(pad) || pad == 0) pad <- 0.1
      xlim <- c(xmin - pad, xmax + pad)
    }
    if (is.null(xbreaks)) {
      break_vals <- scales::pretty_breaks(n = 5)(xlim)
    } else if (is.function(xbreaks)) {
      break_vals <- xbreaks(xlim)
    } else {
      break_vals <- xbreaks
    }
    break_vals <- break_vals[is.finite(break_vals)]
  }
  
  # ---- y ordering + offsets ----
  y_order <- if (is.null(ylabel_order)) dplyr::pull(dplyr::distinct(df, .data[[y_var]]), .data[[y_var]]) else ylabel_order
  df <- dplyr::mutate(df, y_lbl = factor(.data[[y_var]], levels = rev(y_order)), y = as.numeric(y_lbl))
  y_labs <- dplyr::arrange(dplyr::distinct(df, y_lbl, y), y)
  
  bands <- y_labs |>
    dplyr::mutate(idx = dplyr::row_number()) |>
    dplyr::filter(idx %% 2 == 0) |>
    dplyr::transmute(ymin = y - 0.5, ymax = y + 0.5, xmin = xlim[1], xmax = xlim[2])
  
  set_offsets <- tibble::tibble(
    set = factor(legend_labels, levels = legend_labels),
    off = (seq_along(legend_labels) - (n_sets + 1) / 2) * (points_dist / max(n_sets, 1))
  )
  df <- dplyr::left_join(df, set_offsets, by = "set")
  df$y_plot <- df$y + df$off
  
  # ---- legend policy ----
  legend_point_scale <- 0.75
  guide_color <- if (isTRUE(show_color_legend)) ggplot2::guide_legend(nrow = legend_nrow, byrow = TRUE) else "none"
  guide_shape <- if (!isTRUE(show_color_legend) && isTRUE(show_shape_legend)) {
    ggplot2::guide_legend(nrow = legend_nrow, byrow = TRUE,
                          override.aes = list(size = mean(size) * legend_point_scale))
  } else "none"
  
  # ---- plot ----
  p <- ggplot2::ggplot(df, ggplot2::aes(or, y_plot, color = set)) +
    ggplot2::geom_rect(data = bands, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = band_color, color = NA) +
    ggplot2::geom_segment(data = bands, inherit.aes = FALSE,
                          ggplot2::aes(x = xmin, xend = xmax, y = ymin, yend = ymin),
                          color = band_border_color, linewidth = band_border_linewidth) +
    ggplot2::geom_segment(data = bands, inherit.aes = FALSE,
                          ggplot2::aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
                          color = band_border_color, linewidth = band_border_linewidth) +
    ggplot2::scale_y_continuous(breaks = y_labs$y, labels = as.character(y_labs$y_lbl),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lcl,
                                        xmax = ucl, 
                                        group = interaction(y_lbl, set)),
                           orientation = "y", width = 0, linewidth = 0.6)
  
  # dashed reference lines: filter to finite + within xlim to avoid warnings
  vbreaks <- break_vals
  vbreaks <- vbreaks[is.finite(vbreaks)]
  vbreaks <- vbreaks[vbreaks >= xlim[1] & vbreaks <= xlim[2]]
  
  if (effect_type == "or") {
    vbreaks <- vbreaks[vbreaks > 0]
    p <- p +
      ggplot2::geom_vline(xintercept = 1, linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = vbreaks[vbreaks != 1],
                          linetype = "22", linewidth = 0.4, color = "grey50") +
      ggplot2::scale_x_log10(limits = xlim, breaks = break_vals, labels = scales::label_number())
  } else {
    p <- p +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = vbreaks[vbreaks != 0],
                          linetype = "22", linewidth = 0.4, color = "grey50") +
      ggplot2::scale_x_continuous(limits = xlim, breaks = break_vals)
  }
  
  # points + shape scale
  if (!is.null(sign_thresh)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(shape = .sig, size = set, alpha = set,
                                       group = interaction(y_lbl, set))) +
      ggplot2::scale_shape_manual(values = sig_shapes, name = NULL,
                                  breaks = c("sig", "ns"),
                                  labels = c(paste0("P \u2264 ", sign_thresh),
                                             paste0("P > ", sign_thresh)))
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(shape = set, size = set, alpha = set,
                                       group = interaction(y_lbl, set))) +
      ggplot2::scale_shape_manual(values = cohort_shapes, breaks = legend_labels,
                                  drop = FALSE, name = legend_name)
  }
  
  p <- p +
    ggplot2::scale_color_manual(values = cohort_cols, breaks = legend_labels,
                                drop = FALSE, name = legend_name, guide = guide_color) +
    ggplot2::scale_size_manual(values = setNames(size, legend_labels), breaks = legend_labels,
                               drop = FALSE, name = legend_name, guide = "none") +
    ggplot2::scale_alpha_manual(values = setNames(alpha, legend_labels), breaks = legend_labels,
                                drop = FALSE, name = legend_name, guide = "none") +
    ggplot2::guides(color = guide_color, shape = guide_shape, size = "none", alpha = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = xlabel, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 13 * scale) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text  = ggplot2::element_text(size = axis_text_size * scale),
      axis.title = ggplot2::element_text(size = axis_title_size * scale),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = legend_title_size * scale),
      legend.text  = ggplot2::element_text(size = legend_text_size * scale),
      legend.box = if (legend_position %in% c("top", "bottom")) "horizontal" else "vertical",
      plot.margin = ggplot2::margin(10, 30, 10, 10)
    )
  
  if (!is.null(title)) {
    p <- p + ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_text_size * scale, hjust = 0, margin = ggplot2::margin(b = 8)
      )
    )
  }
  
  p
}
