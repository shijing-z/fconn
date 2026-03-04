#' Plot group comparison bar chart
#'
#' Generate grouped bar plots with error bars and optional significance
#' annotations for comparing connectivity metrics between groups.
#'
#' @param data Data frame with one row per subject containing connectivity
#'   values and a grouping column
#' @param conn_vars Character vector of column names in \code{data} to
#'   compare. Each becomes a category on the x-axis
#' @param group Character string naming the grouping column in \code{data}.
#'   Must be a factor or character column with 2 or more levels
#' @param sig Logical. If TRUE (default) and \code{group} has exactly 2
#'   levels, adds significance stars from Welch's t-test above each bar
#'   cluster. Ignored when \code{group} has more than 2 levels
#' @param error_bar Character. Type of error bar to display. One of:
#'   \itemize{
#'     \item \code{"se"} (default): standard error of the mean
#'     \item \code{"sd"}: standard deviation
#'     \item \code{"ci"}: 95 percent confidence interval
#'     \item \code{"none"}: no error bars
#'   }
#' @param baseline Logical. If TRUE (default), draw a dotted horizontal
#'   reference line at y = 0
#' @param clean_labels Logical. If TRUE (default), clean column names for
#'   x-axis display by replacing underscores with hyphens and applying title
#'   case. If FALSE, use raw column names
#' @param title Optional character string for the plot title
#'
#' @return A \code{ggplot} object that can be further customized with
#'   standard ggplot2 syntax.
#'
#' @details
#' The function computes per-group summary statistics (mean and error bars)
#' from the input data frame, then generates a grouped bar plot with one
#' cluster per connectivity variable.
#'
#' Significance stars use Welch's t-test (two-sided, unequal variance)
#' with the following thresholds: *** p < 0.001, ** p < 0.01, * p < 0.05.
#' Stars are only available for exactly 2 groups; set \code{sig = FALSE}
#' or use more than 2 groups to suppress.
#'
#' The function works with any column names, not limited to \code{calc_*}
#' output. Any numeric columns in the data frame can be used as
#' \code{conn_vars}.
#'
#' X-axis labels can always be overridden with
#' \code{+ scale_x_discrete(labels = ...)} regardless of the
#' \code{clean_labels} setting.
#'
#' @examples
#' # Within-network comparison by group (constructed difference for demo)
#' indices <- get_indices(ex_conn_array, roi_include = "schaefer")
#' within_df <- calc_within(ex_conn_array, indices)
#' within_df$within_default <- within_df$within_default +
#'   rep(c(0.3, 0), times = c(5, 5))
#' within_df$group <- rep(c("YA", "OA"), times = c(5, 5))
#'
#' plot_compare(within_df,
#'   conn_vars = c("within_default", "within_vis"),
#'   group = "group"
#' )
#'
#' # ROI-to-network comparison without significance stars
#' ahip_df <- calc_conn(ex_conn_array, indices,
#'   from = "ahip", to = c("default", "cont", "vis")
#' )
#' ahip_df$group <- rep(c("YA", "OA"), times = c(5, 5))
#'
#' plot_compare(ahip_df,
#'   conn_vars = c("ahip_default", "ahip_cont", "ahip_vis"),
#'   group = "group", sig = FALSE
#' )
#'
#' # Custom labels override clean_labels
#' plot_compare(within_df,
#'   conn_vars = c("within_default", "within_vis"),
#'   group = "group"
#' ) + ggplot2::scale_x_discrete(labels = c("DMN", "VIS"))
#'
#' \dontrun{
#' # Full workflow with real data
#' z_mat <- load_matrices("data/conn.mat", type = "zmat", exclude = c(46, 57))
#' indices <- get_indices(z_mat, roi_include = "schaefer")
#' within_df <- calc_within(z_mat, indices)
#' within_df$group <- demo$group
#'
#' plot_compare(within_df,
#'   conn_vars = c("within_default", "within_cont", "within_vis"),
#'   group = "group",
#'   title = "Within-Network Connectivity by Age Group"
#' )
#' }
#'
#' @seealso
#' \code{\link{plot_heatmap}} for connectivity matrix heatmaps.
#' \code{\link{plot_scatter}} for connectivity-behavior scatter plots.
#' \code{\link{calc_within}}, \code{\link{calc_between}},
#'   \code{\link{calc_conn}} for generating connectivity data.
#'
#' @export

plot_compare <- function(
  data,
  conn_vars,
  group,
  sig = TRUE,
  error_bar = c("se", "sd", "ci", "none"),
  baseline = TRUE,
  clean_labels = TRUE,
  title = NULL
) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }

  if (!is.character(conn_vars) || length(conn_vars) == 0) {
    stop(
      "conn_vars must be a character vector with at least one column name.",
      call. = FALSE
    )
  }

  missing_vars <- conn_vars[!conn_vars %in% names(data)]
  if (length(missing_vars) > 0) {
    stop(
      "Column(s) not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.character(group) || length(group) != 1) {
    stop(
      "group must be a single character string naming a column in data.",
      call. = FALSE
    )
  }

  if (!group %in% names(data)) {
    stop(
      "Group column '",
      group,
      "' not found in data.",
      call. = FALSE
    )
  }

  # Check that conn_vars are numeric
  non_numeric <- conn_vars[!sapply(data[conn_vars], is.numeric)]
  if (length(non_numeric) > 0) {
    stop(
      "conn_vars must be numeric columns. Non-numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }

  error_bar <- match.arg(error_bar)

  # Ensure group is a factor
  data[[group]] <- as.factor(data[[group]])
  group_levels <- levels(data[[group]])
  n_groups <- length(group_levels)

  if (n_groups < 2) {
    stop(
      "group column must have at least 2 levels. Found: ",
      n_groups,
      call. = FALSE
    )
  }

  # Pivot to long format
  long_list <- lapply(conn_vars, function(var) {
    data.frame(
      group_var = data[[group]],
      variable = var,
      connectivity = data[[var]],
      stringsAsFactors = FALSE
    )
  })

  plot_long <- do.call(rbind, long_list)
  plot_long$variable <- factor(plot_long$variable, levels = conn_vars)

  # Compute summary statistics
  summary_df <- do.call(
    rbind,
    lapply(
      split(
        plot_long,
        list(plot_long$group_var, plot_long$variable),
        drop = TRUE
      ),
      function(chunk) {
        vals <- chunk$connectivity
        n <- sum(!is.na(vals))
        m <- mean(vals, na.rm = TRUE)
        s <- stats::sd(vals, na.rm = TRUE)

        err <- switch(
          error_bar,
          se = s / sqrt(n),
          sd = s,
          ci = stats::qt(0.975, df = n - 1) * s / sqrt(n),
          none = 0
        )

        data.frame(
          group_var = chunk$group_var[1],
          variable = chunk$variable[1],
          mean_conn = m,
          error = err,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(summary_df) <- NULL
  summary_df$variable <- factor(summary_df$variable, levels = conn_vars)

  # Clean labels
  if (clean_labels) {
    display_labels <- gsub("_", "-", conn_vars)
    display_labels <- sapply(
      display_labels,
      function(label) {
        parts <- strsplit(label, "-")[[1]]
        parts <- tools::toTitleCase(parts)
        paste(parts, collapse = "-")
      },
      USE.NAMES = FALSE
    )
  } else {
    display_labels <- conn_vars
  }

  # Build plot
  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = variable,
      y = mean_conn,
      fill = group_var
    )
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.7
    ) +
    ggplot2::scale_x_discrete(labels = display_labels) +
    ggplot2::labs(
      x = NULL,
      y = "Connectivity",
      fill = "Group"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.text.x = ggplot2::element_text(
        hjust = 0.5,
        vjust = 1,
        face = "bold"
      ),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  # Add error bars
  if (error_bar != "none") {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = mean_conn - error,
          ymax = mean_conn + error
        ),
        position = ggplot2::position_dodge(width = 0.7),
        width = 0.2,
        color = "black"
      )
  }

  # Add baseline
  if (baseline) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted",
        color = "black",
        linewidth = 0.5
      )
  }

  # Add significance stars (2-group only)
  if (sig && n_groups == 2) {
    sig_df <- do.call(
      rbind,
      lapply(conn_vars, function(var) {
        group1 <- data[[var]][data[[group]] == group_levels[1]]
        group2 <- data[[var]][data[[group]] == group_levels[2]]

        p_val <- stats::t.test(group1, group2)$p.value

        star <- if (p_val < 0.001) {
          "***"
        } else if (p_val < 0.01) {
          "**"
        } else if (p_val < 0.05) {
          "*"
        } else {
          ""
        }

        data.frame(
          variable = var,
          star = star,
          stringsAsFactors = FALSE
        )
      })
    )

    sig_df$variable <- factor(sig_df$variable, levels = conn_vars)

    # Position stars above the tallest bar + error bar in each cluster
    max_heights <- do.call(
      rbind,
      lapply(
        split(summary_df, summary_df$variable),
        function(chunk) {
          data.frame(
            variable = chunk$variable[1],
            y_pos = max(chunk$mean_conn + chunk$error, na.rm = TRUE)
          )
        }
      )
    )

    sig_df <- merge(sig_df, max_heights, by = "variable")
    sig_df$variable <- factor(sig_df$variable, levels = conn_vars)

    # Add padding above tallest bar
    y_range <- max(summary_df$mean_conn + summary_df$error, na.rm = TRUE) -
      min(summary_df$mean_conn - summary_df$error, na.rm = TRUE)
    sig_df$y_pos <- sig_df$y_pos + y_range * 0.05

    # Only annotate significant results
    sig_df <- sig_df[sig_df$star != "", ]

    if (nrow(sig_df) > 0) {
      p <- p +
        ggplot2::geom_text(
          data = sig_df,
          ggplot2::aes(
            x = variable,
            y = y_pos,
            label = star
          ),
          inherit.aes = FALSE,
          size = 6,
          fontface = "bold"
        )
    }
  }

  # Add title
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}
