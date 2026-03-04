#' Plot connectivity matrix heatmap
#'
#' Visualize a group-averaged ROI x ROI connectivity matrix as a heatmap with
#' network boundary lines and network-level axis labels. Averages across
#' selected subjects and displays network structure from \code{indices}.
#'
#' @param conn_array 3D numeric array of connectivity values with dimensions
#'   (ROI x ROI x subjects), as returned by \code{\link{load_matrices}}
#' @param indices Named list of integer vectors mapping network names to ROI
#'   index positions, as returned by \code{\link{get_indices}}
#' @param subjects Integer or logical vector to subset subjects (third
#'   dimension). Defaults to all subjects. Logical vectors are converted to
#'   integer indices internally
#' @param diag Character. Controls diagonal and triangle display. One of:
#'   \itemize{
#'     \item \code{"blank"} (default): set diagonal to NA
#'     \item \code{"lower"}: show lower triangle only
#'     \item \code{"upper"}: show upper triangle only
#'   }
#' @param grid Logical. If TRUE (default), draw network boundary lines
#'   on the heatmap. Set to FALSE for a cleaner look without grid lines
#' @param title Optional character string for the plot title
#'
#' @return A \code{ggplot} object that can be further customized with
#'   standard ggplot2 syntax.
#'
#' @details
#' The function averages connectivity values across selected subjects, then
#' plots the resulting matrix using \code{ggplot2::geom_tile()} with a
#' diverging color scale centered at zero.
#'
#' Network boundary lines and axis labels are derived from the \code{indices}
#' object. ROIs are reordered to match the network ordering in \code{indices},
#' so the heatmap reflects any reordering applied by \code{\link{get_indices}}.
#'
#' Which ROIs appear in the heatmap is controlled by \code{indices}. Use
#' \code{roi_include = "schaefer"} in \code{\link{get_indices}} to exclude
#' non-Schaefer ROIs, or \code{roi_include = "all"} to include them as
#' additional blocks with their own boundaries and labels.
#'
#' Works with any Schaefer version (100-1000 parcels, 7 or 17 networks).
#'
#' @examples
#' # Basic heatmap of all subjects
#' indices <- get_indices(ex_conn_array)
#' plot_heatmap(ex_conn_array, indices)
#'
#' # Subset subjects by index
#' plot_heatmap(ex_conn_array, indices, subjects = 1:5)
#'
#' # Subset with logical vector
#' plot_heatmap(ex_conn_array, indices, subjects = c(rep(TRUE, 5), rep(FALSE, 5)))
#'
#' # Lower triangle only
#' plot_heatmap(ex_conn_array, indices, diag = "lower", title = "Lower Triangle")
#'
#' # Without boundary lines
#' plot_heatmap(ex_conn_array, indices, grid = FALSE)
#'
#' # Schaefer-only heatmap (exclude non-Schaefer ROIs)
#' indices_sch <- get_indices(ex_conn_array, roi_include = "schaefer")
#' plot_heatmap(ex_conn_array, indices_sch)
#'
#' \dontrun{
#' # Group comparison workflow
#' z_mat <- load_matrices("data/conn.mat", type = "zmat", exclude = c(46, 57))
#' indices <- get_indices(z_mat)
#'
#' # Young adults
#' plot_heatmap(z_mat, indices, subjects = 1:41, title = "Young Adults")
#'
#' # Older adults using logical vector from demographics
#' plot_heatmap(z_mat, indices, subjects = demo$group == "OA",
#'              title = "Older Adults")
#' }
#'
#' @seealso
#' \code{\link{get_indices}} for generating the \code{indices} input.
#' \code{\link{load_matrices}} for loading connectivity data.
#'
#' @export

plot_heatmap <- function(
  conn_array,
  indices,
  subjects = NULL,
  diag = c("blank", "lower", "upper"),
  grid = TRUE,
  title = NULL
) {
  # Input validation
  if (!is.array(conn_array) || length(dim(conn_array)) != 3) {
    stop(
      "conn_array must be a 3D array (ROI x ROI x subjects). ",
      "Use load_matrices() to load connectivity data.",
      call. = FALSE
    )
  }

  if (!is.list(indices) || is.null(names(indices))) {
    stop(
      "indices must be a named list of ROI index vectors. ",
      "Use get_indices() to generate indices.",
      call. = FALSE
    )
  }

  diag <- match.arg(diag)

  # Handle subject subsetting
  n_total <- dim(conn_array)[3]

  if (is.null(subjects)) {
    subjects <- seq_len(n_total)
  } else if (is.logical(subjects)) {
    if (length(subjects) != n_total) {
      stop(
        "Logical subjects vector length (",
        length(subjects),
        ") must match number of subjects (",
        n_total,
        ").",
        call. = FALSE
      )
    }
    subjects <- which(subjects)
  } else if (is.numeric(subjects)) {
    if (any(subjects < 1) || any(subjects > n_total)) {
      stop(
        "Subject indices must be between 1 and ",
        n_total,
        ".",
        call. = FALSE
      )
    }
  } else {
    stop("subjects must be an integer vector or logical vector.", call. = FALSE)
  }

  if (length(subjects) == 0) {
    stop("No subjects selected.", call. = FALSE)
  }

  # Build ROI order from indices
  roi_order <- unlist(indices, use.names = FALSE)

  # Average across selected subjects
  if (length(subjects) == 1) {
    avg_mat <- conn_array[roi_order, roi_order, subjects]
  } else {
    avg_mat <- apply(
      conn_array[roi_order, roi_order, subjects],
      c(1, 2),
      mean,
      na.rm = TRUE
    )
  }

  n_rois <- length(roi_order)

  # Apply diagonal style
  avg_mat <- switch(
    diag,
    blank = {
      diag(avg_mat) <- NA
      avg_mat
    },
    lower = {
      avg_mat[upper.tri(avg_mat, diag = TRUE)] <- NA
      avg_mat
    },
    upper = {
      avg_mat[lower.tri(avg_mat, diag = TRUE)] <- NA
      avg_mat
    }
  )

  # Reshape to long format
  plot_df <- expand.grid(row = seq_len(n_rois), col = seq_len(n_rois))
  plot_df$value <- as.vector(avg_mat)

  # Compute network boundaries and midpoints
  net_sizes <- lengths(indices)
  net_bounds <- cumsum(net_sizes) + 0.5
  midpoints <- cumsum(net_sizes) - net_sizes / 2
  net_labels <- names(indices)

  # Build heatmap
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#4682B4",
      mid = "white",
      high = "#CD5C5C",
      midpoint = 0,
      na.value = "white",
      name = "Connectivity"
    ) +
    ggplot2::scale_x_continuous(
      breaks = midpoints,
      labels = net_labels,
      expand = c(0, 0),
      position = "bottom"
    ) +
    ggplot2::scale_y_reverse(
      breaks = midpoints,
      labels = net_labels,
      expand = c(0, 0)
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  # Add network boundary lines (skip last boundary = edge of plot)
  boundary_positions <- net_bounds[-length(net_bounds)]

  if (grid && length(boundary_positions) > 0) {
    if (diag == "blank") {
      # Full-spanning lines for symmetric display
      p <- p +
        ggplot2::geom_hline(
          yintercept = boundary_positions,
          linewidth = 0.4,
          color = "black"
        ) +
        ggplot2::geom_vline(
          xintercept = boundary_positions,
          linewidth = 0.4,
          color = "black"
        )
    } else if (diag == "lower") {
      # Clip lines to lower triangle: row >= col
      plot_min <- 0.5
      plot_max <- n_rois + 0.5

      # Horizontal lines: left edge to diagonal
      h_seg <- data.frame(
        x = plot_min,
        xend = boundary_positions,
        y = boundary_positions,
        yend = boundary_positions
      )

      # Vertical lines: diagonal to bottom edge
      v_seg <- data.frame(
        x = boundary_positions,
        xend = boundary_positions,
        y = boundary_positions,
        yend = plot_max
      )

      p <- p +
        ggplot2::geom_segment(
          data = h_seg,
          ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
          linewidth = 0.4,
          color = "black",
          inherit.aes = FALSE
        ) +
        ggplot2::geom_segment(
          data = v_seg,
          ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
          linewidth = 0.4,
          color = "black",
          inherit.aes = FALSE
        )
    } else if (diag == "upper") {
      # Clip lines to upper triangle: col >= row
      plot_min <- 0.5
      plot_max <- n_rois + 0.5

      # Horizontal lines: diagonal to right edge
      h_seg <- data.frame(
        x = boundary_positions,
        xend = plot_max,
        y = boundary_positions,
        yend = boundary_positions
      )

      # Vertical lines: top edge to diagonal
      v_seg <- data.frame(
        x = boundary_positions,
        xend = boundary_positions,
        y = plot_min,
        yend = boundary_positions
      )

      p <- p +
        ggplot2::geom_segment(
          data = h_seg,
          ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
          linewidth = 0.4,
          color = "black",
          inherit.aes = FALSE
        ) +
        ggplot2::geom_segment(
          data = v_seg,
          ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
          linewidth = 0.4,
          color = "black",
          inherit.aes = FALSE
        )
    }
  }

  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}
