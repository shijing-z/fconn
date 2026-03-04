# ==============================================================================
# Basic Functionality Tests
# ==============================================================================

test_that("plot_heatmap exists and is callable", {
  expect_true(is.function(plot_heatmap))
})

test_that("plot_heatmap has required parameters", {
  args <- formals(plot_heatmap)

  expect_true("conn_array" %in% names(args))
  expect_true("indices" %in% names(args))
  expect_true("subjects" %in% names(args))
  expect_true("diag" %in% names(args))
  expect_true("grid" %in% names(args))
  expect_true("title" %in% names(args))
})

test_that("plot_heatmap has correct defaults", {
  args <- formals(plot_heatmap)

  expect_null(args$subjects)
  expect_false(args$grid)
  expect_null(args$title)
})

# ==============================================================================
# Input Validation Tests
# ==============================================================================

test_that("plot_heatmap errors with non-3D array", {
  indices <- get_indices(ex_conn_array)

  # 2D matrix
  expect_error(
    plot_heatmap(matrix(1:9, 3, 3), indices),
    "3D array"
  )

  # Vector
  expect_error(
    plot_heatmap(1:10, indices),
    "3D array"
  )
})

test_that("plot_heatmap errors with invalid indices", {
  # Unnamed list
  expect_error(
    plot_heatmap(ex_conn_array, list(1:5, 6:10)),
    "named list"
  )

  # Not a list
  expect_error(
    plot_heatmap(ex_conn_array, c(1, 2, 3)),
    "named list"
  )
})

test_that("plot_heatmap validates diag parameter", {
  indices <- get_indices(ex_conn_array)

  expect_error(
    plot_heatmap(ex_conn_array, indices, diag = "invalid"),
    "should be one of"
  )
})

# ==============================================================================
# Subject Subsetting Tests
# ==============================================================================

test_that("plot_heatmap errors with out-of-range integer subjects", {
  indices <- get_indices(ex_conn_array)
  n_subj <- dim(ex_conn_array)[3]

  expect_error(
    plot_heatmap(ex_conn_array, indices, subjects = 0),
    "between 1 and"
  )

  expect_error(
    plot_heatmap(ex_conn_array, indices, subjects = n_subj + 1),
    "between 1 and"
  )
})

test_that("plot_heatmap errors with wrong-length logical subjects", {
  indices <- get_indices(ex_conn_array)

  expect_error(
    plot_heatmap(ex_conn_array, indices, subjects = c(TRUE, FALSE)),
    "must match number of subjects"
  )
})

test_that("plot_heatmap errors with invalid subjects type", {
  indices <- get_indices(ex_conn_array)

  expect_error(
    plot_heatmap(ex_conn_array, indices, subjects = "first"),
    "integer vector or logical vector"
  )
})

test_that("plot_heatmap errors with empty subject selection", {
  indices <- get_indices(ex_conn_array)
  n_subj <- dim(ex_conn_array)[3]

  # Logical vector with all FALSE
  expect_error(
    plot_heatmap(ex_conn_array, indices, subjects = rep(FALSE, n_subj)),
    "No subjects selected"
  )
})

test_that("plot_heatmap accepts integer subjects", {
  indices <- get_indices(ex_conn_array)

  p <- plot_heatmap(ex_conn_array, indices, subjects = 1:5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap accepts logical subjects", {
  indices <- get_indices(ex_conn_array)
  n_subj <- dim(ex_conn_array)[3]
  logical_vec <- c(rep(TRUE, 5), rep(FALSE, n_subj - 5))

  p <- plot_heatmap(ex_conn_array, indices, subjects = logical_vec)
  expect_s3_class(p, "ggplot")
})

# ==============================================================================
# Return Value Tests
# ==============================================================================

test_that("plot_heatmap returns a ggplot object", {
  indices <- get_indices(ex_conn_array)

  p <- plot_heatmap(ex_conn_array, indices)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap returns ggplot for each diag option", {
  indices <- get_indices(ex_conn_array)

  expect_s3_class(
    plot_heatmap(ex_conn_array, indices, diag = "blank"),
    "ggplot"
  )
  expect_s3_class(
    plot_heatmap(ex_conn_array, indices, diag = "lower"),
    "ggplot"
  )
  expect_s3_class(
    plot_heatmap(ex_conn_array, indices, diag = "upper"),
    "ggplot"
  )
})

test_that("plot_heatmap returns ggplot with grid = TRUE", {
  indices <- get_indices(ex_conn_array)

  expect_s3_class(
    plot_heatmap(ex_conn_array, indices, grid = TRUE),
    "ggplot"
  )
})

test_that("plot_heatmap uses diverging fill scale", {
  indices <- get_indices(ex_conn_array)
  p <- plot_heatmap(ex_conn_array, indices)

  # Find the fill scale in the scales list
  scale_names <- vapply(
    p$scales$scales,
    function(s) s$aesthetics[1],
    character(1)
  )
  fill_idx <- which(scale_names == "fill")
  expect_length(fill_idx, 1)

  fill_scale <- p$scales$scales[[fill_idx]]
  expect_true(inherits(fill_scale, "ScaleContinuous"))
})

test_that("plot_heatmap has correct axis labels from indices", {
  indices <- get_indices(ex_conn_array, roi_include = "schaefer")
  p <- plot_heatmap(ex_conn_array, indices)

  # X-axis labels should match network names
  x_scale <- p$scales$get_scales("x")
  expect_equal(x_scale$labels, names(indices))

  # Y-axis labels should match network names
  y_scale <- p$scales$get_scales("y")
  expect_equal(y_scale$labels, names(indices))
})

test_that("plot_heatmap legend is labeled Connectivity", {
  indices <- get_indices(ex_conn_array)
  p <- plot_heatmap(ex_conn_array, indices)

  fill_scale <- p$scales$get_scales("fill")
  expect_equal(fill_scale$name, "Connectivity")
})

# ==============================================================================
# Title Tests
# ==============================================================================

test_that("plot_heatmap adds title when provided", {
  indices <- get_indices(ex_conn_array)

  p <- plot_heatmap(ex_conn_array, indices, title = "Test Title")
  expect_equal(p$labels$title, "Test Title")
})

test_that("plot_heatmap has no title when title is NULL", {
  indices <- get_indices(ex_conn_array)

  p <- plot_heatmap(ex_conn_array, indices)
  expect_null(p$labels$title)
})

# ==============================================================================
# Grid Mode Tests
# ==============================================================================

test_that("plot_heatmap with grid = FALSE has fewer layers than grid = TRUE", {
  indices <- get_indices(ex_conn_array)

  p_no_grid <- plot_heatmap(ex_conn_array, indices, grid = FALSE)
  p_grid <- plot_heatmap(ex_conn_array, indices, grid = TRUE)

  expect_true(length(p_grid$layers) > length(p_no_grid$layers))
})

test_that("plot_heatmap grid adds layers for each diag option", {
  indices <- get_indices(ex_conn_array, roi_include = "schaefer")

  blank_no <- length(
    plot_heatmap(ex_conn_array, indices, diag = "blank")$layers
  )
  blank_yes <- length(
    plot_heatmap(ex_conn_array, indices, diag = "blank", grid = TRUE)$layers
  )
  expect_true(blank_yes > blank_no)

  lower_no <- length(
    plot_heatmap(ex_conn_array, indices, diag = "lower")$layers
  )
  lower_yes <- length(
    plot_heatmap(ex_conn_array, indices, diag = "lower", grid = TRUE)$layers
  )
  expect_true(lower_yes > lower_no)

  upper_no <- length(
    plot_heatmap(ex_conn_array, indices, diag = "upper")$layers
  )
  upper_yes <- length(
    plot_heatmap(ex_conn_array, indices, diag = "upper", grid = TRUE)$layers
  )
  expect_true(upper_yes > upper_no)
})

# ==============================================================================
# ROI Inclusion Tests
# ==============================================================================

test_that("plot_heatmap works with Schaefer-only indices", {
  indices <- get_indices(ex_conn_array, roi_include = "schaefer")

  p <- plot_heatmap(ex_conn_array, indices)
  expect_s3_class(p, "ggplot")

  # Network labels should not include ahip or phip
  x_scale <- p$scales$get_scales("x")
  expect_false("ahip" %in% x_scale$labels)
  expect_false("phip" %in% x_scale$labels)
})

test_that("plot_heatmap works with all-ROI indices", {
  indices <- get_indices(ex_conn_array, roi_include = "all")

  p <- plot_heatmap(ex_conn_array, indices)
  expect_s3_class(p, "ggplot")

  # Network labels should include ahip and phip
  x_scale <- p$scales$get_scales("x")
  expect_true("ahip" %in% x_scale$labels)
  expect_true("phip" %in% x_scale$labels)
})

# ==============================================================================
# Edge Case Tests
# ==============================================================================

test_that("plot_heatmap works with single subject", {
  single_subj <- ex_conn_array[,, 1, drop = FALSE]
  indices <- get_indices(single_subj)

  p <- plot_heatmap(single_subj, indices)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap works with two networks", {
  two_net_indices <- list(
    vis = get_indices(ex_conn_array)$vis,
    default = get_indices(ex_conn_array)$default
  )

  p <- plot_heatmap(ex_conn_array, two_net_indices)
  expect_s3_class(p, "ggplot")

  x_scale <- p$scales$get_scales("x")
  expect_equal(x_scale$labels, c("vis", "default"))
})

test_that("plot_heatmap grid works with two networks", {
  two_net_indices <- list(
    vis = get_indices(ex_conn_array)$vis,
    default = get_indices(ex_conn_array)$default
  )

  # Only 1 boundary position (between the 2 networks)
  p <- plot_heatmap(ex_conn_array, two_net_indices, grid = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap works with single network", {
  one_net_indices <- list(vis = get_indices(ex_conn_array)$vis)

  # grid = TRUE but 0 boundary positions — should not error
  p <- plot_heatmap(ex_conn_array, one_net_indices, grid = TRUE)
  expect_s3_class(p, "ggplot")
})
