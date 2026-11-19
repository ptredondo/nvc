generate.connect <- function(edges_input, plot.name, plot.title) {
  
  ## Basic settings
  group_names <- c("RF", "RT", "C", "O", "P", "LT", "LF")
  band_names <- c("delta", "theta", "alpha", "beta", "gamma")
  
  n_groups <- length(group_names)
  nodes_per_group <- length(band_names)
  n_nodes <- n_groups * nodes_per_group
  
  start_angle <- 82 * pi / 180
  group_gap <- 0.1
  group_label_radius <- 1.22
  
  ## Node information
  vertices <- data.frame(
    name = paste0("N", seq_len(n_nodes)),
    group = rep(group_names, each = nodes_per_group),
    band = rep(band_names, times = n_groups),
    stringsAsFactors = FALSE
  )
  
  vertex_labels <- parse(text = vertices$band)
  
  ## Frequency-band colors
  band_colors <- c(
    delta = "#B22222", theta = "#D2691E", alpha = "#DAA520",
    beta = "#228B22", gamma = "#1F4E79"
  )
  
  ## Construct graph
  g <- igraph::graph_from_data_frame(
    d = edges_input, vertices = vertices, directed = FALSE
  )
  
  ## Node and edge properties
  igraph::V(g)$color <- "white"
  igraph::V(g)$frame.color <- band_colors[igraph::V(g)$band]
  igraph::V(g)$frame.width <- 6
  
  igraph::E(g)$color <- band_colors[igraph::E(g)$band]
  igraph::E(g)$width <- igraph::E(g)$weight
  igraph::E(g)$lty <- igraph::E(g)$line_type
  
  ## Circular layout
  available_angle <- 2 * pi - n_groups * group_gap
  node_step <- available_angle / n_nodes
  
  angles <- numeric(n_nodes)
  current_angle <- start_angle
  
  for (i in seq_len(n_nodes)) {
    angles[i] <- current_angle
    current_angle <- current_angle - node_step
    
    if (i %% nodes_per_group == 0) {
      current_angle <- current_angle - group_gap
    }
  }
  
  circular_layout <- cbind(x = cos(angles), y = sin(angles))
  
  vertex_order <- match(igraph::V(g)$name, vertices$name)
  circular_layout <- circular_layout[vertex_order, , drop = FALSE]
  graph_angles <- angles[vertex_order]
  
  ## Group-label positions
  group_center_angles <- numeric(n_groups)
  
  for (j in seq_len(n_groups)) {
    group_indices <- ((j - 1) * nodes_per_group + 1):(j * nodes_per_group)
    
    group_center_angles[j] <- atan2(
      mean(sin(angles[group_indices])), mean(cos(angles[group_indices]))
    )
  }
  
  group_label_x <- group_label_radius * cos(group_center_angles)
  group_label_y <- group_label_radius * sin(group_center_angles)
  
  ## Inward edge curvature
  edge_ends <- igraph::ends(g, igraph::E(g), names = FALSE)
  
  x1 <- circular_layout[edge_ends[, 1], 1]
  y1 <- circular_layout[edge_ends[, 1], 2]
  x2 <- circular_layout[edge_ends[, 2], 1]
  y2 <- circular_layout[edge_ends[, 2], 2]
  
  angle_from <- graph_angles[edge_ends[, 1]]
  angle_to <- graph_angles[edge_ends[, 2]]
  
  angular_distance <- abs(angle_from - angle_to)
  angular_distance <- pmin(angular_distance, 2 * pi - angular_distance)
  
  max_curvature <- 1.00
  min_curvature <- 0.25
  
  curvature_magnitude <- min_curvature +
    (max_curvature - min_curvature) * (1 - angular_distance / pi)
  
  cross_product <- x1 * y2 - y1 * x2
  inward_direction <- sign(cross_product)
  inward_direction[inward_direction == 0] <- 1
  
  igraph::E(g)$curvature <- curvature_magnitude * inward_direction
  
  ## Plot connectivity network
  grDevices::png(plot.name, width = 1000, height = 1000)
  
  graphics::par(mar = c(0, 0, 2, 0), oma = c(0, 0, 0, 0), xpd = TRUE)
  
  graphics::plot(
    g, layout = circular_layout,
    vertex.color = igraph::V(g)$color,
    vertex.frame.color = igraph::V(g)$frame.color,
    vertex.size = igraph::V(g)$size,
    vertex.label = vertex_labels[vertex_order],
    vertex.label.color = "black", vertex.label.cex = 3,
    edge.color = grDevices::adjustcolor(igraph::E(g)$color, alpha.f = 0.85),
    edge.width = igraph::E(g)$width,
    edge.lty = igraph::E(g)$lty,
    edge.curved = igraph::E(g)$curvature,
    edge.arrow.mode = 0,
    asp = 1, rescale = FALSE,
    xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)
  )
  
  graphics::mtext(plot.title, side = 3, cex = 4, font = 2, line = -1.5, adj = 0)
  
  ## Add brain-region labels
  graphics::text(group_label_x, group_label_y, labels = group_names,
                 font = 2, cex = 4)
  
  grDevices::dev.off()
}


