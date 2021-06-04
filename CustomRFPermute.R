#Custom Rfpermute scripts and colors
.confMat <- function(rf) {
  if(!inherits(rf, "randomForest")) {
    stop("'rf' must be a randomForest object or inherit from one")
  }
  cbind(table(rf$y, rf$predicted))
}

plotConfMat2 <- function(rf, title = NULL, plot = TRUE) {
  conf <- .confMat(rf)
  pct.correct <- (100 * sum(diag(conf)) / sum(conf)) %>% 
    round(0) %>% 
    paste0("% correct")
  title <- if(is.null(title)) {
    pct.correct 
  } else {
    paste0(title, " (", pct.correct, ")")
  }
  freq <- rowSums(conf)
  rownames(conf) <- paste0(names(freq), " (", freq, ")")
  
  p <- conf %>% 
    prop.table(1) %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("observed") %>% 
    tidyr::gather("predicted", "prop", -.data$observed) %>% 
    dplyr::mutate(
      observed = factor(.data$observed),
      observed = stats::reorder(.data$observed, dplyr::desc(.data$observed)),
      predicted = factor(.data$predicted)
    ) %>% 
    ggplot2::ggplot(ggplot2::aes_string("predicted", "observed")) +
    ggplot2::geom_tile(ggplot2::aes_string(fill = "prop"), color = "black", size = 2) +
    ggplot2::scale_fill_viridis_c(
      option = "cividis", 
      direction = 1, 
      limits = c(0, 1)
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = "Predicted", y = "True", title = title) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Proportion")) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0),
      panel.background = ggplot2::element_blank()
    )
  
  if(plot) print(p)
  invisible(p)
}

proximityPlot2 <- function(rf, dim.x = 1, dim.y = 2, class.cols = NULL,
                          legend.type = c("legend", "label", "none"),
                          legend.loc = c("top", "bottom", "left", "right"),
                          point.size = 2, circle.size = 8, circle.border = 1, 
                          group.type = c("ellipse", "hull", "contour", "none"),
                          group.alpha = 0.3, ellipse.level = 0.95, 
                          n.contour.grid = 100, label.size = 4, 
                          label.alpha = 0.7, plot = TRUE) {
  
  if(is.null(rf$proximity)) {
    stop("'rf' has no 'proximity' element. rerun with 'proximity = TRUE'")
  }
  
  prox.mds <- stats::cmdscale(1 - rf$proximity, k = max(c(dim.x, dim.y)))
  prox.mds <- prox.mds[, c(dim.x, dim.y)]
  mds.df <- data.frame(prox.mds, class = rf$y, predicted = rf$predicted)
  colnames(mds.df)[1:2] <- c("x", "y")
  
  g <- ggplot2::ggplot(mds.df, ggplot2::aes_(~x, ~y, color = ~class)) 
  
  # Origin axes
  g <- g + 
    ggplot2::geom_hline(yintercept = 0, color = "lightgrey") +
    ggplot2::geom_vline(xintercept = 0, color = "lightgrey")
  
  legend.type <- match.arg(legend.type)
  
  # Group designators
  if(rf$type != "regression") {
    group.type <- match.arg(group.type)
    switch(
      group.type, 
      ellipse = {
        g <- g + ggplot2::stat_ellipse(
          ggplot2::aes_(fill = ~class), 
          geom = "polygon",
          alpha = group.alpha,
          level = ellipse.level,
          show.legend = legend.type == "legend"
        )
      },
      hull = {
        for(cl in unique(mds.df$class)) {
          cl.df <- mds.df[mds.df$class == cl, ]
          i <- grDevices::chull(cl.df$x, cl.df$y)
          g <- g + ggplot2::geom_polygon(
            ggplot2::aes_(fill = ~class), 
            data = cl.df[c(i, i[1]), ], 
            alpha = group.alpha,
            show.legend = legend.type == "legend"
          ) 
        }
      },
      contour = {
        g <- g + ggplot2::geom_density_2d(
          alpha = group.alpha,
          n = n.contour.grid,
          show.legend = legend.type == "legend"
        )
      }
    )
  }
  
  # Points
  if(!is.null(point.size)) {
    g <- g + ggplot2::geom_point(
      size = point.size, 
      show.legend = legend.type == "legend"
    ) 
  }
  
  # Predicted circles
  if(!is.null(circle.size)) {
    g <- g + ggplot2::geom_point(
      ggplot2::aes_(color = ~predicted), 
      shape = 21, 
      size = circle.size,
      stroke = circle.border,
      show.legend = FALSE
    )
  }
  
  # Class colors
  if(is.null(class.cols)) {
    n <- length(rf$classes)
    hues = seq(15, 375, length = n + 1)
    class.cols <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if(is.null(names(class.cols))) names(class.cols) <- rf$classes
  
  # Class labels
  if(legend.type == "label") {
    g <- g + ggplot2::geom_label(
      ggplot2::aes_(label = ~class), 
      data = mds.df %>% 
        dplyr::group_by(.data$class) %>% 
        dplyr::summarize(
          x = mean(.data$x),
          y = mean(.data$y)
        ) %>% 
        dplyr::ungroup(),
      fill = "white",
      alpha = label.alpha,
      size = label.size,
      show.legend = FALSE
    )
  }
  
  # Plot decoration
  g <- g + 
    ggplot2::scale_color_manual(values = class.cols) +
    ggplot2::scale_fill_manual(values = class.cols) +
    ggplot2::labs(
      x = paste0("MDS", dim.x), 
      y = paste0("MDS", dim.y)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = ifelse(
        legend.type == "legend",
        match.arg(legend.loc),
        "none"
      ),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(color = NA, fill = NA),
      #panel.background = ggplot2::element_blank(),
      #panel.border = ggplot2::element_blank(),
      #panel.grid.major = ggplot2::element_blank(),
      #panel.grid.minor = ggplot2::element_blank(),
      #axis.line = ggplot2::element_blank(),
      #axis.ticks = ggplot2::element_blank()
    )
  
  if(plot) print(g)
  invisible(list(prox.mds = prox.mds, g = g))
}