library(ggplot2) # nice plots
library(gridExtra) # nice plots
library(reshape2) # data manipulation

draw <- function(signal_data, type_draw = c("signal", "pca"),
                 output = c("default", "window", "png", "pdf"),
                 dirpath = ".", filename = "%003d", height = 480,
                 width = 640, pdf_onefile = TRUE, ...) {
  # Data initialisation and checks
  type_draw <- match.arg(type_draw)
  output <- match.arg(output)
  fullpath <- paste(file.path(dirpath, filename), output, sep = ".")
  create_file <- TRUE
  create_window_draw <- FALSE

  # Drawing --------------------------------------------------------------------
  # output
  switch(output,
    default = {
      create_file <- FALSE
    },
    window = {
      create_window_draw <- TRUE
      create_file <- FALSE
    },
    png = {
      grDevices::png(fullpath, width, height)
    },
    pdf = {
      grDevices::pdf(fullpath,
        width = width / 72, height = height / 72,
        onefile = pdf_onefile
      )
    },
    {
      stop("Unknown output type.")
    }
  )

  # Drawing type (signal/spectrum or PCA)
  funs <- list(signal = draw_signal, pca = DrawPCA)
  if (type_draw %in% names(funs)) {
    fun <- funs[[type_draw]]
  } else {
    stop(paste("Unknown type:", type_draw))
  }

  # Plot finalisation ----------------------------------------------
  if (is.vector(signal_data)) {
    signal_data <- vec2mat(signal_data)
  }
  fun(signal_data, create_window_draw = create_window_draw, ...)
  if (create_file) {
    grDevices::dev.off()
  }
}

#####   draw_signal
draw_signal <- function(signal_data,
                        subtype = c("stacked", "together", "separate",
                                    "diffmean", "diffmedian", "diffwith"),
                        re_im_mod_arg = c(TRUE, FALSE, FALSE, FALSE),
                        vertical = TRUE, xlab = "rowname", row_names = NULL,
                        row = 1, num_stacked = 4, main = NULL,
                        create_window_drawsignal) {
  # nticks

  # Data initialisation and checks

  subtype <- match.arg(subtype)
  vec <- is.vector(signal_data)
  if (vec) {
    signal_data <- vec2mat(signal_data)
  }

  n <- nrow(signal_data)
  m <- ncol(signal_data)

  if (n < num_stacked) {
    num_stacked <- n
  }

  scale <- colnames(signal_data)

  num_plot <- sum(re_im_mod_arg)

  variable <- rowname <- value <- NULL # only for R CMD check

  # Drawing array
  if (num_plot <= 0) {
    stop("Nothing selected in ReImModArg.")
  } else if (num_plot <= 2) {
    if (vertical) {
      nrow <- num_plot
      ncol <- 1
    } else {
      nrow <- 1
      ncol <- num_plot
    }
  } else {
    nrow <- 2
    ncol <- 2
  }

  # row_names
  if (is.null(row_names)) {
    row_names <- rownames(signal_data)
    if (is.null(row_names)) {
      row_names <- 1:n
    }
  } else {
    if (!is.vector(row_names)) {
      stop("row_names is not a vector")
    }
    if (length(row_names) != n) {
      stop(paste("row_names has length", length(row_names),
                 "and there are", n, "FIDs."))
    }
  }

  if (n == 1) {
    row_names <- deparse(substitute(signal_data))
  }

  elements <- list()
  if (re_im_mod_arg[1]) {
    elements[["Re"]] <- Re(signal_data)
    rownames(elements[["Re"]]) <- row_names
  }
  if (re_im_mod_arg[2]) {
    elements[["Im"]] <- Im(signal_data)
    rownames(elements[["Im"]]) <- row_names
  }
  if (re_im_mod_arg[3]) {
    elements[["Mod"]] <- Mod(signal_data)
    rownames(elements[["Mod"]]) <- row_names
  }
  if (re_im_mod_arg[4]) {
    elements[["Arg"]] <- Arg(signal_data)
    rownames(elements[["Arg"]]) <- row_names
  }

  # Drawing --------------------------------------------------------------------
  y <- x <- NULL # only for R CMD check

  # SEPARATE or STACKED ===============
  if (subtype == "separate" || subtype == "stacked") {
    i <- 1
    while (i <= n) {
      if (create_window_drawsignal) {
        grDevices::dev.new(noRStudioGD = TRUE)
      }
      if (subtype == "separate") {
        # The other uses gridExtra to do that
        graphics::par(mfrow = c(nrow, ncol))
      }
      plots <- list()
      if (subtype == "separate") {
        last <- i
      } else {
        last <- min(i + num_stacked - 1, n)
      }
      for (name in names(elements)) {
        if (subtype == "separate") {
          if (n == 1) {
            df <- data.frame(x = as.numeric(scale), y = elements[[name]])
          } else {
            df <- data.frame(x = as.numeric(scale), y = elements[[name]][i, ])
          }

          plots[[name]] <- ggplot2::ggplot(data = df,
                                           ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_line(size = 1) +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = xlab, y = name) +
            ggplot2::ggtitle(row_names[i]) +
            ggplot2::theme_bw()

          if ((df[1, "x"] - df[(dim(df)[1]), "x"]) > 0) {
            plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
          }
        } else {
          if (n == 1) {
            melted <- data.frame(
              rowname = rep(name, m),
              variable = as.numeric(scale), value = elements[[name]][i, ]
            )
          } else if (last == i) {
            melted <- data.frame(
              rowname = rep(rownames(elements[[name]])[i], m),
              variable = as.numeric(scale), value = elements[[name]][i, ]
            )
          } else {
            melted <- reshape2::melt(elements[[name]][i:last, ],
                                     varnames = c("rowname", "variable"))
          }


          plots[[name]] <- ggplot2::ggplot(data = melted,
                                           ggplot2::aes(x = variable, y = value)) +
            ggplot2::geom_line(size = 0.3) +
            ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = xlab, y = name) +
            ggplot2::ggtitle(label = main) +
            ggplot2::theme_bw()

          if ((melted[1, "variable"] - melted[(dim(melted)[1]), "variable"]) > 0) {
            plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
          }
        }
      }

      if (subtype == "stacked") {
        do.call(gridExtra::grid.arrange,
                c(plots, list(nrow = nrow, ncol = ncol)))
      }

      i <- last + 1
    }
  } else if (subtype %in% c("together", "diffmean", "diffmedian", "diffwith")) {
    # TOGHETER or DIFFMEAN or DIFFMEDIAN or DIFFWITH ===============

    rainbow_colors <- grDevices::rainbow(n)

    if (create_window_drawsignal) {
      grDevices::dev.new(noRStudioGD = TRUE)
    }
    graphics::par(mfrow = c(nrow, ncol))

    plots <- list()

    # Loop for Re, Im, Mod and Arg
    for (name in names(elements)) {
      # Get this part of the signal
      element <- elements[[name]]

      # Express the signal according to a reference if asked by `subtype'
      if (subtype == "diffmean") {
        element <- sweep(element, MARGIN = 2, colMeans(element), `-`)
      } else if (subtype == "diffmedian") {
        element <- sweep(element, MARGIN = 2,
                         matrixStats::colMedians(element), `-`)
      } else if (subtype == "diffwith") {
        element <- sweep(element, MARGIN = 2, element[row, ], `-`)
        if (row == 1 && n > 1) {
          # Since we use plot on the first row and lines on the following, the y

          # causes problems
          tmp <- element[1, ]
          element[1, ] <- element[2, ]
          element[2, ] <- tmp
        }
      }
      melted <- reshape2::melt(elements[[name]][i:last, ],
                               varnames = c("rowname", "variable"))

      plots[[name]] <- ggplot2::ggplot(melted, ggplot2::aes(
        x = variable,
        y = value, group = rowname, colour = rowname
      )) +
        ggplot2::geom_line() +
        ggplot2::labs(x = xlab, y = name) +
        ggplot2::scale_colour_discrete(name = NULL) +
        ggplot2::ggtitle(main)

      if ((melted[1, "variable"] - melted[(dim(melted)[1]), "variable"]) >
            0) {
        plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
      }

      do.call(gridExtra::grid.arrange, c(plots, list(
        nrow = nrow,
        ncol = ncol
      )))
    }
  }
}
