require(ggplot2)
require(gridExtra)
require(reshape2)


Draw <- function(Signal_data, type.draw = c("signal", "pca"), output = c(
                     "default",
                     "window", "png", "pdf"
                 ), dirpath = ".", filename = "%003d", height = 480,
                 width = 640, pdf.onefile = TRUE, ...) {
    # Data initialisation and checks ----------------------------------------------
    type.draw <- match.arg(type.draw)
    output <- match.arg(output)
    fullpath <- paste(file.path(dirpath, filename), output, sep = ".")
    createFile <- TRUE
    createWindow <- FALSE

    # Drawing --------------------------------------------------------------------
    # output
    switch(output,
        default = {
            createFile <- FALSE
        },
        window = {
            createWindow <- TRUE
            createFile <- FALSE
        },
        png = {
            grDevices::png(fullpath, width, height)
        },
        pdf = {
            grDevices::pdf(fullpath,
                width = width / 72, height = height / 72,
                onefile = pdf.onefile
            )
        },
        {
            stop("Unknown output type.")
        }
    )

    # Drawing type (signal/spectrum or PCA)
    funs <- list(signal = DrawSignal, pca = DrawPCA)
    if (type.draw %in% names(funs)) {
        fun <- funs[[type.draw]]
    } else {
        stop(paste("Unknown type:", type.draw))
    }

    # Plot finalisation ----------------------------------------------
    if (is.vector(Signal_data)) {
        Signal_data <- vec2mat(Signal_data)
    }
    fun(Signal_data, createWindow = createWindow, ...)
    if (createFile) {
        grDevices::dev.off()
    }
}



#####   DrawSignal

DrawSignal <- function(Signal_data, subtype = c(
                           "stacked", "together",
                           "separate", "diffmean", "diffmedian", "diffwith"
                       ),
                       ReImModArg = c(TRUE, FALSE, FALSE, FALSE), vertical = T,
                       xlab = "rowname", RowNames = NULL, row = 1, num.stacked = 4,
                       main = NULL, createWindow) {
    # nticks

    # Data initialisation and checks ----------------------------------------------

    subtype <- match.arg(subtype)
    vec <- is.vector(Signal_data)
    if (vec) {
        Signal_data <- vec2mat(Signal_data)
    }

    n <- nrow(Signal_data)
    m <- ncol(Signal_data)

    if (n < num.stacked) {
        num.stacked <- n
    }

    scale <- colnames(Signal_data)

    num.plot <- sum(ReImModArg)

    Var <- rowname <- value <- NULL # only for R CMD check

    # Drawing array
    if (num.plot <= 0) {
        stop("Nothing selected in ReImModArg.")
    } else if (num.plot <= 2) {
        if (vertical) {
            nrow <- num.plot
            ncol <- 1
        } else {
            nrow <- 1
            ncol <- num.plot
        }
    } else {
        nrow <- 2
        ncol <- 2
    }

    # RowNames
    if (is.null(RowNames)) {
        RowNames <- rownames(Signal_data)
        if (is.null(RowNames)) {
            RowNames <- 1:n
        }
    } else {
        if (!is.vector(RowNames)) {
            stop("RowNames is not a vector")
        }
        if (length(RowNames) != n) {
            stop(paste("RowNames has length", length(RowNames), "and there are", n, "FIDs."))
        }
    }

    if (n == 1) {
        RowNames <- deparse(substitute(Signal_data))
    }

    elements <- list()
    if (ReImModArg[1]) {
        elements[["Re"]] <- Re(Signal_data)
        rownames(elements[["Re"]]) <- RowNames
    }
    if (ReImModArg[2]) {
        elements[["Im"]] <- Im(Signal_data)
        rownames(elements[["Im"]]) <- RowNames
    }
    if (ReImModArg[3]) {
        elements[["Mod"]] <- Mod(Signal_data)
        rownames(elements[["Mod"]]) <- RowNames
    }
    if (ReImModArg[4]) {
        elements[["Arg"]] <- Arg(Signal_data)
        rownames(elements[["Arg"]]) <- RowNames
    }




    # Drawing --------------------------------------------------------------------

    y <- x <- NULL # only for R CMD check


    # SEPARATE or STACKED ===============
    if (subtype == "separate" | subtype == "stacked") {
        i <- 1
        while (i <= n) {
            if (createWindow) {
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
                last <- min(i + num.stacked - 1, n)
            }
            for (name in names(elements)) {
                if (subtype == "separate") {
                    if (n == 1) {
                        df <- data.frame(x = as.numeric(scale), y = elements[[name]])
                    } else {
                        df <- data.frame(x = as.numeric(scale), y = elements[[name]][i, ])
                    }

                    plots[[name]] <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y)) +
                        ggplot2::geom_line(size = 1) +
                        ggplot2::theme(legend.position = "none") +
                        ggplot2::labs(x = xlab, y = name) +
                        ggplot2::ggtitle(RowNames[i]) +
                        ggplot2::theme_bw()

                    if ((df[1, "x"] - df[(dim(df)[1]), "x"]) > 0) {
                        plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
                    }
                } else {
                    if (n == 1) {
                        melted <- data.frame(
                            rowname = rep(name, m),
                            Var = as.numeric(scale), value = elements[[name]][i, ]
                        )
                    } else if (last == i) {
                        melted <- data.frame(
                            rowname = rep(rownames(elements[[name]])[i], m),
                            Var = as.numeric(scale), value = elements[[name]][i, ]
                        )
                    } else {
                        melted <- reshape2::melt(elements[[name]][i:last, ],
                            varnames = c("rowname", "Var")
                        )
                    }


                    plots[[name]] <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var, y = value)) +
                        ggplot2::geom_line(size = 0.3) +
                        ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
                        ggplot2::theme(legend.position = "none") +
                        ggplot2::labs(x = xlab, y = name) +
                        ggplot2::ggtitle(label = main) +
                        ggplot2::theme_bw()

                    if ((melted[1, "Var"] - melted[(dim(melted)[1]), "Var"]) > 0) {
                        plots[[name]] <- plots[[name]] + ggplot2::scale_x_reverse()
                    }
                }
            }

            if (subtype == "stacked") {
                do.call(gridExtra::grid.arrange, c(plots, list(nrow = nrow, ncol = ncol)))
            }

            i <- last + 1
        }
    } else if (subtype %in% c("together", "diffmean", "diffmedian", "diffwith")) {
        # TOGHETER or DIFFMEAN or DIFFMEDIAN or DIFFWITH ===============

        rainbow_colors <- grDevices::rainbow(n)

        if (createWindow) {
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
                element <- sweep(element, MARGIN = 2, matrixStats::colMedians(element), `-`)
            } else if (subtype == "diffwith") {
                element <- sweep(element, MARGIN = 2, element[row, ], `-`)
                if (row == 1 & n > 1) {
                    # Since we use plot on the first row and lines on the following, the y
                    # scale is calculated at the first row so if the first row is all 0, it
                    # causes problems
                    tmp <- element[1, ]
                    element[1, ] <- element[2, ]
                    element[2, ] <- tmp
                }
            }


            melted <- reshape2::melt(elements[[name]], varnames = c("rowname", "Var"))


            plots[[name]] <- ggplot2::ggplot(melted, ggplot2::aes(
                x = Var,
                y = value, group = rowname, colour = rowname
            )) +
                ggplot2::geom_line() +
                ggplot2::labs(x = xlab, y = name) +
                ggplot2::scale_colour_discrete(name = NULL) +
                ggplot2::ggtitle(main)

            if ((melted[1, "Var"] - melted[(dim(melted)[1]), "Var"]) >
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
