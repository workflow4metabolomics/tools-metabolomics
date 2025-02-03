viridis <- function(n, alpha = 1, begin = 0, end = 1, direction = 1, option = "D") {
    if (begin < 0 | begin > 1 | end < 0 | end > 1) {
        stop("begin and end must be in [0,1]")
    }
    if (abs(direction) != 1) {
        stop("direction must be 1 or -1")
    }
    if (direction == -1) {
        tmp <- begin
        begin <- end
        end <- tmp
    }
    option <- switch(EXPR = option,
        A = "A",
        magma = "A",
        B = "B",
        inferno = "B",
        C = "C",
        plasma = "C",
        D = "D",
        viridis = "D",
        E = "E",
        cividis = "E",
        {
            warning(paste0("Option '", option, "' does not exist. Defaulting to 'viridis'."))
            "D"
        }
    )
    map <- viridisLite::viridis.map[viridisLite::viridis.map$opt ==
        option, ]
    map_cols <- grDevices::rgb(map$R, map$G, map$B)
    fn_cols <- grDevices::colorRamp(map_cols,
        space = "Lab",
        interpolate = "spline"
    )
    cols <- fn_cols(seq(begin, end, length.out = n)) / 255
    grDevices::rgb(cols[, 1], cols[, 2], cols[, 3], alpha = alpha)
}
