## ==========================
# Internal functions
## ==========================

# beginTreatment
begin_treatment <- function(name, signal_data = NULL, signal_info = NULL,
                            force_real = FALSE) {
  cat("Begin", name, "\n")


  # Formatting the signal_data and signal_info -----------------------
  vec <- is.vector(signal_data)
  if (vec) {
    signal_data <- vec2mat(signal_data)
  }
  if (is.vector(signal_info)) {
    signal_info <- vec2mat(signal_info)
  }
  if (!is.null(signal_data)) {
    if (!is.matrix(signal_data)) {
      stop("signal_data is not a matrix.")
    }
    if (!is.complex(signal_data) && !is.numeric(signal_data)) {
      stop("signal_data contains non-numerical values.")
    }
  }
  if (!is.null(signal_info) && !is.matrix(signal_info)) {
    stop("signal_info is not a matrix.")
  }


  original_data <- signal_data

  # Extract the real part of the spectrum ---------------------------
  if (force_real) {
    if (is.complex(signal_data)) {
      signal_data <- Re(signal_data)
    } else {
      # The signal is numeric Im(signal_data) is zero anyway so let's avoid
      # using complex(real=...,imaginary=0) which would give a complex signal
      # in end_treatment()
      force_real <- FALSE
    }
  }

  # Return the formatted data and metadata entries --------------------

  return(list(start = proc.time(), vec = vec, force_real = force_real,
    original_data = original_data, signal_data = signal_data,
    signal_info = signal_info
  ))
}

# end_treatment
end_treatment <- function(name, begin_info, signal_data) {
  # begin_info: object outputted from beginTreatment
  # Formatting the entries and printing process time -----------------------
  end_time <- proc.time() # record it as soon as possible
  start_time <- begin_info[["start"]]
  delta_time <- end_time - start_time
  delta <- delta_time[]
  cat("End", name, "\n")
  cat("It lasted", round(delta["user.self"], 3),
    "s user time,", round(delta["sys.self"], 3),
    "s system time and", round(delta["elapsed"], 3), "s elapsed time.\n"
  )

  if (begin_info[["force_real"]]) {
    # The imaginary part is left untouched
    i <- complex(real = 0, imaginary = 1)
    signal_data <- signal_data + i * Im(begin_info[["original_data"]])
  }

  if (begin_info[["vec"]]) {
    signal_data <- signal_data[1, ]
  }

  # Return the formatted data and metadata entries --------------------
  signal_data
}

# check_arg
check_arg <- function(arg, checks, can_be_null = FALSE) {
  check_list <- list(bool = c(is.logical, "a boolean"),
    int = c(function(x) {
      x %% 1 == 0
    }, "an integer"),
    num = c(is.numeric, "a numeric"),
    str = c(is.character, "a string"),
    pos = c(function(x) {
      x > 0
    }, "positive"),
    pos0 = c(function(x) {
      x >= 0
    }, "positive or zero"),
    l1 = c(function(x) {
      length(x) == 1
    }, "of length 1")
  )
  if (is.null(arg)) {
    if (!can_be_null) {
      stop(deparse(substitute(arg)), " is null.")
    }
  } else {
    if (is.matrix(arg)) {
      stop(deparse(substitute(arg)), " is not scalar.")
    }
    for (c in checks) {
      if (!check_list[[c]][[1]](arg)) {
        stop(deparse(substitute(arg)), " is not ", check_list[[c]][[2]], ".")
      }
    }
  }
}

# get_arg
get_arg <- function(arg, info, argname, can_be_absent = FALSE) {
  if (is.null(arg)) {
    start <- paste("impossible to get argument", argname,
                   "it was not given directly and")
    if (!is.matrix(info)) {
      stop(paste(start, "the info matrix was not given"))
    }
    if (!(argname %in% colnames(info))) {
      if (can_be_absent) {
        return(NULL)
      } else {
        stop(paste(start, "is not in the info matrix"))
      }
    }
    if (nrow(info) < 1) {
      stop(paste(start, "the info matrix has no row"))
    }
    arg <- info[1, argname]
    if (is.na(arg)) {
      stop(paste(start, "it is NA in the info matrix"))
    }
  }
  return(arg)
}

# binary_search
binary_search <- function(a, target, lower = TRUE) {
  # search the index i in a such that a[i] == target
  # if it doesn't exists and lower,
  # it searches the closer a[i] such that a[i] < target
  # if !lower, it seraches the closer a[i] such that a[i] > target
  # a should be monotone but can be increasing or decreasing

  # if a is increasing INVARIANT: a[amin] < target < a[amax]
  n <- length(a)
  if ((a[n] - target) * (a[n] - a[1]) <= 0) {
    return(n)
  }
  if ((a[1] - target) * (a[n] - a[1]) >= 0) {
    return(1)
  }
  amin <- 1
  amax <- n
  while (amin + 1 < amax) {
    amid <- floor((amin + amax) / 2)
    if ((a[amid] - target) * (a[amax] - a[amid]) < 0) {
      amin <- amid
    } else if ((a[amid] - target) * (a[amax] - a[amid]) > 0) {
      amax <- amid
    } else {
      # target
      return(amid)
    }
  }
  if (xor(lower, a[amin] > a[amax])) {
    # If increasing and we want the lower, we take amin
    # If decreasing and we want the bigger, we take amin too
    return(amin)
  } else {
    return(amax)
  }
}

# Interpol
interpol <- function(t, y) {
  # y: sample
  # t : warping function

  m <- length(y)
  # because if t > m-1, y[ti+1] will be NA when we compute g
  valid <- 1 <= t & t <= m - 1 # FIXME it was '<' in Bubble v2
  s <- (1:m)[valid]
  ti <- floor(t[s])
  tr <- t[s] - ti
  g <- y[ti + 1] - y[ti]
  f <- y[ti] + tr * g
  list(f = f, s = s, g = g)
}

# vec2mat
vec2mat <- function(vec) {
  matrix(vec, nrow = 1, dimnames = list(c(1), names(vec)))
}

# binary_search
binary_search <- function(a, target, lower = TRUE) {
  # search the index i in a such that a[i] == target
  # if it doesn't exists and lower,
  # it searches the closer a[i] such that a[i] < target
  # if !lower, it seraches the closer a[i] such that a[i] > target
  # a should be monotone but can be increasing or decreasing

  # if a is increasing INVARIANT: a[amin] < target < a[amax]
  n <- length(a)
  if ((a[n] - target) * (a[n] - a[1]) <= 0) {
    return(n)
  }
  if ((a[1] - target) * (a[n] - a[1]) >= 0) {
    return(1)
  }
  amin <- 1
  amax <- n
  while (amin + 1 < amax) {
    amid <- floor((amin + amax) / 2)
    if ((a[amid] - target) * (a[amax] - a[amid]) < 0) {
      amin <- amid
    } else if ((a[amid] - target) * (a[amax] - a[amid]) > 0) {
      amax <- amid
    } else {
      # target
      return(amid)
    }
  }
  if (xor(lower, a[amin] > a[amax])) {
    # If increasing and we want the lower, we take amin
    # If decreasing and we want the bigger, we take amin too
    return(amin)
  } else {
    return(amax)
  }
}

# indexInterval
index_interval <- function(a, from, to, inclusive = TRUE) {
  # If inclusive and from <= to, we need to take the lower
  # If not inclusive and from > to, we need to take the lower too
  lower_from <- xor(inclusive, from > to)
  from_index <- binary_search(a, from, lower_from)
  to_index <- binary_search(a, to, !lower_from)
  return(from_index:to_index)
}

## ==========================
# GroupDelayCorrection
## ==========================
group_delay_correction <- function(fid_data, fid_info = NULL,
                                   group_delay = NULL) {
  # Data initialisation and checks
  begin_info <- begin_treatment("GroupDelayCorrection", fid_data, fid_info)
  fid_data <- begin_info[["signal_data"]]
  dimension_names <- dimnames(fid_data)
  fid_info <- begin_info[["signal_info"]]
  check_arg(group_delay, c("num", "pos0"), can_be_null = TRUE)
  # if fid_info and group_delay are NULL, get_arg will generate an error

  group_delay <- get_arg(group_delay, fid_info, "GRPDLY", can_be_absent = TRUE)

  if (is.null(group_delay)) {
    # See DetermineBrukerDigitalFilter.m in matNMR MATLAB library
    group_delay_matrix <- matrix(
      c(44.75, 46, 46.311, 33.5, 36.5, 36.53, 66.625, 48, 47.87, 59.0833,
        50.1667, 50.229, 68.5625, 53.25, 53.289, 60.375, 69.5, 69.551,
        69.5313, 72.25, 71.6, 61.0208, 70.1667, 70.184, 70.0156, 72.75,
        72.138, 61.3438, 70.5, 70.528, 70.2578, 73, 72.348, 61.5052, 70.6667,
        70.7, 70.3789, 72.5, 72.524, 61.5859, 71.3333, NA, 70.4395, 72.25, NA,
        61.6263, 71.6667, NA, 70.4697, 72.125, NA, 61.6465, 71.8333, NA,
        70.4849, 72.0625, NA, 61.6566, 71.9167, NA, 70.4924, 72.0313, NA
      ),
      nrow = 21,
      ncol = 3, byrow = TRUE, dimnames = list(
        c(2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384,
          512, 768, 1024, 1536, 2048
        ),
        c(10, 11, 12)
      )
    )
    decim <- fid_info[1, "DECIM"]
    dspfvs <- fid_info[1, "DSPFVS"]
    if (!(toString(decim) %in% rownames(group_delay_matrix))) {
      stop(paste("Invalid DECIM", decim, "it should be one of",
                 rownames(group_delay_matrix)))
    }
    if (!(toString(dspfvs) %in% colnames(group_delay_matrix))) {
      stop(paste("Invalid DSPFVS", dspfvs, "it should be one of",
                 colnames(group_delay_matrix)))
    }
    group_delay <- group_delay_matrix[toString(decim), toString(dspfvs)]
    if (is.na(group_delay)) {
      stop(paste("Invalid DECIM", decim, "for DSPFVS", dspfvs))
    }
  }
  m <- ncol(fid_data)
  n <- nrow(fid_data)

  # GroupDelayCorrection
  # We do the shifting in the Fourier domain: the shift can be non-integer.
  # That way we automatically have the circular behaviour of the shift and the
  # interpolation if it is non-integer.
  spectrum <- t(stats::mvfft(t(fid_data)))

  p <- ceiling(m / 2)
  new_index <- c((p + 1):m, 1:p)
  spectrum <- spectrum[, new_index]
  spectrum <- matrix(data = spectrum, ncol = m, nrow = n)

  omega <- (0:(m - 1)) / m
  i <- complex(real = 0, imaginary = 1)

  if (n > 1) {
    spectrum <- sweep(spectrum, MARGIN = 2,
                      exp(i * group_delay * 2 * pi * omega), `*`)
    spectrum <- spectrum[, new_index]
  } else {
    spectrum <- spectrum * exp(i * group_delay * 2 * pi * omega)
    spectrum <- spectrum[new_index]
    spectrum <- matrix(data = spectrum, ncol = m, nrow = n)
  }

  fid_data <- t(stats::mvfft(t(spectrum), inverse = TRUE)) / m
  colnames(fid_data) <- dimension_names[[2]]
  rownames(fid_data) <- dimension_names[[1]]

  # Data finalisation
  return(end_treatment("GroupDelayCorrection", begin_info, fid_data))
}

## ==========================
# SolventSuppression
## ==========================
solvent_suppression <- function(fid_data, lambda_ss = 1e+06, ptw_ss = TRUE,
                                plot_solvent = FALSE, return_solvent = FALSE) {
  # Data initialisation and checks
  begin_info <- begin_treatment("SolventSuppression", fid_data)
  fid_data <- begin_info[["signal_data"]]
  check_arg(ptw_ss, c("bool"))
  check_arg(lambda_ss, c("num", "pos0"))

  # difsm function definition for the smoother
  if (ptw_ss) {
    # Use of the function in ptw that smoothes signals with a finite
    # difference penalty of order 2
    difsm <- ptw::difsm
  } else {
    # Or manual implementation based on sparse matrices for large data series
    # (cf. Eilers, 2003. 'A perfect smoother')
    difsm <- function(y, d = 2, lambda) {
      m <- length(y)
      # Sparse identity matrix m x m
      e <- Matrix::Diagonal(m)
      dd <- Matrix::diff(e, differences = d)
      a <- e + lambda_ss * Matrix::t(dd) %*% dd
      # base::chol does not take into account that a is sparse
      # and is extremely slow
      c <- Matrix::chol(a)
      x <- Matrix::solve(c, Matrix::solve(Matrix::t(c), y))
      return(as.numeric(x))
    }
  }

  # Solvent Suppression ----------------------------------------------
  n <- dim(fid_data)[1]
  if (return_solvent) {
    solvent_re <- fid_data
    solvent_im <- fid_data
  }
  for (i in 1:n) {
    fid_re <- Re(fid_data[i, ])
    fid_im <- Im(fid_data[i, ])
    solvent_re <- difsm(y = fid_re, lambda = lambda_ss)
    solvent_im <- difsm(y = fid_im, lambda = lambda_ss)

    if (plot_solvent) {
      m <- length(fid_re)
      graphics::plot(1:m, fid_re, type = "l", col = "red")
      graphics::lines(1:m, solvent_re, type = "l", col = "blue")
      graphics::plot(1:m, fid_im, type = "l", col = "red")
      graphics::lines(1:m, solvent_im, type = "l", col = "blue")
    }
    fid_re <- fid_re - solvent_re
    fid_im <- fid_im - solvent_im
    fid_data[i, ] <- complex(real = fid_re, imaginary = fid_im)
    if (return_solvent) {
      solvent_re[i, ] <- solvent_re
      solvent_im[i, ] <- solvent_im
    }
  }

  # Data finalisation ----------------------------------------------
  fid_data <- end_treatment("SolventSuppression", begin_info, fid_data)
  if (return_solvent) {
    return(list(fid_data = fid_data, solvent_re = solvent_re,
                solvent_im = solvent_im))
  } else {
    return(fid_data)
  }
}


## ==========================
# Apodization
# =============================
apodization <- function(fid_data, fid_info = NULL, dt = NULL,
                        type_apod = c("exp", "cos2", "blockexp",
                                      "blockcos2", "gauss", "hanning",
                                      "hamming"),
                        phase = 0, rect_ratio = 1 / 2,
                        gauss_lb = 1, exp_lb = 1, plot_window = FALSE,
                        return_factor = FALSE) {
  # Data initialisation and checks
  begin_info <- begin_treatment("Apodization", fid_data, fid_info)
  fid_data <- begin_info[["signal_data"]]
  fid_info <- begin_info[["signal_info"]]
  # Data check
  type_apod <- match.arg(type_apod)
  check_arg(dt, c("num", "pos"), can_be_null = TRUE)
  check_arg(phase, c("num"))

  # Apodization ----------------------------------------------
  dt <- get_arg(dt, fid_info, "DT") # Dwell Time
  m <- ncol(fid_data)
  t <- (1:m) * dt # Time
  rect_size <- ceiling(rect_ratio * m)
  gauss_lb <- (gauss_lb / (sqrt(8 * log(2))))
  # Define the types of apodization:
  switch(type_apod,
    exp = {
      # exponential
      factor <- exp(-exp_lb * t)
    },
    cos2 = {
      # cos2
      c <- cos((1:m) * pi / (2 * m) - phase * pi / 2)
      factor <- c * c
    },
    blockexp = {
      # block and exponential
      factor <- c(rep.int(1, rect_size), rep.int(0, m - rect_size))
      # | rect_size | 1 ___________ | \ 0 \____
      factor[(rect_size + 1):m] <- exp(-exp_lb * t[1:(m - rect_size)])
    },
    blockcos2 = {
      # block and cos^2
      factor <- c(rep.int(1, rect_size), rep.int(0, m - rect_size))
      c <- cos((1:(m - rect_size)) * pi / (2 * (m - rect_size)))
      factor[(rect_size + 1):m] <- c * c
    },
    gauss = {
      # gaussian
      factor <- exp(-(gauss_lb * t)^2 / 2)
      factor <- factor / max(factor)
    },
    hanning = {
      # Hanning
      factor <- 0.5 + 0.5 * cos((1:m) * pi / m - phase * pi)
    },
    hamming = {
      # Hamming
      factor <- 0.54 + 0.46 * cos((1:m) * pi / m - phase * pi)
    }
  )
  if (plot_window) {
    graphics::plot(1:m, factor, "l")
  }

  # Apply the apodization factor on the spectra
  fid_data <- sweep(fid_data, MARGIN = 2, factor, `*`)

  # Data finalisation ----------------------------------------------
  fid_data <- end_treatment("Apodization", begin_info, fid_data)
  if (return_factor) {
    return(list(fid_data = fid_data, factor = factor))
  } else {
    return(fid_data)
  }
}


## ====================================================
# FourierTransform
## ====================================================

fftshift_1d_2d <- function(x) {
  vec <- FALSE
  if (is.vector(x)) {
    x <- vec2mat(x)
    vec <- TRUE
  }
  m <- dim(x)[2]
  p <- ceiling(m / 2)
  new_index <- c((p + 1):m, 1:p)
  y <- x[, new_index, drop = vec]
}

# FourierTransform
fourier_transform <- function(fid_data, fid_info = NULL, sw_h = NULL,
                              sw = NULL, o1 = NULL, reverse_axis = TRUE) {
  # Data initialisation and checks
  begin_info <- begin_treatment("FourierTransform", fid_data, fid_info)
  fid_data <- begin_info[["signal_data"]]
  fid_info <- begin_info[["signal_info"]]

  m <- ncol(fid_data)
  n <- nrow(fid_data)

  if (is.null(sw_h)) {
    sw_h <- get_arg(sw_h, fid_info, "SW_h")
  }

  if (is.null(sw)) {
    # Sweep Width in ppm (semi frequency scale in ppm)
    sw <- get_arg(sw, fid_info, "SW")
  }

  if (is.null(o1)) {
    o1 <- get_arg(o1, fid_info, "O1")
  }

  check_arg(reverse_axis, c("bool"))

  # Fourier Transformation
  # mvfft does the unnormalized fourier transform (see ?mvfft),
  # so we need divide by m.  It does not matter a lot in our case
  # since the spectrum will be normalized.

  # FT
  raw_spect_data <- fftshift_1d_2d(t(stats::mvfft(t(fid_data))))
  # recover the frequencies values
  f <- ((0:(m - 1)) - floor(m / 2)) * fid_info[1, "SW_h"] / (m - 1)

  if (reverse_axis == TRUE) {
    revind <- rev(1:m)
    # reverse the spectrum
    raw_spect_data <- raw_spect_data[, revind]
  }

  raw_spect_data <- matrix(raw_spect_data, nrow = n, ncol = m)
  colnames(raw_spect_data) <- f
  rownames(raw_spect_data) <- rownames(fid_data)

  # PPM conversion
  # The Sweep Width has to be the same since the column names are the same
  ppm_interval <- sw / (m - 1)

  o1_index <- round((m + 1) / 2 + o1 * (m - 1) / sw_h)

  end <- o1_index - m
  start <- o1_index - 1
  ppm_scale <- (start:end) * ppm_interval
  raw_spect_data <- matrix(raw_spect_data,
    nrow = n,
    ncol = -(end - start) + 1,
    dimnames = list(rownames(raw_spect_data), ppm_scale)
  )

  # Data finalisation ----------------------------------------------
  return(end_treatment("FourierTransform", begin_info, raw_spect_data))
}

## ====================================================
#   InternalReferencing
## ====================================================
internal_referencing <- function(spectrum_data, fid_info,
                                 method = c("max", "thres"),
                                 range = c("nearvalue", "all", "window"),
                                 ppm_value = 0,
                                 direction = "left",
                                 shift_handling = c("zerofilling", "cut",
                                                    "NAfilling", "circular"),
                                 c = 2, pc = 0.02,
                                 fromto_rc = NULL, ppm_ir = TRUE,
                                 rowindex_graph = NULL) {

  # Data initialisation and checks
  begin_info <- begin_treatment("InternalReferencing", spectrum_data, fid_info)
  spectrum_data <- begin_info[["signal_data"]]
  fid_info <- begin_info[["signal_info"]]

  ######## Check input arguments
  range <- match.arg(range)
  shift_handling <- match.arg(shift_handling)
  method <- match.arg(method)
  plots <- NULL

  check_arg(ppm_ir, c("bool"))
  check_arg(unlist(fromto_rc), c("num"), can_be_null = TRUE)
  check_arg(pc, c("num"))
  check_arg(ppm_value, c("num"))
  check_arg(rowindex_graph, "num", can_be_null = TRUE)

  # fromto_rc : if range == "window",
  # fromto_rc defines the spectral window where to search for the peak
  if (!is.null(fromto_rc)) {
    diff <- diff(unlist(fromto_rc))[seq_along(diff(unlist(fromto_rc)))
                                    %% 2 != 0]
    for (i in seq_len(diff)) {
      if (diff[i] >= 0) {
        fromto <- c(fromto_rc[[i]][2], fromto_rc[[i]][1])
        fromto_rc[[i]] <- fromto
      }
    }
  }

  # find_tmsp_peak function ----------------------------------------------
  # If method == "tresh", find_tmsp_peak will find the position of the first
  # peak (from left or right) which is higher than a predefined threshold
  # and is computed as: c*(cumulated_mean/cumulated_sd)
  find_tmsp_peak <- function(ft, c = 2, direction = "left") {
    ft <- Re(ft) # extraction de la partie réelle
    n <- length(ft)
    if (direction == "left") {
      newindex <- rev(1:n)
      ft <- rev(ft)
    }
    thres <- 99999
    i <- 1000 # Start at point 1000 to find the peak
    vect <- ft[1:i]

    while (vect[i] <= (c * thres)) {
      cumsd <- stats::sd(vect)
      cummean <- mean(vect)
      thres <- cummean + 3 * cumsd
      i <- i + 1
      vect <- ft[1:i]
    }
    if (direction == "left") {
      v <- newindex[i]
    } else {
      v <- i
    }

    if (is.na(v)) {
      warning("No peak found, need to lower the threshold.")
      return(NA)
    } else {
      # recherche dans les 1% de points suivants du max trouve
      # pour etre au sommet du pic
      d <- which.max(ft[v:(v + n * 0.01)])
      new_peak <- v + d - 1 # nouveau pic du TMSP si d > 0
      if (names(which.max(ft[v:(v + n * 0.01)])) !=
            names(which.max(ft[v:(v + n * 0.03)]))) {
        # recherche dans les 3% de points suivants du max trouve
        # pour eviter un fauxpositif
        warning("the TMSP peak might be located further away,
                increase the threshold to check.")
      }
      return(new_peak)
    }
  }

  # Define the search zone  ----------------------------------------
  n <- nrow(spectrum_data)
  m <- ncol(spectrum_data)

  # The Sweep Width (sw) has to be the same since the column names are the same
  sw <- fid_info[1, "SW"] # Sweep Width in ppm
  ppm_interval <- sw / (m - 1) # size of a ppm interval

  # range: How the search zone is defined ("all", "nearvalue" or "window")
  if (range == "all") {
    data <- spectrum_data
  } else {
    # range = "nearvalue" or "window"
    # Need to define colindex (column indexes) to apply indexInterval on it
    if (range == "nearvalue") {
      # automatic fromto values in ppm
      fromto_rc <- list(c(-(sw * pc) / 2 + ppm_value,
                          (sw * pc) / 2 + ppm_value))
      colindex <- as.numeric(colnames(spectrum_data))
    } else {
      # fromto_rc is already user-defined
      if (ppm_ir == TRUE) {
        colindex <- as.numeric(colnames(spectrum_data))
      } else {
        colindex <- 1:m
      }
    }

    # index intervals taking into account the different elements
    # in the list fromto_rc
    interval <- vector("list", length(fromto_rc))
    for (i in seq_along(fromto_rc)) {
      interval[[i]] <- index_interval(colindex,
                                      from = fromto_rc[[i]][1],
                                      to = fromto_rc[[i]][2],
                                      inclusive = TRUE)
    }

    # define Data as the cropped spectrum including the index intervals
    # outside the research zone, the intensities are set to the minimal
    # intensity of the research zone
    if (n > 1) {
      data <- apply(Re(spectrum_data[, unlist(interval)]), 1,
                    function(x) rep(min(x), m))
      data <- t(data)
      data[, unlist(interval)] <- Re(spectrum_data[, unlist(interval)])
    } else {
      data <- rep(min(Re(spectrum_data)), m)
      data[unlist(interval)] <- Re(spectrum_data[unlist(interval)])
    }
  }

  # Apply the peak location search method ('thres' or 'max') on spectra
  if (method == "thres") {
    tmsp_peaks <- apply(data, 1, find_tmsp_peak, c = c, direction = direction)
  } else { # method == "max
    tmsp_peaks <- apply(Re(data), 1, which.max)
  }

  # Shift spectra according to the TMSPpeaks found
  # Depends on the shift_handling
  # TMSPpeaks is a column index
  maxpeak <- max(tmsp_peaks) # max accross spectra
  minpeak <- min(tmsp_peaks) # min accross spectra

  if (shift_handling %in% c("zerofilling", "NAfilling", "cut")) {
    fill <- NA
    if (shift_handling == "zerofilling") {
      fill <- 0
    }

    start <- maxpeak - 1
    end <- minpeak - m

    # ppm values of each interval for the whole
    # spectral range of the spectral matrix
    ppm_scale <- (start:end) * ppm_interval

    # check if ppm_value is in the ppm_scale interval
    if (ppm_value < min(ppm_scale) || ppm_value > max(ppm_scale)) {
      warning(
              "ppm_value = ", ppm_value, " is not in the ppm interval [",
              round(min(ppm_scale), 2), ",", round(max(ppm_scale), 2),
              "], and is set to its default ppm_value 0")
      ppm_value <- 0
    }

    # if ppm_value != 0, ppm_scale is adapted
    ppm_scale <- ppm_scale + ppm_value

    # create the spectral matrix with realigned spectra
    spectrum_data_calib <- matrix(fill,
      nrow = n, ncol = -(end - start) + 1,
      dimnames = list(rownames(spectrum_data), ppm_scale)
    )

    # fills in spectrum_data_calib with shifted spectra
    for (i in 1:n) {
      shift <- (1 - tmsp_peaks[i]) + start
      spectrum_data_calib[i, (1 + shift):(m + shift)] <- spectrum_data[i, ]
    }

    if (shift_handling == "cut") {
      spectrum_data_calib <- as.matrix(stats::na.omit(t(spectrum_data_calib)))
      spectrum_data_calib <- t(spectrum_data_calib)
      base::attr(spectrum_data_calib, "na.action") <- NULL
    }
  } else {
    # circular
    start <- 1 - maxpeak
    end <- m - maxpeak

    ppm_scale <- (start:end) * ppm_interval

    # check if ppm_value in is the ppm_scale interval
    if (ppm_value < min(ppm_scale) || ppm_value > max(ppm_scale)) {
      warning(
        "ppm_value = ", ppm_value, " is not in the ppm interval [",
        round(min(ppm_scale), 2), ",", round(max(ppm_scale), 2), "],
        and is set to its default ppm_value 0"
      )
      ppm_value <- 0
    }

    # if ppm_value != 0, ppm_scale is adapted
    ppm_scale <- ppm_scale + ppm_value

    # create the spectral matrix with realigned spectra
    spectrum_data_calib <- matrix(
      nrow = n, ncol = end - start + 1,
      dimnames = list(rownames(spectrum_data), ppm_scale)
    )

    # fills in spectrum_data_calib with shifted spectra
    for (i in 1:n) {
      shift <- (maxpeak - tmsp_peaks[i])
      spectrum_data_calib[i, (1 + shift):m] <- spectrum_data[i, 1:(m - shift)]
      if (shift > 0) {
        spectrum_data_calib[i, 1:shift] <- spectrum_data[i, (m - shift + 1):m]
      }
    }
  }

  # Plot of the spectra (depending on rowindex_graph)
  ppm <- xstart <- value <- xend <- legend <- NULL # only for R CMD check

  # with the search zone for TMSP and the location of the peaks just found
  if (!is.null(rowindex_graph)) {
    if (range == "window") {
      if (ppm_ir == TRUE) {
        fromto <- fromto_rc
      } else {
        fromto <- list()
        idcol <- as.numeric(colnames(spectrum_data))
        for (i in seq_along(fromto_rc)) {
          fromto[[i]] <- as.numeric(colnames(spectrum_data))[fromto_rc[[i]]]
        }
      }
    } else {
      fromto <- fromto_rc
    }

    # TMSPloc in ppm
    tmsp_loc <- as.numeric(colnames(spectrum_data))[tmsp_peaks[rowindex_graph]]

    # num plot per window
    num_stacked <- 6

    # rectanglar bands of color for the search zone
    rects <- data.frame(
      xstart = sapply(fromto, function(x) x[[1]]),
      xend = sapply(fromto, function(x) x[[2]]),
      legend = "Peak search zone and location"
    )

    # vlines for TMSP peak
    addlines <- data.frame(
                           rowname = rownames(spectrum_data)[rowindex_graph],
                           tmsp_loc)

    nn <- length(rowindex_graph)
    i <- 1
    j <- 1
    plots <- vector(mode = "list", length = ceiling(nn / num_stacked))

    while (i <= nn) {
      last <- min(i + num_stacked - 1, nn)
      melted <- reshape2::melt(Re(spectrum_data[i:last, ]),
        varnames = c("rowname", "ppm")
      )

      plots[[j]] <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_line(data = melted, ggplot2::aes(x = ppm, y = value)) +
        ggplot2::geom_rect(data = rects, ggplot2::aes(
          xmin = xstart, xmax = xend,
          ymin = -Inf, ymax = Inf, fill = legend
        ), alpha = 0.4) +
        ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_vline(
          data = addlines, ggplot2::aes(xintercept = tmsp_loc),
          color = "red", show.legend = TRUE
        ) +
        ggplot2::ggtitle("Peak search zone and location") +
        ggplot2::theme(legend.position = "top",
                       legend.text = ggplot2::element_text())

      if ((melted[1, "ppm"] - melted[(dim(melted)[1]), "ppm"]) > 0) {
        plots[[j]] <- plots[[j]] + ggplot2::scale_x_reverse()
      }

      i <- last + 1
      j <- j + 1
    }
    plots
  }

  # Return the results ----------------------------------------------
  spectrum_data <- end_treatment("InternalReferencing",
                                 begin_info, spectrum_data_calib)

  if (is.null(plots)) {
    return(spectrum_data)
  } else {
    return(list(spectrum_data = spectrum_data, plots = plots))
  }
}

## ====================================================
# ZeroOrderPhaseCorrection
## ====================================================

zero_order_phase_correction <- function(spectrum_data,
                                        type_zopc = c("rms", "manual", "max"),
                                        plot_rms = NULL, return_angle = FALSE,
                                        create_window = TRUE,
                                        angle = NULL, plot_spectra = FALSE,
                                        ppm_zopc = TRUE,
                                        exclude_zopc = list(c(5.1, 4.5))) {

  # Data initialisation and checks
  # Entry arguments definition:
  # plot_rms : graph of rms criterion return_angle : if TRUE, returns avector of
  # optimal angles createWindow : for plot_rms plots angle :
  # If angle is not NULL, spectra are rotated according to the angle
  # vector values plot_spectra : if TRUE, plot rotated spectra
  begin_info <- begin_treatment("ZeroOrderPhaseCorrection", spectrum_data)
  spectrum_data <- begin_info[["signal_data"]]
  n <- nrow(spectrum_data)
  m <- ncol(spectrum_data)

  rnames <- rownames(spectrum_data)

  # Check input arguments
  type_zopc <- match.arg(type_zopc)
  check_arg(ppm_zopc, c("bool"))
  check_arg(unlist(exclude_zopc), c("num"), can_be_null = TRUE)

  # type_zopc in c("max", "rms") -----------------------------------------
  if (type_zopc %in% c("max", "rms")) {
    # angle is found by optimization
    # rms function to be optimised
    rms <- function(ang, y, meth = c("max", "rms")) {
      roty <- y * exp(complex(real = 0, imaginary = ang)) # spectrum rotation
      rey <- Re(roty)

      if (meth == "rms") {
        rey_pos <- rey[rey >= 0] # select positive intensities
        pos_ss <- sum((rey_pos)^2, na.rm = TRUE) # SS for positive intensities
        ss <- sum((rey)^2, na.rm = TRUE) #  SS for all intensities
        # criterion : SS for positive values / SS for all intensities
        return(pos_ss / ss)
      } else {
        maxi <- max(rey, na.rm = TRUE)
        return(maxi)
      }
    }

    # Define the interval where to search for (by defining data)
    if (is.null(exclude_zopc)) {
      data <- spectrum_data
    } else {
      # if ppm_zopc == TRUE, then exclude_zopc is in the colnames values,
      # else, in the column index
      if (ppm_zopc == TRUE) {
        colindex <- as.numeric(colnames(spectrum_data))
      } else {
        colindex <- 1:m
      }

      # Second check for the argument exclude_zopc
      diff <- diff(unlist(exclude_zopc))[seq_along(diff(unlist(exclude_zopc)))
                                         %% 2 != 0]
      for (i in seq_along(diff)) {
        if (ppm_zopc == TRUE && diff[i] >= 0) {
          stop(paste("Invalid region removal because from <= to in ppm_zopc"))
        } else if (ppm_zopc == FALSE && diff[i] <= 0) {
          stop(paste("Invalid region removal because from
                     >= to in column index"))
        }
      }

      interval <- vector("list", length(exclude_zopc))
      for (i in seq_along(exclude_zopc)) {
        interval[[i]] <- index_interval(colindex,
          from = exclude_zopc[[i]][1],
          to = exclude_zopc[[i]][2], inclusive = TRUE
        )
      }

      vector <- rep(1, m)
      vector[unlist(interval)] <- 0
      if (n > 1) {
        # Cropped_spectrum
        data <- sweep(spectrum_data, MARGIN = 2, FUN = "*", vector)
      } else {
        data <- spectrum_data * vector
      } # Cropped_spectrum
    }

    # angles computation
    angles <- c()
    for (k in 1:n) {
      # The function is rms is periodic (period 2pi) and it seems that
      # there is a phase x such that rms is unimodal (i.e. decreasing)
      # then increasing on the interval [x; x+2pi]. However, if we do
      # the optimization for example on [x-pi; x+pi], instead of being
      # decreasing then increasing, it might be increasing then decreasing
      # in which case optimize, thinking it is a valley will have to choose
      # between the left or the right of this hill and if it chooses wrong,
      # it will end up at like x-pi while the minimum is close to x+pi.

      # Supposing that rms is unimodal, the classical 1D unimodal optimization
      # will work in either [-pi;pi] or [0;2pi] (this is not easy to be
      # convinced by that I agree) and we can check which one it is simply
      # by the following trick

      f0 <- rms(0, data[k, ], meth = type_zopc)
      fpi <- rms(pi, data[k, ], meth = type_zopc)
      if (f0 < fpi) {
        interval <- c(-pi, pi)
      } else {
        interval <- c(0, 2 * pi)
      }

      # graphs of rms criteria
      # rms should not plot anything now, only when called by optimize
      debug_plot <- FALSE
      if (!is.null(plot_rms) && rnames[k] %in% plot_rms) {
        x <- seq(min(interval), max(interval), length.out = 100)
        y <- rep(1, 100)
        for (K in (1:100)) {
          y[K] <- rms(x[K], data[k, ], meth = type_zopc)
        }
        if (create_window == TRUE) {
          grDevices::dev.new(noRStudioGD = FALSE)
        }
        graphics::plot(x, y,
          main = paste(
            "Criterion maximization \n",
            rownames(data)[k]
          ), ylim = c(0, 1.1),
          ylab = "positiveness criterion", xlab = "angle "
        )
        debug_plot <- TRUE
      }

      # Best angle
      best <- stats::optimize(rms,
        interval = interval, maximum = TRUE,
        y = data[k, ], meth = type_zopc
      )
      ang <- best[["maximum"]]

      if (debug_plot) {
        graphics::abline(v = ang, col = "black")
        graphics::text(x = (ang + 0.1 * ang), y = (y[ang] - 0.1 * y[ang]),
                       labels = round(ang, 3))
      }

      # spectrum rotation
      spectrum_data[k, ] <- spectrum_data[k, ] * exp(complex(real = 0,
                                                             imaginary = ang))
      angles <- c(angles, ang)
    }
  } else {
    # type_zopc is "manual"
    # if angles is already specified and no optimisation is needed
    angles <- angle

    if (!is.vector(angle)) {
      stop("angle is not a vector")
    }

    if (!is.numeric(angle)) {
      stop("angle is not a numeric")
    }

    if (length(angle) != n) {
      stop(paste("angle has length", length(angle), "and there are", n,
                 "spectra to rotate."))
    }
    for (k in 1:n) {
      spectrum_data[k, ] <- spectrum_data[k, ] *
        exp(complex(real = 0, imaginary = -angle[k]))
    }
  }

  #  Draw spectra
  if (plot_spectra == TRUE) {
    nn <- ceiling(n / 4)
    i <- 1
    for (k in 1:nn) {
      if (create_window == TRUE) {
        grDevices::dev.new(noRStudioGD = FALSE)
      }
      graphics::par(mfrow = c(4, 2))
      while (i <= n) {
        last <- min(i + 4 - 1, n)
        graphics::plot(Re(spectrum_data[i, ]),
          type = "l", ylab = "intensity",
          xlab = "Index", main = paste0(rownames(spectrum_data)[i],
                                        " - Real part")
        )
        graphics::plot(Im(spectrum_data[i, ]),
          type = "l", ylab = "intensity",
          xlab = "Index", main = paste0(rownames(spectrum_data)[i],
                                        " - Imaginary part")
        )
        i <- i + 1
      }
      i <- last + 1
    }
  }

  # Data finalisation ----------------------------------------------
  spectrum_data <- end_treatment("ZeroOrderPhaseCorrection",
                                 begin_info, spectrum_data)
  if (return_angle) {
    return(list(spectrum_data = spectrum_data, angles = angles))
  } else {
    return(spectrum_data)
  }
}


## ====================================================
# Baseline Correction
## ====================================================
baseline_correction <- function(spectrum_data, ptw_bc = TRUE, max_iter = 42,
                                lambda_bc = 1e+07, p_bc = 0.05, eps = 1e-08,
                                ppm_bc = TRUE, exclude_bc = list(c(5.1, 4.5)),
                                return_baseline = FALSE) {

  # Data initialisation ----------------------------------------------
  begin_info <- begin_treatment("BaselineCorrection",
                                spectrum_data, force_real = TRUE)
  spectrum_data <- begin_info[["signal_data"]]
  p <- p_bc
  lambda <- lambda_bc
  n <- dim(spectrum_data)[1]
  m <- dim(spectrum_data)[2]

  # Data check
  check_arg(ptw_bc, c("bool"))
  check_arg(max_iter, c("int", "pos"))
  check_arg(lambda, c("num", "pos0"))
  check_arg(p_bc, c("num", "pos0"))
  check_arg(eps, c("num", "pos0"))
  check_arg(return_baseline, c("bool"))
  check_arg(ppm_bc, c("bool"))
  check_arg(unlist(exclude_bc), c("num"), can_be_null = TRUE)

  # Define the interval where to search for (by defining Data)
  if (is.null(exclude_bc)) {
    exclude_index <- NULL
  } else {
    # if ppm_bc == TRUE, then exclude_bc is in the colnames values, else,
    # in the column index
    if (ppm_bc == TRUE) {
      colindex <- as.numeric(colnames(spectrum_data))
    } else {
      colindex <- 1:m
    }

    interval <- vector("list", length(exclude_bc))
    for (i in seq_len(exclude_bc)) {
      interval[[i]] <- index_interval(colindex,
        from = exclude_bc[[i]][1],
        to = exclude_bc[[i]][2], inclusive = TRUE
      )
    }
    exclude_index <- unlist(interval)
  }

  # Baseline Correction implementation definition
  # 2 Ways: either use the function asysm from the ptw package or by
  # built-in functions
  if (ptw_bc) {
    asysm <- ptw::asysm
  } else {
    difsmw <- function(y, lambda, w, d) {
      # Weighted smoothing with a finite difference penalty cf Eilers, 2003.
      # (A perfect smoother)
      # y: signal to be smoothed
      # lambda: smoothing parameter
      # w: weights (use0 zeros for missing values)
      # d: order of differences in penalty (generally 2)
      m <- length(y)
      w_diag <- Matrix::Diagonal(x = w)
      e <- Matrix::Diagonal(m)
      d_diff <- Matrix::diff(e, differences = d)
      c_chol <- Matrix::chol(w_diag + lambda * t(d_diff) %*% d_diff)
      x <- Matrix::solve(c_chol, Matrix::solve(t(c_chol), w * y))
      as.numeric(x)
    }
    asysm <- function(y, lambda, p, eps, exclude_index) {
      # Baseline estimation with asymmetric least squares
      # y: signal
      # lambda: smoothing parameter (generally 1e5 to 1e8)
      # p: asymmetry parameter (generally 0.001)
      # d: order of differences in penalty (generally 2)
      # eps: 1e-8 in ptw package
      m <- length(y)
      w <- rep(1, m)
      i <- 1
      repeat {
        z <- difsmw(y, lambda, w, d = 2)
        w0 <- w
        p_vect <- rep((1 - p), m) # if y <= z + eps
        p_vect[y > z + eps | y < 0] <- p # if y > z + eps | y < 0
        if (!is.null(exclude_index)) {
          p_vect[exclude_index] <- 0 # if exclude area
        }

        w <- p_vect

        if (sum(abs(w - w0)) == 0) {
          break
        }
        i <- i + 1
        if (i > max_iter) {
          warning("cannot find Baseline estimation in asysm")
          break
        }
      }
      z
    }
  }

  # Baseline estimation
  baseline <- matrix(NA, nrow = nrow(spectrum_data), ncol = ncol(spectrum_data))

  if (ptw_bc) {
    baseline <- apply(spectrum_data, 1, asysm,
                      lambda = lambda, p = p,
                      eps = eps)
  } else {
    baseline <- apply(spectrum_data, 1, asysm,
                      lambda = lambda, p = p,
                      eps = eps, exclude_index = exclude_index)
  }

  spectrum_data <- spectrum_data - t(baseline)

  # Data finalisation ----------------------------------------------
  spectrum_data <- end_treatment("BaselineCorrection", begin_info,
                                 spectrum_data)

  if (return_baseline) {
    return(list(spectrum_data = spectrum_data, baseline = baseline))
  } else {
    return(spectrum_data)
  }
}



## ====================================================
# negative_values_zeroing
## ====================================================

negative_values_zeroing <- function(spectrum_data) {
  # Data initialisation and checks -------------------------------
  begin_info <- begin_treatment("negative_values_zeroing", spectrum_data,
                                force_real = TRUE)
  spectrum_data <- begin_info[["signal_data"]]

  # negative_values_zeroing ----------------------------------------------
  spectrum_data[spectrum_data < 0] <- 0

  # Data finalisation ----------------------------------------------
  end_treatment("negative_values_zeroing", begin_info, spectrum_data)
}
