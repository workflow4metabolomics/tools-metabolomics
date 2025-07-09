## ==========================
# Internal functions
## ==========================

# beginTreatment
beginTreatment <- function(name, Signal_data = NULL, Signal_info = NULL,
                           force.real = FALSE) {
  cat("Begin", name, "\n")


  # Formatting the Signal_data and Signal_info -----------------------

  vec <- is.vector(Signal_data)
  if (vec) {
    Signal_data <- vec2mat(Signal_data)
  }
  if (is.vector(Signal_info)) {
    Signal_info <- vec2mat(Signal_info)
  }
  if (!is.null(Signal_data)) {
    if (!is.matrix(Signal_data)) {
      stop("Signal_data is not a matrix.")
    }
    if (!is.complex(Signal_data) && !is.numeric(Signal_data)) {
      stop("Signal_data contains non-numerical values.")
    }
  }
  if (!is.null(Signal_info) && !is.matrix(Signal_info)) {
    stop("Signal_info is not a matrix.")
  }


  Original_data <- Signal_data

  # Extract the real part of the spectrum ---------------------------

  if (force.real) {
    if (is.complex(Signal_data)) {
      Signal_data <- Re(Signal_data)
    } else {
      # The signal is numeric Im(Signal_data) is zero anyway so let's avoid
      # using complex(real=...,imaginary=0) which would give a complex signal
      # in endTreatment()
      force.real <- FALSE
    }
  }


  # Return the formatted data and metadata entries --------------------

  return(list(
    start = proc.time(), vec = vec, force.real = force.real,
    Original_data = Original_data, Signal_data = Signal_data, Signal_info = Signal_info
  ))
}

# endTreatment
endTreatment <- function(name, begin_info, Signal_data) {
  # begin_info: object outputted from beginTreatment


  # Formatting the entries and printing process time -----------------------
  end_time <- proc.time() # record it as soon as possible
  start_time <- begin_info[["start"]]
  delta_time <- end_time - start_time
  delta <- delta_time[]
  cat("End", name, "\n")
  cat(
    "It lasted", round(delta["user.self"], 3), "s user time,", round(delta["sys.self"], 3),
    "s system time and", round(delta["elapsed"], 3), "s elapsed time.\n"
  )


  if (begin_info[["force.real"]]) {
    # The imaginary part is left untouched
    i <- complex(real = 0, imaginary = 1)
    Signal_data <- Signal_data + i * Im(begin_info[["Original_data"]])
  }

  if (begin_info[["vec"]]) {
    Signal_data <- Signal_data[1, ]
  }

  # Return the formatted data and metadata entries --------------------
  return(Signal_data)
}

# checkArg
checkArg <- function(arg, checks, can.be.null = FALSE) {
  check.list <- list(
    bool = c(is.logical, "a boolean"),
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
    if (!can.be.null) {
      stop(deparse(substitute(arg)), " is null.")
    }
  } else {
    if (is.matrix(arg)) {
      stop(deparse(substitute(arg)), " is not scalar.")
    }
    for (c in checks) {
      if (!check.list[[c]][[1]](arg)) {
        stop(deparse(substitute(arg)), " is not ", check.list[[c]][[2]], ".")
      }
    }
  }
}

# getArg
getArg <- function(arg, info, argname, can.be.absent = FALSE) {
  if (is.null(arg)) {
    start <- paste("impossible to get argument", argname, "it was not given directly and")
    if (!is.matrix(info)) {
      stop(paste(start, "the info matrix was not given"))
    }
    if (!(argname %in% colnames(info))) {
      if (can.be.absent) {
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

# binarySearch
binarySearch <- function(a, target, lower = TRUE) {
  # search the index i in a such that a[i] == target
  # if it doesn't exists and lower, it searches the closer a[i] such that a[i] < target
  # if !lower, it seraches the closer a[i] such that a[i] > target
  # a should be monotone but can be increasing or decreasing

  # if a is increasing INVARIANT: a[amin] < target < a[amax]
  N <- length(a)
  if ((a[N] - target) * (a[N] - a[1]) <= 0) {
    return(N)
  }
  if ((a[1] - target) * (a[N] - a[1]) >= 0) {
    return(1)
  }
  amin <- 1
  amax <- N
  while (amin + 1 < amax) {
    amid <- floor((amin + amax) / 2)
    if ((a[amid] - target) * (a[amax] - a[amid]) < 0) {
      amin <- amid
    } else if ((a[amid] - target) * (a[amax] - a[amid]) > 0) {
      amax <- amid
    } else {
      # a[amid] == a[amax] or a[amid] == target In both cases, a[amid] ==
      # target
      return(amid)
    }
  }
  if (xor(lower, a[amin] > a[amax])) {
    # (lower && a[amin] < a[amax]) || (!lower && a[min] > a[max])
    # If increasing and we want the lower, we take amin
    # If decreasing and we want the bigger, we take amin too
    return(amin)
  } else {
    return(amax)
  }
}

# Interpol
Interpol <- function(t, y) {
  # y: sample
  # t : warping function

  m <- length(y)
  # t <= m-1
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
  return(matrix(vec, nrow = 1, dimnames = list(c(1), names(vec))))
}

# binarySearch
binarySearch <- function(a, target, lower = TRUE) {
  # search the index i in a such that a[i] == target
  # if it doesn't exists and lower, it searches the closer a[i] such that a[i] < target
  # if !lower, it seraches the closer a[i] such that a[i] > target
  # a should be monotone but can be increasing or decreasing

  # if a is increasing INVARIANT: a[amin] < target < a[amax]
  N <- length(a)
  if ((a[N] - target) * (a[N] - a[1]) <= 0) {
    return(N)
  }
  if ((a[1] - target) * (a[N] - a[1]) >= 0) {
    return(1)
  }
  amin <- 1
  amax <- N
  while (amin + 1 < amax) {
    amid <- floor((amin + amax) / 2)
    if ((a[amid] - target) * (a[amax] - a[amid]) < 0) {
      amin <- amid
    } else if ((a[amid] - target) * (a[amax] - a[amid]) > 0) {
      amax <- amid
    } else {
      # a[amid] == a[amax] or a[amid] == target In both cases, a[amid] ==
      # target
      return(amid)
    }
  }
  if (xor(lower, a[amin] > a[amax])) {
    # (lower && a[amin] < a[amax]) || (!lower && a[min] > a[max])
    # If increasing and we want the lower, we take amin
    # If decreasing and we want the bigger, we take amin too
    return(amin)
  } else {
    return(amax)
  }
}


# indexInterval
indexInterval <- function(a, from, to, inclusive = TRUE) {
  # If inclusive and from <= to, we need to take the lower
  # If not inclusive and from > to, we need to take the lower too
  lowerFrom <- xor(inclusive, from > to)
  fromIndex <- binarySearch(a, from, lowerFrom)
  toIndex <- binarySearch(a, to, !lowerFrom)
  return(fromIndex:toIndex)
}



## ==========================
# GroupDelayCorrection
## ==========================
GroupDelayCorrection <- function(Fid_data, Fid_info = NULL, group_delay = NULL) {
  # Data initialisation and checks ----------------------------------------------

  begin_info <- beginTreatment("GroupDelayCorrection", Fid_data, Fid_info)
  Fid_data <- begin_info[["Signal_data"]]
  dimension_names <- dimnames(Fid_data)
  Fid_info <- begin_info[["Signal_info"]]
  checkArg(group_delay, c("num", "pos0"), can.be.null = TRUE)
  # if Fid_info and group_delay are NULL, getArg will generate an error

  group_delay <- getArg(group_delay, Fid_info, "GRPDLY", can.be.absent = TRUE)

  if (is.null(group_delay)) {
    # See DetermineBrukerDigitalFilter.m in matNMR MATLAB library
    group_delay_matrix <- matrix(
      c(
        44.75, 46, 46.311, 33.5, 36.5, 36.53, 66.625,
        48, 47.87, 59.0833, 50.1667, 50.229, 68.5625, 53.25, 53.289, 60.375,
        69.5, 69.551, 69.5313, 72.25, 71.6, 61.0208, 70.1667, 70.184, 70.0156,
        72.75, 72.138, 61.3438, 70.5, 70.528, 70.2578, 73, 72.348, 61.5052, 70.6667,
        70.7, 70.3789, 72.5, 72.524, 61.5859, 71.3333, NA, 70.4395, 72.25, NA,
        61.6263, 71.6667, NA, 70.4697, 72.125, NA, 61.6465, 71.8333, NA, 70.4849,
        72.0625, NA, 61.6566, 71.9167, NA, 70.4924, 72.0313, NA
      ),
      nrow = 21,
      ncol = 3, byrow = TRUE, dimnames = list(
        c(
          2, 3, 4, 6, 8, 12, 16, 24,
          32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048
        ),
        c(10, 11, 12)
      )
    )
    decim <- Fid_info[1, "DECIM"]
    dspfvs <- Fid_info[1, "DSPFVS"]
    if (!(toString(decim) %in% rownames(group_delay_matrix))) {
      stop(paste("Invalid DECIM", decim, "it should be one of", rownames(group_delay_matrix)))
    }
    if (!(toString(dspfvs) %in% colnames(group_delay_matrix))) {
      stop(paste("Invalid DSPFVS", dspfvs, "it should be one of", colnames(group_delay_matrix)))
    }
    group_delay <- group_delay_matrix[toString(decim), toString(dspfvs)]
    if (is.na(group_delay)) {
      stop(paste("Invalid DECIM", decim, "for DSPFVS", dspfvs))
    }
  }
  m <- ncol(Fid_data)
  n <- nrow(Fid_data)

  # GroupDelayCorrection ----------------------------------------------

  # We do the shifting in the Fourier domain because the shift can be non-integer.
  # That way we automatically have the circular behaviour of the shift and the
  # interpolation if it is non-integer.

  Spectrum <- t(stats::mvfft(t(Fid_data)))

  # Spectrum <- FourierTransform(Fid_data, Fid_info)
  p <- ceiling(m / 2)
  new_index <- c((p + 1):m, 1:p)
  Spectrum <- Spectrum[, new_index]
  Spectrum <- matrix(data = Spectrum, ncol = m, nrow = n)

  Omega <- (0:(m - 1)) / m
  i <- complex(real = 0, imaginary = 1)

  if (n > 1) {
    Spectrum <- sweep(Spectrum, MARGIN = 2, exp(i * group_delay * 2 * pi * Omega), `*`)
    Spectrum <- Spectrum[, new_index]
  } else {
    Spectrum <- Spectrum * exp(i * group_delay * 2 * pi * Omega)
    Spectrum <- Spectrum[new_index]
    Spectrum <- matrix(data = Spectrum, ncol = m, nrow = n)
  }


  Fid_data <- t(stats::mvfft(t(Spectrum), inverse = TRUE)) / m
  colnames(Fid_data) <- dimension_names[[2]]
  rownames(Fid_data) <- dimension_names[[1]]

  # Data finalisation ----------------------------------------------

  return(endTreatment("GroupDelayCorrection", begin_info, Fid_data))
}

## ==========================
# SolventSuppression
## ==========================
SolventSuppression <- function(Fid_data, lambda.ss = 1e+06, ptw.ss = TRUE,
                               plotSolvent = F, returnSolvent = F) {
  # Data initialisation and checks ----------------------------------------------

  begin_info <- beginTreatment("SolventSuppression", Fid_data)
  Fid_data <- begin_info[["Signal_data"]]
  checkArg(ptw.ss, c("bool"))
  checkArg(lambda.ss, c("num", "pos0"))


  # difsm function definition for the smoother -----------------------------------

  if (ptw.ss) {
    # Use of the function in ptw that smoothes signals with a finite difference
    # penalty of order 2
    difsm <- ptw::difsm
  } else {
    # Or manual implementation based on sparse matrices for large data series (cf.
    # Eilers, 2003. 'A perfect smoother')
    difsm <- function(y, d = 2, lambda) {
      m <- length(y)
      # Sparse identity matrix m x m
      E <- Matrix::Diagonal(m)
      D <- Matrix::diff(E, differences = d)
      A <- E + lambda.ss * Matrix::t(D) %*% D
      # base::chol does not take into account that A is sparse and is extremely slow
      C <- Matrix::chol(A)
      x <- Matrix::solve(C, Matrix::solve(Matrix::t(C), y))
      return(as.numeric(x))
    }
  }

  # Solvent Suppression ----------------------------------------------

  n <- dim(Fid_data)[1]
  if (returnSolvent) {
    SolventRe <- Fid_data
    SolventIm <- Fid_data
  }
  for (i in 1:n) {
    FidRe <- Re(Fid_data[i, ])
    FidIm <- Im(Fid_data[i, ])
    solventRe <- difsm(y = FidRe, lambda = lambda.ss)
    solventIm <- difsm(y = FidIm, lambda = lambda.ss)

    if (plotSolvent) {
      m <- length(FidRe)
      graphics::plot(1:m, FidRe, type = "l", col = "red")
      graphics::lines(1:m, solventRe, type = "l", col = "blue")
      graphics::plot(1:m, FidIm, type = "l", col = "red")
      graphics::lines(1:m, solventIm, type = "l", col = "blue")
    }
    FidRe <- FidRe - solventRe
    FidIm <- FidIm - solventIm
    Fid_data[i, ] <- complex(real = FidRe, imaginary = FidIm)
    if (returnSolvent) {
      SolventRe[i, ] <- solventRe
      SolventIm[i, ] <- solventIm
    }
  }


  # Data finalisation ----------------------------------------------

  Fid_data <- endTreatment("SolventSuppression", begin_info, Fid_data)
  if (returnSolvent) {
    return(list(Fid_data = Fid_data, SolventRe = SolventRe, SolventIm = SolventIm))
  } else {
    return(Fid_data)
  }
}


## ==========================
# Apodization
# =============================
Apodization <- function(Fid_data, Fid_info = NULL, DT = NULL,
                        type.apod = c(
                          "exp", "cos2", "blockexp", "blockcos2",
                          "gauss", "hanning", "hamming"
                        ), phase = 0, rectRatio = 1 / 2,
                        gaussLB = 1, expLB = 1, plotWindow = F, returnFactor = F) {
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("Apodization", Fid_data, Fid_info)
  Fid_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]
  # Data check
  type.apod <- match.arg(type.apod)
  checkArg(DT, c("num", "pos"), can.be.null = TRUE)
  checkArg(phase, c("num"))

  # Apodization ----------------------------------------------
  DT <- getArg(DT, Fid_info, "DT") # Dwell Time
  m <- ncol(Fid_data)
  t <- (1:m) * DT # Time
  rectSize <- ceiling(rectRatio * m)
  gaussLB <- (gaussLB / (sqrt(8 * log(2))))
  # Define the types of apodization:
  switch(type.apod,
    exp = {
      # exponential
      Factor <- exp(-expLB * t)
    },
    cos2 = {
      # cos^2
      c <- cos((1:m) * pi / (2 * m) - phase * pi / 2)
      Factor <- c * c
    },
    blockexp = {
      # block and exponential
      Factor <- c(rep.int(1, rectSize), rep.int(0, m - rectSize))
      # | rectSize | 1 ___________ | \ 0 \____
      Factor[(rectSize + 1):m] <- exp(-expLB * t[1:(m - rectSize)])
    },
    blockcos2 = {
      # block and cos^2
      Factor <- c(rep.int(1, rectSize), rep.int(0, m - rectSize))
      c <- cos((1:(m - rectSize)) * pi / (2 * (m - rectSize)))
      Factor[(rectSize + 1):m] <- c * c
    },
    gauss = {
      # gaussian
      Factor <- exp(-(gaussLB * t)^2 / 2)
      Factor <- Factor / max(Factor)
    },
    hanning = {
      # Hanning
      Factor <- 0.5 + 0.5 * cos((1:m) * pi / m - phase * pi)
    },
    hamming = {
      # Hamming
      Factor <- 0.54 + 0.46 * cos((1:m) * pi / m - phase * pi)
    }
  )
  if (plotWindow) {
    graphics::plot(1:m, Factor, "l")
    # dev.off() # device independent, it is the responsability of the
    # caller to do it
  }
  # Apply the apodization factor on the spectra
  Fid_data <- sweep(Fid_data, MARGIN = 2, Factor, `*`)

  # Data finalisation ----------------------------------------------
  Fid_data <- endTreatment("Apodization", begin_info, Fid_data)
  if (returnFactor) {
    return(list(Fid_data = Fid_data, Factor = Factor))
  } else {
    return(Fid_data)
  }
}


## ====================================================
# FourierTransform
## ====================================================


# fftshift1D2D
fftshift1D2D <- function(x) {
  vec <- F
  if (is.vector(x)) {
    x <- vec2mat(x)
    vec <- T
  }
  m <- dim(x)[2]
  p <- ceiling(m / 2)
  new_index <- c((p + 1):m, 1:p)
  y <- x[, new_index, drop = vec]
}

# FourierTransform
FourierTransform <- function(Fid_data, Fid_info = NULL, SW_h = NULL, SW = NULL, O1 = NULL, reverse.axis = TRUE) {
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("FourierTransform", Fid_data, Fid_info)
  Fid_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]

  m <- ncol(Fid_data)
  n <- nrow(Fid_data)

  if (is.null(SW_h)) {
    SW_h <- getArg(SW_h, Fid_info, "SW_h")
  }

  if (is.null(SW)) {
    SW <- getArg(SW, Fid_info, "SW") # Sweep Width in ppm (semi frequency scale in ppm)
  }


  if (is.null(O1)) {
    O1 <- getArg(O1, Fid_info, "O1")
  }


  checkArg(reverse.axis, c("bool"))

  # Fourier Transformation ----------------------------------------------
  # mvfft does the unnormalized fourier transform (see ?mvfft), so we need divide
  # by m.  It does not matter a lot in our case since the spectrum will be
  # normalized.

  # FT
  RawSpect_data <- fftshift1D2D(t(stats::mvfft(t(Fid_data))))
  # recover the frequencies values
  f <- ((0:(m - 1)) - floor(m / 2)) * Fid_info[1, "SW_h"] / (m - 1)

  if (reverse.axis == TRUE) {
    revind <- rev(1:m)
    RawSpect_data <- RawSpect_data[, revind] # reverse the spectrum
  }

  RawSpect_data <- matrix(RawSpect_data, nrow = n, ncol = m)
  colnames(RawSpect_data) <- f
  rownames(RawSpect_data) <- rownames(Fid_data)

  # PPM conversion ----------------------------------------------

  # The Sweep Width has to be the same since the column names are the same

  ppmInterval <- SW / (m - 1)

  O1index <- round((m + 1) / 2 + O1 * (m - 1) / SW_h)

  end <- O1index - m
  start <- O1index - 1
  ppmScale <- (start:end) * ppmInterval
  RawSpect_data <- matrix(RawSpect_data,
    nrow = n, ncol = -(end - start) + 1, dimnames =
      list(rownames(RawSpect_data), ppmScale)
  )


  # Data finalisation ----------------------------------------------
  return(endTreatment("FourierTransform", begin_info, RawSpect_data))
}

## ====================================================
#   InternalReferencing
## ====================================================

InternalReferencing <- function(Spectrum_data, Fid_info, method = c("max", "thres"),
                                range = c("nearvalue", "all", "window"), ppm.value = 0,
                                direction = "left", shiftHandling = c(
                                  "zerofilling", "cut",
                                  "NAfilling", "circular"
                                ), c = 2, pc = 0.02, fromto.RC = NULL,
                                ppm.ir = TRUE, rowindex_graph = NULL) {
  # Data initialisation and checks ----------------------------------------------

  begin_info <- beginTreatment("InternalReferencing", Spectrum_data, Fid_info)
  Spectrum_data <- begin_info[["Signal_data"]]
  Fid_info <- begin_info[["Signal_info"]]


  ######## Check input arguments

  range <- match.arg(range)
  shiftHandling <- match.arg(shiftHandling)
  method <- match.arg(method)
  plots <- NULL

  checkArg(ppm.ir, c("bool"))
  checkArg(unlist(fromto.RC), c("num"), can.be.null = TRUE)
  checkArg(pc, c("num"))
  checkArg(ppm.value, c("num"))
  checkArg(rowindex_graph, "num", can.be.null = TRUE)

  # fromto.RC : if range == "window",
  # fromto.RC defines the spectral window where to search for the peak
  if (!is.null(fromto.RC)) {
    diff <- diff(unlist(fromto.RC))[1:length(diff(unlist(fromto.RC))) %% 2 != 0]
    for (i in 1:length(diff)) {
      if (diff[i] >= 0) {
        fromto <- c(fromto.RC[[i]][2], fromto.RC[[i]][1])
        fromto.RC[[i]] <- fromto
      }
    }
  }


  # findTMSPpeak function ----------------------------------------------
  # If method == "tresh", findTMSPpeak will find the position of the first
  # peak (from left or right) which is higher than a predefined threshold
  # and is computed as: c*(cumulated_mean/cumulated_sd)
  findTMSPpeak <- function(ft, c = 2, direction = "left") {
    ft <- Re(ft) # extraction de la partie réelle
    N <- length(ft)
    if (direction == "left") {
      newindex <- rev(1:N)
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
      # recherche dans les 1% de points suivants du max trouve pour etre au sommet du
      # pic
      d <- which.max(ft[v:(v + N * 0.01)])
      new.peak <- v + d - 1 # nouveau pic du TMSP si d > 0

      if (names(which.max(ft[v:(v + N * 0.01)])) != names(which.max(ft[v:(v + N * 0.03)]))) {
        # recherche dans les 3% de points suivants du max trouve pour eviter un faux
        # positif
        warning("the TMSP peak might be located further away, increase the threshold to check.")
      }
      return(new.peak)
    }
  }


  # Define the search zone  ----------------------------------------

  n <- nrow(Spectrum_data)
  m <- ncol(Spectrum_data)

  # The Sweep Width (SW) has to be the same since the column names are the same
  SW <- Fid_info[1, "SW"] # Sweep Width in ppm
  ppmInterval <- SW / (m - 1) # size of a ppm interval

  # range: How the search zone is defined ("all", "nearvalue" or "window")
  if (range == "all") {
    Data <- Spectrum_data
  } else { # range = "nearvalue" or "window"
    # Need to define colindex (column indexes) to apply indexInterval on it

    if (range == "nearvalue") {
      fromto.RC <- list(c(-(SW * pc) / 2 + ppm.value, (SW * pc) / 2 + ppm.value)) # automatic fromto values in ppm
      colindex <- as.numeric(colnames(Spectrum_data))
    } else {
      # range == "window"
      # fromto.RC is already user-defined
      if (ppm.ir == TRUE) {
        colindex <- as.numeric(colnames(Spectrum_data))
      } else {
        colindex <- 1:m
      }
    }

    # index intervals taking into account the different elements in the list fromto.RC
    Int <- vector("list", length(fromto.RC))
    for (i in 1:length(fromto.RC)) {
      Int[[i]] <- indexInterval(colindex,
        from = fromto.RC[[i]][1],
        to = fromto.RC[[i]][2], inclusive = TRUE
      )
    }

    # define Data as the cropped spectrum including the index intervals
    # outside the research zone, the intensities are set to the minimal
    # intensity of the research zone

    if (n > 1) {
      Data <- apply(Re(Spectrum_data[, unlist(Int)]), 1, function(x) rep(min(x), m))
      Data <- t(Data)
      Data[, unlist(Int)] <- Re(Spectrum_data[, unlist(Int)])
    } else {
      Data <- rep(min(Re(Spectrum_data)), m)
      Data[unlist(Int)] <- Re(Spectrum_data[unlist(Int)])
    }
  }


  # Apply the peak location search method ('thres' or 'max') on spectra
  # -----------------------------------------------------------------------

  if (method == "thres") {
    TMSPpeaks <- apply(Data, 1, findTMSPpeak, c = c, direction = direction)
  } else { # method == "max
    TMSPpeaks <- apply(Re(Data), 1, which.max)
  }


  # Shift spectra according to the TMSPpeaks found --------------------------------
  # Depends on the shiftHandling

  # TMSPpeaks is a column index
  maxpeak <- max(TMSPpeaks) # max accross spectra
  minpeak <- min(TMSPpeaks) # min accross spectra


  if (shiftHandling %in% c("zerofilling", "NAfilling", "cut")) {
    fill <- NA
    if (shiftHandling == "zerofilling") {
      fill <- 0
    }

    start <- maxpeak - 1
    end <- minpeak - m

    # ppm values of each interval for the whole spectral range of the spectral matrix
    ppmScale <- (start:end) * ppmInterval

    # check if ppm.value is in the ppmScale interval
    if (ppm.value < min(ppmScale) | ppm.value > max(ppmScale)) {
      warning(
        "ppm.value = ", ppm.value, " is not in the ppm interval [",
        round(min(ppmScale), 2), ",", round(max(ppmScale), 2), "], and is set to its default ppm.value 0"
      )
      ppm.value <- 0
    }

    # if ppm.value != 0, ppmScale is adapted
    ppmScale <- ppmScale + ppm.value

    # create the spectral matrix with realigned spectra
    Spectrum_data_calib <- matrix(fill,
      nrow = n, ncol = -(end - start) + 1,
      dimnames = list(rownames(Spectrum_data), ppmScale)
    )

    # fills in Spectrum_data_calib with shifted spectra
    for (i in 1:n) {
      shift <- (1 - TMSPpeaks[i]) + start
      Spectrum_data_calib[i, (1 + shift):(m + shift)] <- Spectrum_data[i, ]
    }

    if (shiftHandling == "cut") {
      Spectrum_data_calib <- as.matrix(stats::na.omit(t(Spectrum_data_calib)))
      Spectrum_data_calib <- t(Spectrum_data_calib)
      base::attr(Spectrum_data_calib, "na.action") <- NULL
    }
  } else {
    # circular
    start <- 1 - maxpeak
    end <- m - maxpeak

    ppmScale <- (start:end) * ppmInterval

    # check if ppm.value in is the ppmScale interval
    if (ppm.value < min(ppmScale) | ppm.value > max(ppmScale)) {
      warning(
        "ppm.value = ", ppm.value, " is not in the ppm interval [",
        round(min(ppmScale), 2), ",", round(max(ppmScale), 2), "], and is set to its default ppm.value 0"
      )
      ppm.value <- 0
    }

    # if ppm.value != 0, ppmScale is adapted
    ppmScale <- ppmScale + ppm.value

    # create the spectral matrix with realigned spectra
    Spectrum_data_calib <- matrix(
      nrow = n, ncol = end - start + 1,
      dimnames = list(rownames(Spectrum_data), ppmScale)
    )

    # fills in Spectrum_data_calib with shifted spectra
    for (i in 1:n) {
      shift <- (maxpeak - TMSPpeaks[i])
      Spectrum_data_calib[i, (1 + shift):m] <- Spectrum_data[i, 1:(m - shift)]
      if (shift > 0) {
        Spectrum_data_calib[i, 1:shift] <- Spectrum_data[i, (m - shift + 1):m]
      }
    }
  }




  # Plot of the spectra (depending on rowindex_graph) ---------------------------------------------------

  ppm <- xstart <- value <- xend <- Legend <- NULL # only for R CMD check


  # with the search zone for TMSP and the location of the peaks just found
  if (!is.null(rowindex_graph)) {
    if (range == "window") {
      if (ppm.ir == TRUE) {
        fromto <- fromto.RC
      } else {
        fromto <- list()
        idcol <- as.numeric(colnames(Spectrum_data))
        for (i in 1:length(fromto.RC)) {
          fromto[[i]] <- as.numeric(colnames(Spectrum_data))[fromto.RC[[i]]]
        }
      }
    } else {
      fromto <- fromto.RC
    }

    # TMSPloc in ppm
    TMSPloc <- as.numeric(colnames(Spectrum_data))[TMSPpeaks[rowindex_graph]]

    # num plot per window
    num.stacked <- 6

    # rectanglar bands of color for the search zone
    rects <- data.frame(
      xstart = sapply(fromto, function(x) x[[1]]),
      xend = sapply(fromto, function(x) x[[2]]),
      Legend = "Peak search zone and location"
    )

    # vlines for TMSP peak
    addlines <- data.frame(rowname = rownames(Spectrum_data)[rowindex_graph], TMSPloc)

    nn <- length(rowindex_graph)
    i <- 1
    j <- 1
    plots <- vector(mode = "list", length = ceiling(nn / num.stacked))

    while (i <= nn) {
      last <- min(i + num.stacked - 1, nn)

      melted <- reshape2::melt(Re(Spectrum_data[i:last, ]),
        varnames = c("rowname", "ppm")
      )

      plots[[j]] <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_line(
          data = melted,
          ggplot2::aes(x = ppm, y = value)
        ) +
        ggplot2::geom_rect(data = rects, ggplot2::aes(
          xmin = xstart, xmax = xend,
          ymin = -Inf, ymax = Inf, fill = Legend
        ), alpha = 0.4) +
        ggplot2::facet_grid(rowname ~ ., scales = "free_y") +
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_vline(
          data = addlines, ggplot2::aes(xintercept = TMSPloc),
          color = "red", show.legend = TRUE
        ) +
        ggplot2::ggtitle("Peak search zone and location") +
        ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text())



      if ((melted[1, "ppm"] - melted[(dim(melted)[1]), "ppm"]) > 0) {
        plots[[j]] <- plots[[j]] + ggplot2::scale_x_reverse()
      }

      i <- last + 1
      j <- j + 1
    }

    plots
  }


  # Return the results ----------------------------------------------
  Spectrum_data <- endTreatment("InternalReferencing", begin_info, Spectrum_data_calib)

  if (is.null(plots)) {
    return(Spectrum_data)
  } else {
    return(list(Spectrum_data = Spectrum_data, plots = plots))
  }
}

## ====================================================
# ZeroOrderPhaseCorrection
## ====================================================

ZeroOrderPhaseCorrection <- function(Spectrum_data, type.zopc = c("rms", "manual", "max"),
                                     plot_rms = NULL, returnAngle = FALSE, createWindow = TRUE,
                                     angle = NULL, plot_spectra = FALSE,
                                     ppm.zopc = TRUE, exclude.zopc = list(c(5.1, 4.5))) {
  # Data initialisation and checks ----------------------------------------------

  # Entry arguments definition:
  # plot_rms : graph of rms criterion returnAngle : if TRUE, returns avector of
  # optimal angles createWindow : for plot_rms plots angle : If angle is not NULL,
  # spectra are rotated according to the angle vector values
  # plot_spectra : if TRUE, plot rotated spectra



  begin_info <- beginTreatment("ZeroOrderPhaseCorrection", Spectrum_data)
  Spectrum_data <- begin_info[["Signal_data"]]
  n <- nrow(Spectrum_data)
  m <- ncol(Spectrum_data)

  rnames <- rownames(Spectrum_data)

  # Check input arguments
  type.zopc <- match.arg(type.zopc)
  checkArg(ppm.zopc, c("bool"))
  checkArg(unlist(exclude.zopc), c("num"), can.be.null = TRUE)


  # type.zopc in c("max", "rms") -----------------------------------------
  if (type.zopc %in% c("max", "rms")) {
    # angle is found by optimization

    # rms function to be optimised
    rms <- function(ang, y, meth = c("max", "rms")) {
      # if (debug_plot) { graphics::abline(v=ang, col='gray60') }
      roty <- y * exp(complex(real = 0, imaginary = ang)) # spectrum rotation
      Rey <- Re(roty)

      if (meth == "rms") {
        ReyPos <- Rey[Rey >= 0] # select positive intensities
        POSss <- sum((ReyPos)^2, na.rm = TRUE) # SS for positive intensities
        ss <- sum((Rey)^2, na.rm = TRUE) #  SS for all intensities
        return(POSss / ss) # criterion : SS for positive values / SS for all intensities
      } else {
        maxi <- max(Rey, na.rm = TRUE)
        return(maxi)
      }
    }


    # Define the interval where to search for (by defining Data)
    if (is.null(exclude.zopc)) {
      Data <- Spectrum_data
    } else {
      # if ppm.zopc == TRUE, then exclude.zopc is in the colnames values, else, in the column
      # index
      if (ppm.zopc == TRUE) {
        colindex <- as.numeric(colnames(Spectrum_data))
      } else {
        colindex <- 1:m
      }

      # Second check for the argument exclude.zopc
      diff <- diff(unlist(exclude.zopc))[1:length(diff(unlist(exclude.zopc))) %% 2 != 0]
      for (i in 1:length(diff)) {
        if (ppm.zopc == TRUE & diff[i] >= 0) {
          stop(paste("Invalid region removal because from <= to in ppm.zopc"))
        } else if (ppm.zopc == FALSE & diff[i] <= 0) {
          stop(paste("Invalid region removal because from >= to in column index"))
        }
      }


      Int <- vector("list", length(exclude.zopc))
      for (i in 1:length(exclude.zopc)) {
        Int[[i]] <- indexInterval(colindex,
          from = exclude.zopc[[i]][1],
          to = exclude.zopc[[i]][2], inclusive = TRUE
        )
      }

      vector <- rep(1, m)
      vector[unlist(Int)] <- 0
      if (n > 1) {
        Data <- sweep(Spectrum_data, MARGIN = 2, FUN = "*", vector) # Cropped_Spectrum
      } else {
        Data <- Spectrum_data * vector
      } # Cropped_Spectrum
    }


    # angles computation
    Angle <- c()
    for (k in 1:n)
    {
      # The function is rms is periodic (period 2pi) and it seems that there is a phase
      # x such that rms is unimodal (i.e. decreasing then increasing) on the interval
      # [x; x+2pi].  However, if we do the optimization for example on [x-pi; x+pi],
      # instead of being decreasing then increasing, it might be increasing then
      # decreasing in which case optimize, thinking it is a valley will have to choose
      # between the left or the right of this hill and if it chooses wrong, it will end
      # up at like x-pi while the minimum is close to x+pi.

      # Supposing that rms is unimodal, the classical 1D unimodal optimization will
      # work in either [-pi;pi] or [0;2pi] (this is not easy to be convinced by that I
      # agree) and we can check which one it is simply by the following trick

      f0 <- rms(0, Data[k, ], meth = type.zopc)
      fpi <- rms(pi, Data[k, ], meth = type.zopc)
      if (f0 < fpi) {
        interval <- c(-pi, pi)
      } else {
        interval <- c(0, 2 * pi)
      }

      # graphs of rms criteria
      debug_plot <- F # rms should not plot anything now, only when called by optimize
      if (!is.null(plot_rms) && rnames[k] %in% plot_rms) {
        x <- seq(min(interval), max(interval), length.out = 100)
        y <- rep(1, 100)
        for (K in (1:100)) {
          y[K] <- rms(x[K], Data[k, ], meth = type.zopc)
        }
        if (createWindow == TRUE) {
          grDevices::dev.new(noRStudioGD = FALSE)
        }
        graphics::plot(x, y,
          main = paste(
            "Criterion maximization \n",
            rownames(Data)[k]
          ), ylim = c(0, 1.1),
          ylab = "positiveness criterion", xlab = "angle "
        )
        debug_plot <- T
      }

      # Best angle
      best <- stats::optimize(rms,
        interval = interval, maximum = TRUE,
        y = Data[k, ], meth = type.zopc
      )
      ang <- best[["maximum"]]


      if (debug_plot) {
        graphics::abline(v = ang, col = "black")
        graphics::text(x = (ang + 0.1 * ang), y = (y[ang] - 0.1 * y[ang]), labels = round(ang, 3))
      }

      # Spectrum rotation
      Spectrum_data[k, ] <- Spectrum_data[k, ] * exp(complex(real = 0, imaginary = ang))
      Angle <- c(Angle, ang)
    }
  } else {
    # type.zopc is "manual" -------------------------------------------------------
    # if Angle is already specified and no optimisation is needed
    Angle <- angle

    if (!is.vector(angle)) {
      stop("angle is not a vector")
    }

    if (!is.numeric(angle)) {
      stop("angle is not a numeric")
    }

    if (length(angle) != n) {
      stop(paste("angle has length", length(angle), "and there are", n, "spectra to rotate."))
    }
    for (k in 1:n) {
      Spectrum_data[k, ] <- Spectrum_data[k, ] * exp(complex(real = 0, imaginary = -angle[k]))
    }
  }


  #  Draw spectra
  if (plot_spectra == TRUE) {
    nn <- ceiling(n / 4)
    i <- 1
    for (k in 1:nn) {
      if (createWindow == TRUE) {
        grDevices::dev.new(noRStudioGD = FALSE)
      }
      graphics::par(mfrow = c(4, 2))
      while (i <= n) {
        last <- min(i + 4 - 1, n)
        graphics::plot(Re(Spectrum_data[i, ]),
          type = "l", ylab = "intensity",
          xlab = "Index", main = paste0(rownames(Spectrum_data)[i], " - Real part")
        )
        graphics::plot(Im(Spectrum_data[i, ]),
          type = "l", ylab = "intensity",
          xlab = "Index", main = paste0(rownames(Spectrum_data)[i], " - Imaginary part")
        )
        i <- i + 1
      }
      i <- last + 1
    }
  }


  # Data finalisation ----------------------------------------------

  Spectrum_data <- endTreatment("ZeroOrderPhaseCorrection", begin_info, Spectrum_data)
  if (returnAngle) {
    return(list(Spectrum_data = Spectrum_data, Angle = Angle))
  } else {
    return(Spectrum_data)
  }
}


## ====================================================
# Baseline Correction
## ====================================================
BaselineCorrection <- function(Spectrum_data, ptw.bc = TRUE, maxIter = 42,
                               lambda.bc = 1e+07, p.bc = 0.05, eps = 1e-08,
                               ppm.bc = TRUE, exclude.bc = list(c(5.1, 4.5)),
                               returnBaseline = F) {
  # Data initialisation ----------------------------------------------
  begin_info <- beginTreatment("BaselineCorrection", Spectrum_data, force.real = T)
  Spectrum_data <- begin_info[["Signal_data"]]
  p <- p.bc
  lambda <- lambda.bc
  n <- dim(Spectrum_data)[1]
  m <- dim(Spectrum_data)[2]


  # Data check
  checkArg(ptw.bc, c("bool"))
  checkArg(maxIter, c("int", "pos"))
  checkArg(lambda, c("num", "pos0"))
  checkArg(p.bc, c("num", "pos0"))
  checkArg(eps, c("num", "pos0"))
  checkArg(returnBaseline, c("bool"))
  checkArg(ppm.bc, c("bool"))
  checkArg(unlist(exclude.bc), c("num"), can.be.null = TRUE)

  # Define the interval where to search for (by defining Data)
  if (is.null(exclude.bc)) {
    exclude_index <- NULL
  } else {
    # if ppm.bc == TRUE, then exclude.bc is in the colnames values, else, in the column
    # index
    if (ppm.bc == TRUE) {
      colindex <- as.numeric(colnames(Spectrum_data))
    } else {
      colindex <- 1:m
    }

    Int <- vector("list", length(exclude.bc))
    for (i in 1:length(exclude.bc)) {
      Int[[i]] <- indexInterval(colindex,
        from = exclude.bc[[i]][1],
        to = exclude.bc[[i]][2], inclusive = TRUE
      )
    }
    exclude_index <- unlist(Int)
  }

  # Baseline Correction implementation definition ----------------------

  # 2 Ways: either use the function asysm from the ptw package or by
  # built-in functions
  if (ptw.bc) {
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
      W <- Matrix::Diagonal(x = w)
      E <- Matrix::Diagonal(m)
      D <- Matrix::diff(E, differences = d)
      C <- Matrix::chol(W + lambda * t(D) %*% D)
      x <- Matrix::solve(C, Matrix::solve(t(C), w * y))
      return(as.numeric(x))
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
        # w <- p * (y > z + eps | y < 0) + (1 - p) * (y <= z + eps)

        if (sum(abs(w - w0)) == 0) {
          break
        }
        i <- i + 1
        if (i > maxIter) {
          warning("cannot find Baseline estimation in asysm")
          break
        }
      }
      return(z)
    }
  }

  # Baseline estimation ----------------------------------------------
  Baseline <- matrix(NA, nrow = nrow(Spectrum_data), ncol = ncol(Spectrum_data))

  # for (k in 1:n) {
  # Baseline[k, ] <- asysm(y = Spectrum_data[k, ], lambda = lambda, p = p, eps = eps)

  if (ptw.bc) {
    Baseline <- apply(Spectrum_data, 1, asysm,
      lambda = lambda, p = p,
      eps = eps
    )
  } else {
    Baseline <- apply(Spectrum_data, 1, asysm,
      lambda = lambda, p = p,
      eps = eps, exclude_index = exclude_index
    )
  }


  Spectrum_data <- Spectrum_data - t(Baseline)
  # }

  # Data finalisation ----------------------------------------------
  Spectrum_data <- endTreatment("BaselineCorrection", begin_info, Spectrum_data) # FIXME create removeImaginary filter ??

  if (returnBaseline) {
    return(list(Spectrum_data = Spectrum_data, Baseline = Baseline))
  } else {
    return(Spectrum_data)
  }
}



## ====================================================
# NegativeValuesZeroing
## ====================================================

NegativeValuesZeroing <- function(Spectrum_data) {
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("NegativeValuesZeroing", Spectrum_data, force.real = T)
  Spectrum_data <- begin_info[["Signal_data"]]

  # NegativeValuesZeroing ----------------------------------------------
  Spectrum_data[Spectrum_data < 0] <- 0

  # Data finalisation ----------------------------------------------
  return(endTreatment("NegativeValuesZeroing", begin_info, Spectrum_data))
}
