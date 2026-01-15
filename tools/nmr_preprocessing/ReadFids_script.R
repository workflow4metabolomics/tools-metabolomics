##########################################################################
#
#   Read FIDs in Bruker format
#
##########################################################################
# Function to convert a vector to a matrix
vec2mat <- function(vec) {
  return(matrix(vec, nrow = 1, dimnames = list(c(1), names(vec))))
}

# Function to read FID files
# Format of files depends on version of software used for spectra aquisition
# < or > 4.x.x for NEO version
# Read 1D FID using Bruker XWinNMR and TopSpin format.  It is inspired of the
# matNMR matlab library which deals with 2D FID and also other formats
# Based on aquisition parameters stored in the acqus file
read_fid <- function(path) {
  param_file <- file.path(path, "acqus")
  params <- read_params(param_file, c("TD", "BYTORDA", "DIGMOD",
                                      "DECIM", "DSPFVS", "SW_h", "SW", "O1"))

  # Version of the TopSpin software
  line <- readLines(param_file)[1]
  version <- as.numeric(substr(strsplit(line, " ")[[1]][5], 1, 1))
  neo <- (version > 3)

  # Group delay first order phase correction: given directly from version 20
  if (params[["DSPFVS"]] >= 20) {
    grpdly <- read_params(param_file, c("GRPDLY"))
    params[["GRPDLY"]] <- grpdly[["GRPDLY"]]
  }
  TD <- params[["TD"]] # nolint: object_name_linter.

  endianness <- if (params$BYTORDA) {
    "big"
  } else {
    "little"
  }
  if (TD %% 2 != 0) {
    stop(paste("Only even numbers are allowed for size in TD because
               it is complex datawith the real and imaginary part for 
               each element.", "The TD value is in the", param_file, "file"))
  }

  # Interpret params Dwell Time, time between 2 data points in the FID
  params[["DT"]] <- 1 / (2 * params[["SW_h"]])

  # Read fid depend on version of aquisition software
  fid_file <- file.path(path, "fid")
  if (neo) {
    fid_on_disk <- readBin(fid_file, what = "double", n = TD, size = NA_integer,
                           signed = TRUE, endian = .Platform$endian)
  } else {
    fid_on_disk <- readBin(fid_file, what = "int", n = TD, size = 4L,
                           endian = endianness)
  }

  # Real size that is on disk (it should be equal to TD2,
  # except for TopSpin/Bruker
  # (which is our case) according to matNMR as just discussed
  td_on_disk <- length(fid_on_disk)
  if (td_on_disk < TD) {
    warning("Size is smaller than expected, the rest is filled with zero
            so the size is the same for every fid")
    fid_good_size <- sapply(vector("list", length = TD), function(x) 0)
    fid_good_size[1:td_on_disk] <- fid_on_disk
  } else if (td_on_disk > TD) {
    warning("Size is bigger than expected, the rest ignored so the size
            is the same for every fid")
    fid_good_size <- fid_on_disk(1:TD)
  } else {
    fid_good_size <- fid_on_disk
  }

  fid_re_part <- fid_good_size[seq(from = 1, to = TD, by = 2)]
  fid_im_part <- fid_good_size[seq(from = 2, to = TD, by = 2)]
  fid <- complex(real = fid_re_part, imaginary = fid_im_part)

  return(list(fid = fid, params = params))
}

# Function to obtain path to files
get_dirs_containing_fid <- function(path) {
  subdirs <- dir(path, full.names = TRUE)
  if (length(subdirs) > 0) {
    cond <- sapply(subdirs, function(x) {
      content <- dir(x)
      # subdirs must contain fid, acqu and acqus files
      return("fid" %in% content && "acqu" %in% content && "acqus" %in% content)
    })
    subdirs <- subdirs[cond]
  }
  return(subdirs)
}

# Function to xxx
begin_treatment <- function(name, signal_data = NULL, signal_info = NULL,
                            force_real = FALSE) {
  cat("Begin", name, "\n")

  # Formatting the signal_data and signal_info
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
    if (is.vector(signal_info)) {
      signal_info <- vec2mat(signal_info)
    }
  }
  if (!is.null(signal_info) && !is.matrix(signal_info)) {
    stop("signal_info is not a matrix.")
  }

  original_data <- signal_data

  # Extract the real part of the spectrum
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

  # Return the formatted data and metadata entries
  return(list(start = proc.time(), vec = vec, force_real = force_real,
              original_data = original_data, signal_data = signal_data,
              signal_info = signal_info))
}

# Function to get information on reading end
end_treatment <- function(name, begin_info, signal_data) {
  end_time <- proc.time() # record it as soon as possible
  start_time <- begin_info[["start"]]
  delta_time <- end_time - start_time
  delta <- delta_time[]
  cat("End", name, "\n")
  cat("It lasted", round(delta["user.self"], 3), "s user time,",
      round(delta["sys.self"], 3), "s system time and",
      round(delta["elapsed"], 3), "s elapsed time.\n")
  if (begin_info[["force.real"]]) {
    # The imaginary part is left untouched
    i <- complex(real = 0, imaginary = 1)
    signal_data <- signal_data + i * Im(begin_info[["original_data"]])
  }
  if (begin_info[["vec"]]) {
    signal_data <- signal_data[1, ]
  }
  return(signal_data)
}

# Function to check arguments
check_arg <- function(arg, checks, can_be_null = FALSE) {
    check.list <- list(bool = c(is.logical, "a boolean"), int = c(function(x) { # nolint
    x %% 1 == 0
  }, "an integer"), num = c(is.numeric, "a numeric"),
  str = c(is.character, "a string"), pos = c(function(x) {
    x > 0
  }, "positive"), pos0 = c(function(x) {
    x >= 0
  }, "positive or zero"), l1 = c(function(x) {
    length(x) == 1
  }, "of length 1"))
  if (is.null(arg)) {
    if (!can_be_null) {
      stop(deparse(substitute(arg)), " is null.")
    }
  }
}


# Function to get arguments
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

# Function to get title and name samples
# Get the name of the signal from the title file or
# from the name of the subdirectory
get_title <- function(path, l, subdirs) {
  title <- NULL
  title_file <- file.path(file.path(file.path(path, "pdata"), "1"), "title")
  if (file.exists(title_file)) {
    lines <- readLines(title_file, warn = FALSE)
    if (length(lines) >= 1) {
      first_line <- gsub("^\\s+|\\s+$", "", lines[l])
      if (nchar(first_line) >= 1) {
        title <- first_line
      } else {
        warning(paste("The", l, "line of the title file is blank for 
                      directory ", path, "and the (sub)dirs names
                      are used instead"))
      }
    } else {
      warning(paste("Title file doesn't exists for directory ",
                    path, "\n the (sub)dirs names are  used instead"))
    }
  } else {
    warning(paste("Title file doesn't exists for directory ", path,
                  "\n the (sub)dirs names are  used instead"))
  }
  if (is.null(title)) {
    if (subdirs) {
      separator <- .Platform$file.sep
      path_elem <- strsplit(path, separator)[[1]]
      title <- paste(path_elem[length(path_elem) - 1],
                     path_elem[length(path_elem)], sep = "_")
    } else {
      title <- basename(path)
    }
  }
  return(title)
}

# Function to read parameter values for fid_info in the read_fids function
read_params <- function(file, params_name) {
  is_digit <- function(c) {
    return(suppressWarnings(!is.na(as.numeric(c))))
  }
  lines <- readLines(file)
  params <- sapply(params_name, function(x) NULL)

  for (paramName in params_name) {
    # Find line with the parameter I add a '$' '=' in the pattern so that for
    # example 'TD0' is not found where I look for 'TD', LOCSW + WBSW when I look
    # for 'SW'
    pattern <- paste("\\$", paramName, "=", sep = "")
    occurences <- grep(pattern, lines)
    if (length(occurences) == 0L) {
      stop(paste(file, "has no field", pattern))
    }
    lines <- readLines(file)
    params <- sapply(params_name, function(x) NULL)

    for (paramName in params_name) {
      # Find the line with the parameter I add a '$' '=' in the pattern
      # so that for example 'TD0' is not found where I look for 'TD'
      # and LOCSW and WBSW when I look
      # for 'SW'
      pattern <- paste("\\$", paramName, "=", sep = "")
      occurences <- grep(pattern, lines)
      if (length(occurences) == 0L) {
        stop(paste(file, "has no field", pattern))
      }
      if (length(occurences) > 1L) {
        warning(paste(file, "has more that one field", pattern,
                      " I take the first one"))
      }
      line <- lines[occurences[1]]

      # Cut beginning and end of the line '##$TD= 65536' -> '65536'
      igual <- as.numeric(regexpr("=", line))

      first <- igual
      while (first <= nchar(line) && !is_digit(substr(line, first, first))) {
        first <- first + 1
      }
      last <- nchar(line)
      while (last > 0 && !is_digit(substr(line, last, last))) {
        last <- last - 1
      }
      params[paramName] <- as.numeric(substr(line, first, last))
    }
    line <- lines[occurences[1]]

    # Cut beginning and end of the line '##$TD= 65536' -> '65536'
    igual <- as.numeric(regexpr("=", line))

    first <- igual
    while (first <= nchar(line) && !is_digit(substr(line, first, first))) {
      first <- first + 1
    }
    last <- nchar(line)
    while (last > 0 && !is_digit(substr(line, last, last))) {
      last <- last - 1
    }
    params[paramName] <- as.numeric(substr(line, first, last))
  }
  return(params)
}

# Function to read all fid's in the directory
read_fids <- function(path, l = 1, subdirs = FALSE, dirs_names = FALSE) {
  # Data initialisation and checks
  begin_info <- begin_treatment("read_fids")
  check_arg(path, c("str"))
  check_arg(l, c("pos"))
  if (file.exists(path) == FALSE) {
    stop(paste("Invalid path:", path))
  }

  # Extract the FIDs and their info
  if (subdirs == FALSE) {
    fid_dirs <- get_dirs_containing_fid(path)
    n <- length(fid_dirs)
    if (n == 0L) {
      stop(paste("No valid fid in", path))
    }
    if (dirs_names) {
      separator <- .Platform$file.sep
      path_elem <- strsplit(fid_dirs, separator)
      fid_names <- sapply(path_elem, function(x) x[[length(path_elem[[1]])]])
    } else {
      fid_names <- sapply(X = fid_dirs, FUN = get_title, l = l,
                          subdirs = subdirs, USE.NAMES = FALSE)
    }

    for (i in 1:n) {
      fid_list <- read_fid(fid_dirs[i])
      fid <- fid_list[["fid"]]
      info <- fid_list[["params"]]
      m <- length(fid)
      if (i == 1) {
        fid_data <- matrix(nrow = n, ncol = m, dimnames = list(fid_names,
                                                               info[["DT"]] *
                                                                 (0:(m - 1))))
        fid_info <- matrix(nrow = n, ncol = length(info),
                           dimnames = list(fid_names, names(info)))
      }
      fid_data[i, ] <- fid
      fid_info[i, ] <- unlist(info)
    }
  } else {
    maindirs <- dir(path, full.names = TRUE) # subdirectories
    fid_data <- numeric()
    fid_info <- numeric()

    fid_dirs <- c()
    for (j in maindirs) {
      fd <- get_dirs_containing_fid(j) # recoved FIDs from subdirectories
      n <- length(fd)
      if (n > 0L) {
        fid_dirs <- c(fid_dirs, fd)
      } else {
        warning(paste("No valid fid in", j))
      }
    }

    if (dirs_names == TRUE) {
      if (length(fid_dirs) != length(dir(path))) {
        # at least one subdir contains more than 1 FID
        separator <- .Platform$file.sep
        path_elem <- strsplit(fid_dirs, separator)
        fid_names <- sapply(path_elem, function(x) {
          paste(x[[length(path_elem[[1]]) - 1]],
                x[[length(path_elem[[1]])]], sep = "_")
        })
      } else {
        fid_names <- dir(path)
      }
    } else {
      fid_names <- sapply(X = fid_dirs, FUN = get_title,
                          l = l, subdirs = subdirs, USE.NAMES = FALSE)
    }

    for (i in seq_len(fid_names)) {
      fid_list <- read_fid(fid_dirs[i])
      fid <- fid_list[["fid"]]
      info <- fid_list[["params"]]
      m <- length(fid)
      if (i == 1) {
        fid_data <- matrix(nrow = length(fid_names), ncol = m,
                           dimnames = list(fid_names,
                                           info[["DT"]] * (0:(m - 1))))
        fid_info <- matrix(nrow = length(fid_names),
                           ncol = length(info),
                           dimnames = list(fid_names, names(info)))
      }
      print(paste("i=", i, "fid_data=", ncol(fid_data), "fid=", length(fid)))

      fid_data[i, ] <- fid
      fid_info[i, ] <- unlist(info)
    }
  }

  # Check for non-unique IDs
  nonnunique_ids <- sum(duplicated(row.names(fid_data)))
  cat("dim fid_data: ", dim(fid_data), "\n")
  cat("IDs: ", rownames(fid_data), "\n")
  cat("non-unique IDs?", nonnunique_ids, "\n")
  if (nonnunique_ids > 0) {
    warning("There are duplicated IDs: ", fid_data[duplicated(fid_data)])
  }

  # Return the results
  return(list(fid_data = end_treatment("read_fids",
                                       begin_info, fid_data),
              fid_info = fid_info))
}