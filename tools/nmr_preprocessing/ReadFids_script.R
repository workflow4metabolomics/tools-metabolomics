################################################################################################
#
#   Read FIDs in Bruker format
#
#
################################################################################################


# vec2mat ==============================================================================
vec2mat <- function(vec) {
  return(matrix(vec, nrow = 1, dimnames = list(c(1), names(vec)))) 
}


# ReadFid ==============================================================================
ReadFid <- function(path) {
  
  # Read 1D FID using Bruker XWinNMR and TopSpin format.  It is inspired of the
  # matNMR matlab library which deals with 2D FID and also other formats
  # Read also the parameters in the acqus file
  
  paramFile <- file.path(path, "acqus")
  # BYTEORDA: 0 -> Little Endian 1 -> Big Endian
  params <- readParams(paramFile, c("TD", "BYTORDA", "DIGMOD", "DECIM", "DSPFVS", 
                                    "SW_h", "SW", "O1"))
  console <- readParamsChar(paramFile, c("INSTRUM"))
  neo <- grepl("neo", tolower(console))
  
  if (params[["DSPFVS"]] >= 20) {
    # The group delay first order phase correction is given directly from version 20
    grpdly <- readParams(paramFile, c("GRPDLY"))
    params[["GRPDLY"]] <- grpdly[["GRPDLY"]]
  }
  TD <- params[["TD"]]

  endianness <- if (params$BYTORDA) 
    "big" else "little"
  if (TD%%2 != 0) {
    stop(paste("Only even numbers are allowed for size in TD because it is complex 
               data with the real and imaginary part for each element.", 
               "The TD value is in the", paramFile, "file"))
  }
  
  # Interpret params Dwell Time, time between 2 data points in the FID
  params[["DT"]] <- 1/(2 * params[["SW_h"]])
  
  # Read fid
  fidFile <- file.path(path, "fid")
  if(neo)
    fidOnDisk <- readBin(fidFile, what = "double", n = TD, size = NA_integer_, signed = TRUE, endian = .Platform$endian)
  else
    fidOnDisk <- readBin(fidFile, what = "int", n = TD, size = 4L, endian = endianness)
  
  # Real size that is on disk (it should be equal to TD2, except for TopSpin/Bruker
  # (which is our case) according to matNMR as just discussed
  TDOnDisk <- length(fidOnDisk)
  if (TDOnDisk < TD) {
    warning("Size is smaller than expected, the rest is filled with zero so the size is the same for every fid")
    fidGoodSize <- sapply(vector("list", length = TD), function(x) 0)
    fidGoodSize[1:TDOnDisk] <- fidOnDisk
    
  } else if (TDOnDisk > TD) {
    warning("Size is bigger than expected, the rest ignored so the size is the same for every fid")
    fidGoodSize <- fidOnDisk(1:TD)
    
  } else {
    fidGoodSize <- fidOnDisk
  }
  
  fidRePart <- fidGoodSize[seq(from = 1, to = TD, by = 2)]
  fidImPart <- fidGoodSize[seq(from = 2, to = TD, by = 2)]
  fid <- complex(real = fidRePart, imaginary = fidImPart)
  
  return(list(fid = fid, params = params))
}




# getDirsContainingFid ==============================================================================
getDirsContainingFid <- function(path) {
  subdirs <- dir(path, full.names = TRUE)
  if (length(subdirs) > 0) {
    cond <- sapply(subdirs, function(x)  {
      content <- dir(x)
      # subdirs must contain fid, acqu and acqus files
      return("fid" %in% content && "acqu" %in% content && "acqus" %in% content)
    })
    subdirs <- subdirs[cond]
  }
  return(subdirs)
}






# beginTreatment ==============================================================================

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
  
  return(list(start = proc.time(), vec = vec, force.real = force.real, 
              Original_data = Original_data, Signal_data = Signal_data, Signal_info = Signal_info))
}


# endTreatment ==============================================================================

endTreatment <- function(name, begin_info, Signal_data) {
  end_time = proc.time() # record it as soon as possible
  start_time = begin_info[["start"]]
  delta_time = end_time - start_time
  delta = delta_time[]
  cat("End", name, "\n")
  cat("It lasted",
      round(delta["user.self"], 3), "s user time,",
      round(delta["sys.self"] , 3), "s system time and",
      round(delta["elapsed"]  , 3), "s elapsed time.\n")
  if (begin_info[["force.real"]]) {
    # The imaginary part is left untouched
    i <- complex(real=0, imaginary=1)
    Signal_data = Signal_data + i * Im(begin_info[["Original_data"]])
  }
  if (begin_info[["vec"]]) {
    Signal_data = Signal_data[1,]
  }
  return(Signal_data)
}


# checkArg ==============================================================================

checkArg <- function(arg, checks, can.be.null=FALSE) {
  check.list <- list(bool=c(is.logical, "a boolean"),
                     int =c(function(x){x%%1==0}, "an integer"),
                     num =c(is.numeric, "a numeric"),
                     str =c(is.character, "a string"),
                     pos =c(function(x){x>0}, "positive"),
                     pos0=c(function(x){x>=0}, "positive or zero"),
                     l1 =c(function(x){length(x)==1}, "of length 1")
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


# getArg ==============================================================================

getArg <- function(arg, info, argname, can.be.absent=FALSE) {
  if (is.null(arg)) {
    start <- paste("impossible to get argument", argname, "it was not given directly and");
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
    arg <- info[1,argname]
    if (is.na(arg)) {
      stop(paste(start, "it is NA in the info matrix"))
    }
  }
  return(arg)
}



# getTitle ==============================================================================

# Get the name of the signal from the title file or fromt the name of the subdirectory
# Get the name of the signal from the title file or fromt the name of the subdirectory

getTitle <- function(path, l, subdirs) {
  title <- NULL
  title_file <- file.path(file.path(file.path(path, "pdata"), "1"), "title")
  if (file.exists(title_file)) {
    lines <- readLines(title_file, warn = FALSE)
    if (length(lines) >= 1)  {
      first_line <- gsub("^\\s+|\\s+$", "", lines[l])
      if (nchar(first_line) >= 1)  {
        title <- first_line
      } else {
        warning(paste("The", l ,"line of the title file is blank for directory ", 
                      path, "and the (sub)dirs names are used instead"))
      }
    } else {
      warning(paste("The title file is empty for directory ", path, "and the (sub)dirs names are used instead"))
    }
  } else {
    warning(paste("Title file doesn't exists for directory ", path, "\n the (sub)dirs names are  used instead"))
  }
  if (is.null(title)) {
    if(subdirs) {
      separator <- .Platform$file.sep
      path_elem <- strsplit(path,separator)[[1]]
      title <- paste(path_elem[length(path_elem)-1], path_elem[length(path_elem)], sep = "_")
    } else{title <- basename(path)} 
  }
  return(title)
}




# readParams ==============================================================================
# Read parameter values for Fid_info in the ReadFids function
readParams <- function(file, paramsName) {
  
  isDigit <- function(c) {
    return(suppressWarnings(!is.na(as.numeric(c))))
  }
  lines <- readLines(file)
  params <- sapply(paramsName, function(x) NULL)
  
  for (paramName in paramsName)  {
    # Find the line with the parameter I add a '$' '=' in the pattern so that for
    # example 'TD0' is not found where I look for 'TD' and LOCSW and WBSW when I look
    # for 'SW'
    pattern <- paste("\\$", paramName, "=", sep = "")
    occurences <- grep(pattern, lines)
    if (length(occurences) == 0L)  {
      stop(paste(file, "has no field", pattern))
    }
    if (length(occurences) > 1L) {
      warning(paste(file, "has more that one field", pattern, " I take the first one"))
    }
    line <- lines[occurences[1]]
    
    # Cut beginning and end of the line '##$TD= 65536' -> '65536'
    igual <- as.numeric(regexpr("=", line))
    
    first <- igual
    while (first <= nchar(line) & !isDigit(substr(line, first, first))) {
      first <- first + 1
    }
    last <- nchar(line)
    while (last > 0 & !isDigit(substr(line, last, last)))  {
      last <- last - 1
    }
    params[paramName] <- as.numeric(substr(line, first, last))
  }
  return(params)
}

readParamsChar <- function(file, paramsName) {
  lines <- readLines(file)
  params <- sapply(paramsName, function(x) NULL)
  
  for (paramName in paramsName)  {
    # Find the line with the parameter I add a '$' '=' in the pattern so that for
    # example 'TD0' is not found where I look for 'TD' and LOCSW and WBSW when I look
    # for 'SW'
    pattern <- paste("\\$", paramName, "=", sep = "")
    occurences <- grep(pattern, lines)
    if (length(occurences) == 0L)  {
      stop(paste(file, "has no field", pattern))
    }
    if (length(occurences) > 1L) {
      warning(paste(file, "has more that one field", pattern, " I take the first one"))
    }
    line <- lines[occurences[1]]
    
    # Cut beginning and end of the line '##$TD= 65536' -> '65536'
    igual <- as.numeric(regexpr("=", line))
    
    first <- igual
    while (first <= nchar(line) & substr(line, first, first) != "<") {
      first <- first + 1
    }
    first <- first + 1
    last <- nchar(line)
    while (last > 0 & substr(line, last, last) != ">")  {
      last <- last - 1
    }
    last <- last - 1
    params[paramName] <- substr(line, first, last)
  }
  return(params)
}

# ReadFids ==============================================================================

ReadFids <- function(path, l = 1, subdirs = FALSE, dirs.names = FALSE) {
  
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("ReadFids")
  checkArg(path, c("str"))
  checkArg(l, c("pos"))
  if (file.exists(path) == FALSE) {
    stop(paste("Invalid path:", path))
  }
  
  
  # Extract the FIDs and their info ----------------------------------------------
  
  if (subdirs == FALSE) {
    fidDirs <- getDirsContainingFid(path)
    n <- length(fidDirs)
    if (n == 0L)  {
      stop(paste("No valid fid in", path))
    }
    if (dirs.names) {
      separator <- .Platform$file.sep
      path_elem <- strsplit(fidDirs,separator)
      fidNames <- sapply(path_elem, function(x) x[[length(path_elem[[1]])]])
    }else {fidNames <- sapply(X = fidDirs, FUN = getTitle, l = l, subdirs = subdirs,  USE.NAMES = F)}
    
    for (i in 1:n)  {
      fidList <- ReadFid(fidDirs[i])
      fid <- fidList[["fid"]]
      info <- fidList[["params"]]
      m <- length(fid)
      if (i == 1)  {
        Fid_data <- matrix(nrow = n, ncol = m, dimnames = list(fidNames, 
                                                               info[["DT"]] * (0:(m - 1))))
        Fid_info <- matrix(nrow = n, ncol = length(info), dimnames = list(fidNames, 
                                                                          names(info)))
      }
      Fid_data[i, ] <- fid
      Fid_info[i, ] <- unlist(info)
    }
    
  } else  {
    maindirs <- dir(path, full.names = TRUE) # subdirectories
    Fid_data <- numeric()
    Fid_info <- numeric()
    
    fidDirs <- c()
    for (j in maindirs) {
      fd <- getDirsContainingFid(j) # recoved FIDs from subdirectories
      n <- length(fd)
      if (n > 0L)  {
        fidDirs <- c(fidDirs, fd)
      } else {warning(paste("No valid fid in",j ))}
    }
    
    if (dirs.names==TRUE) {
      if (length(fidDirs)!= length(dir(path))) { # at least one subdir contains more than 1 FID
        separator <- .Platform$file.sep
        path_elem <- strsplit(fidDirs,separator)
        fidNames <- sapply(path_elem, function(x) paste(x[[length(path_elem[[1]])-1]],
                                                        x[[length(path_elem[[1]])]], sep = "_"))
      }else {fidNames <- dir(path)}
      
    } else {fidNames <- sapply(X = fidDirs, FUN = getTitle, l = l, subdirs = subdirs, USE.NAMES = F)}
    
    for (i in 1:length(fidNames))  {
      fidList <- ReadFid(fidDirs[i])
      fid <- fidList[["fid"]]
      info <- fidList[["params"]]
      m <- length(fid)
      if (i == 1)  {
        Fid_data <- matrix(nrow = length(fidNames), ncol = m, dimnames = list(fidNames, 
                                                                              info[["DT"]] * (0:(m - 1))))
        Fid_info <- matrix(nrow = length(fidNames), ncol = length(info), dimnames = list(fidNames, 
                                                                                         names(info)))
      }
      print(paste("i=", i, "Fid_data=", ncol(Fid_data), "fid=", length(fid)))

      Fid_data[i, ] <- fid
      Fid_info[i, ] <- unlist(info)
    }
  }
  
  # Check for non-unique IDs ----------------------------------------------
  NonnuniqueIds <- sum(duplicated(row.names(Fid_data)))
  cat("dim Fid_data: ", dim(Fid_data), "\n")
  cat("IDs: ", rownames(Fid_data), "\n")
  cat("non-unique IDs?", NonnuniqueIds, "\n")
  if (NonnuniqueIds > 0) {
    warning("There are duplicated IDs: ", Fid_data[duplicated(Fid_data)])
  }
  
  
  # Return the results ----------------------------------------------
  return(list(Fid_data = endTreatment("ReadFids", begin_info, Fid_data), Fid_info = Fid_info))
}

