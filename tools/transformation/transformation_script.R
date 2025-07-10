## Etienne Thevenot
## W4M Core Development Team
## etienne.thevenot@cea.fr
## 2015-04-25

transformF <- function(datMN,
                       metC) {
    
    ## options
    
    optStrAsFacL <- options()[["stringsAsFactors"]]
    options(stringsAsFactors = FALSE)
    
    ## checking		     
    
    if(length(which(datMN < 0))) {
        cat("\nThe 'dataMatrix' contains negative values\n")
        sink()
        stop("The 'dataMatrix' contains negative values", call. = FALSE)
    }
    
    ## Number of missing values
    nasN <- length(which(is.na(datMN)))
    cat("\nMissing values in the 'dataMatrix': ",
        nasN,
        " (",
        round(nasN / cumprod(dim(datMN))[2] * 100),
        "%)\n",
        sep="")
    
    ## Number of zero values
    zerN <- length(which(datMN == 0))
    cat("\nZero values in the 'dataMatrix': ",
        zerN,
        " (",
        round(zerN / cumprod(dim(datMN))[2] * 100),
        "%)\n",
        sep="")
    
    ## transformation
    
    switch(metC,
           log2 = {

               cat("\n'log2' transformation\n", sep="")

               trfMN <- log2(1 + datMN)

           },
           log10 = {

               cat("\n'log10' transformation\n", sep="")

               trfMN <- log10(1 + datMN)

           },
           sqrt = {

               cat("\n'Square root' transformation\n", sep="")
               
               trfMN <- sqrt(datMN)


           }) ## end of method


    ## returning

    options(stringsAsFactors=optStrAsFacL)

    return(trfMN)

} ## end of transformF
