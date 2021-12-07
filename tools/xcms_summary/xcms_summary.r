#!/usr/bin/env Rscript



# ----- ARGUMENTS BLACKLIST -----
#xcms.r
argBlacklist <- c("zipfile", "singlefile_galaxyPath", "singlefile_sampleName", "xfunction", "xsetRdataOutput", "sampleMetadataOutput", "ticspdf", "bicspdf", "rplotspdf")
#CAMERA.r
argBlacklist <- c(argBlacklist, "dataMatrixOutput", "variableMetadataOutput", "new_file_path")


# ----- PACKAGE -----
cat("\tSESSION INFO\n")

#Import the different functions
source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE); base_dir <- dirname(substring(argv[grep("--file=", argv)], 8)); source(paste(base_dir, fname, sep = "/"))
}
source_local("lib.r")

pkgs <- c("CAMERA", "batch")
loadAndDisplayPackages(pkgs)
cat("\n\n");


# ----- FUNCTION -----
writehtml <- function(...) {
  cat(..., "\n", file = htmlOutput, append = TRUE, sep = "")
}
writeraw <- function(htmlOutput, object, open = "at") {
    log_file <- file(htmlOutput, open = open)
    sink(log_file)
    sink(log_file, type = "output")
        print(object)
    sink()
    close(log_file)
}
getSampleNames <- function(xobject) {
    if (class(xobject) == "xcmsSet")
        return(sampnames(xobject))
    if (class(xobject) == "XCMSnExp")
        return(xobject@phenoData@data$sample_name)
}
getFilePaths <- function(xobject) {
    if (class(xobject) == "xcmsSet")
        return(xobject@filepaths)
    if (class(xobject) == "XCMSnExp")
        return(fileNames(xobject))
}
equalParams <- function(param1, param2) {
    writeraw("param1.txt", param1, open = "wt")
    writeraw("param2.txt", param2, open = "wt")
    return(tools::md5sum("param1.txt") == tools::md5sum("param2.txt"))
}


# ----- ARGUMENTS -----

args <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects


# ----- ARGUMENTS PROCESSING -----

#image is an .RData file necessary to use xset variable given by previous tools
load(args$image);

htmlOutput <- "summary.html"
if (!is.null(args$htmlOutput)) htmlOutput <- args$htmlOutput;

user_email <- NULL
if (!is.null(args$user_email)) user_email <- args$user_email;

# if the RData come from XCMS 1.x
if (exists("xset")) {
    xobject <- xset
    # retrocompatability
    if (!exists("sampleNamesList")) sampleNamesList <- list("sampleNamesMakeNames" = make.names(sampnames(xobject)))
}
# if the RData come from CAMERA
if (exists("xa")) {
    xobject <- xa@xcmsSet
    if (!exists("sampleNamesList")) sampleNamesList <- list("sampleNamesMakeNames" = make.names(xa@xcmsSet@phenoData$sample_name))
}
# if the RData come from XCMS 3.x
if (exists("xdata")) {
    xobject <- xdata
    if (!exists("sampleNamesList")) sampleNamesList <- list("sampleNamesMakeNames" = make.names(xdata@phenoData@data$sample_name))
}

if (!exists("xobject")) stop("You need at least a xdata, a xset or a xa object.")


# ----- MAIN PROCESSING INFO -----
writehtml("<!DOCTYPE html>")
writehtml("<HTML lang='en'>")

writehtml("<HEAD>")
    writehtml("<meta http-equiv='Content-Type' content='text/html; charset=UTF-8' />")

    writehtml("<title>[W4M] XCMS analysis summary</title>")

    writehtml("<style>")
        writehtml("table, tr, td, th { border: 1px solid #000000; border-collapse:collapse; }")
        writehtml("td, th { padding: 5px; padding-right: 12px; }")
        writehtml("th { background: #898989; text-align:left;color: white;}")
        writehtml("h2 { color: #FFA212; }")
        writehtml("ul li { margin-bottom:10px; }")
    writehtml("</style>")
writehtml("</HEAD>")

writehtml("<BODY>")
    writehtml("<div><h1>___ XCMS analysis summary using Workflow4Metabolomics ___</h1>")
    # to pass the planemo shed_test
    if (user_email != "test@bx.psu.edu") {
        if (!is.null(user_email)) writehtml("By: ", user_email, " - ")
        writehtml("Date: ", format(Sys.time(), "%y%m%d-%H:%M:%S"))
    }
    writehtml("</div>")

    writehtml("<h2>Samples used:</h2>")
    writehtml("<div><table>")
        if (all(getSampleNames(xobject) == sampleNamesList$sampleNamesMakeNames)) {
            sampleNameHeaderHtml <- paste0("<th>sample</th>")
            sampleNameHtml <- paste0("<td>", getSampleNames(xobject), "</td>")
        } else {
            sampleNameHeaderHtml <- paste0("<th>sample</th><th>sample renamed</th>")
            sampleNameHtml <- paste0("<td>", getSampleNames(xobject), "</td><td>", sampleNamesList$sampleNamesMakeNames, "</td>")
        }

        if (!exists("md5sumList")) {
            md5sumHeaderHtml <- ""
            md5sumHtml <- ""
            md5sumLegend <- ""
        } else if (is.null(md5sumList$removalBadCharacters)) {
            md5sumHeaderHtml <- paste0("<th>md5sum<sup>*</sup></th>")
            md5sumHtml <- paste0("<td>", md5sumList$origin, "</td>")
            md5sumLegend <- "<br/><sup>*</sup>The program md5sum is designed to verify data integrity. So you can check if the data were uploaded correctly or if the data were changed during the process."
        } else {
            md5sumHeaderHtml <- paste0("<th>md5sum<sup>*</sup></th><th>md5sum<sup>**</sup> after bad characters removal</th>")
            md5sumHtml <- paste0("<td>", md5sumList$origin, "</td><td>", md5sumList$removalBadCharacters, "</td>")
            md5sumLegend <- "<br/><sup>*</sup>The program md5sum is designed to verify data integrity. So you can check if the data were uploaded correctly or if the data were changed during the process.<br/><sup>**</sup>Because some bad characters (eg: accent) were removed from your original file, the checksum have changed too.<br/>"
        }

        writehtml("<tr>", sampleNameHeaderHtml, "<th>filename</th>", md5sumHeaderHtml, "</tr>")
        writehtml(paste0("<tr>", sampleNameHtml, "<td>", getFilePaths(xobject), "</td>", md5sumHtml, "</tr>"))

    writehtml("</table>")
    writehtml(md5sumLegend)
    writehtml("</div>")

    writehtml("<h2>Function launched:</h2>")
    writehtml("<div><table>")
        writehtml("<tr><th>timestamp<sup>***</sup></th><th>function</th><th>argument</th><th>value</th></tr>")
        # XCMS 3.x
        if (class(xobject) == "XCMSnExp") {
            xcmsFunction <- NULL
            params <- NULL
            for (processHistoryItem in processHistory(xobject)) {
                if ((xcmsFunction == processType(processHistoryItem)) && equalParams(params, processParam(processHistoryItem)))
                    next
                timestamp <- processDate(processHistoryItem)
                xcmsFunction <- processType(processHistoryItem)
                params <- processParam(processHistoryItem)
                writehtml("<tr><td>", timestamp, "</td><td>", xcmsFunction, "</td><td colspan='2'><pre>")
                writeraw(htmlOutput, params)
                writehtml("</pre></td></tr>")
            }
        }
        # CAMERA and retrocompatability XCMS 1.x
        if (exists("listOFlistArguments")) {
            for (tool in names(listOFlistArguments)) {
                listOFlistArgumentsDisplay <- listOFlistArguments[[tool]][!(names(listOFlistArguments[[tool]]) %in% argBlacklist)]

                timestamp <- strsplit(tool, "_")[[1]][1]
                xcmsFunction <- strsplit(tool, "_")[[1]][2]
                writehtml("<tr><td rowspan='", length(listOFlistArgumentsDisplay), "'>", timestamp, "</td><td rowspan='", length(listOFlistArgumentsDisplay), "'>", xcmsFunction, "</td>")
                line_begin <- ""
                for (arg in names(listOFlistArgumentsDisplay)) {
                    writehtml(line_begin, "<td>", arg, "</td><td>", unlist(listOFlistArgumentsDisplay[arg][1]), "</td></tr>")
                    line_begin <- "<tr>"
                }
            }
        }
    writehtml("</table>")
    writehtml("<br/><sup>***</sup>timestamp format: DD MM dd hh:mm:ss YYYY or yymmdd-hh:mm:ss")
    writehtml("</div>")

    if (class(xobject) == "XCMSnExp") {
        writehtml("<h2>Informations about the XCMSnExp object:</h2>")

        writehtml("<div><pre>")
            writeraw(htmlOutput, xobject)
        writehtml("</pre></div>")
    }

    writehtml("<h2>Informations about the xcmsSet object:</h2>")

    writehtml("<div><pre>")
        # Get the legacy xcmsSet object
        xset <- getxcmsSetObject(xobject)
        writeraw(htmlOutput, xset)
    writehtml("</pre></div>")

    # CAMERA
    if (exists("xa")) {
        writehtml("<h2>Informations about the CAMERA object:</h2>")

        writehtml("<div>")
            writehtml("Number of pcgroup: ", length(xa@pspectra))
        writehtml("</div>")
    }

    writehtml("<h2>Citations:</h2>")
    writehtml("<div><ul>")
        writehtml("<li>To cite the <b>XCMS</b> package in publications use:")
            writehtml("<ul>")
            writehtml("<li>", "Smith, C.A. and Want, E.J. and O'Maille, G. and Abagyan,R. and Siuzdak, G.XCMS: Processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching and identification, Analytical Chemistry, 78:779-787 (2006)", "</li>")
            writehtml("<li>", "Ralf Tautenhahn, Christoph Boettcher, Steffen Neumann: Highly sensitive feature detection for high resolution LC/MS BMC Bioinformatics, 9:504 (2008)", "</li>")
            writehtml("<li>", "H. Paul Benton, Elizabeth J. Want and Timothy M. D. Ebbels Correction of mass calibration gaps in liquid chromatography-mass spectrometry metabolomics data Bioinformatics, 26:2488 (2010)", "</li>")
            writehtml("</ul>")
        writehtml("</li>")

        writehtml("<li>To cite the <b>CAMERA</b> package in publications use:")
            writehtml("<ul>")
            writehtml("<li>", "Kuhl, C., Tautenhahn, R., Boettcher, C., Larson, T. R. and Neumann,S. CAMERA: an integrated strategy for compound spectra extraction and annotation of liquid chromatography/mass spectrometry data sets. Analytical Chemistry, 84:283-289 (2012)", "</li>")
            writehtml("</ul>")
        writehtml("</li>")

        writehtml("<li>To cite the <b>Workflow4Metabolimics (W4M)</b> project in publications use:")
            writehtml("<ul>")
            writehtml("<li>", "Franck Giacomoni, Gildas Le Corguillé, Misharl Monsoor, Marion Landi, Pierre Pericard, Mélanie Pétéra, Christophe Duperier, Marie Tremblay-Franco, Jean-François Martin, Daniel Jacob, Sophie Goulitquer, Etienne A. Thévenot and Christophe Caron (2014). Workflow4Metabolomics: A collaborative research infrastructure for computational metabolomics. Bioinformatics  doi:10.1093/bioinformatics/btu813", "</li>")
            writehtml("</ul>")
        writehtml("</li>")
    writehtml("</ul></div>")

writehtml("</BODY>")

writehtml("</HTML>")
