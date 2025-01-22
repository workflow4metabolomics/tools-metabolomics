#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file
# For questions: Antoine Gravot (Protocole conception) and Misharl Monsoor  (for galaxy wrapper and R script).

# Load the different libraries
library(batch) # necessary for parseCommandArgs function
library(reshape) # necessary for using melt function
library(MASS) # necessary for using the write.matrix()
# interpretation of arguments given in command line as an R list of objects
list_arguments <- parseCommandArgs(evaluate = FALSE)

cat(
    "\nJob starting time:\n", format(Sys.time(), "%a %d %b %Y %X"),
    "\n\n--------------------------------------------------------------------",
    "\nParameters used in 'Metabolites Correlation Analysis':\n\n"
)
print(list_arguments)
cat("--------------------------------------------------------------------\n\n")

# The main function of this script that will execute all the other functions


main_function <- function(sorting, variable_metadata, data_matrix, sample_metadata, corrdel, param_correlation, param_cytoscape, matrix_corr, user_matrix_corr, corr_method) {
    if (sorting == 1) {
        #### Executing the sorting function####
        cat("\nExecuting the sorting function\n")
        # Read the tsv annotateDiffreport file and don't modify the columns name (check.names = FALSE)
        variable_metadata_input <- read.csv(variable_metadata, header = T, sep = "\t", dec = ".", check.names = FALSE)
        data_matrix_input <- read.csv(data_matrix, header = T, sep = "\t", dec = ".", check.names = FALSE)
        first_column_variable <- toString(names(variable_metadata_input)[1])
        first_column_datamatrix <- toString(names(data_matrix_input)[1])
        # @TODO merge
        input_tsv <- cbind(variable_metadata_input, data_matrix_input[, !(colnames(data_matrix_input) %in% c(first_column_datamatrix))])
        # Load the sample.info from the xcmsSet
        sample_metadata_info_tsv <- read.table(sample_metadata, header = T, sep = "\t", dec = ",", check.names = FALSE)
        # Extract the samples name from the sample.info file generated from the xcmsSet step in ABIMS Workflow4Metabo.
        samples_name <- as.vector(t(sample_metadata_info_tsv[1]))
        output_tsv <- sorting(input_tsv, samples_name)
        # output now if corrdel == 0
        if (corrdel == 0) {
            output_tsv_vm <- output_tsv[, which(!(colnames(output_tsv) %in% samples_name))]
            write.table(output_tsv_vm[, c(3:ncol(output_tsv_vm), 2, 1)], sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, file = "sorted_table.tsv")
        }
    }
    if (corrdel == 1) {
        cat("\nExecuting the corr_matrix_del function\n")
        corr_matrix_del(output_tsv, samples_name, param_correlation, param_cytoscape, corr_method)
    }

    if (matrix_corr == 1) {
        cat("\nExecuting the corr_matrix function\n")
        corr_matrix_user(user_matrix_corr, param_cytoscape, corr_method)
    }
}

# The sorting function will sort the dataframe by rt column.
# Then it creates a "signal_moy" column which contains the mean values of the signal values of the sample by row.
# It finally creates a table tsv format "sorted_table.tsv".

sorting <- function(input_tsv, samples_name) {
    # Sort by rt column
    new_input <- input_tsv[order(input_tsv$rt), ]
    # Compute the mean operation of all the signal values of the sample by row, and put the results in a new column "signal_moy"
    new_input["signal_moy"] <- data.frame(Means = rowMeans(new_input[, colnames(new_input) %in% samples_name]))
    # Rearrange the data frame in order to have the columns "signal_moy" and "pcgroup" at the beginning of the table
    new_input <- cbind(new_input["signal_moy"], new_input["pcgroup"], new_input[, !(colnames(new_input) %in% c("signal_moy", "pcgroup"))])
    # Sort the "signal_moy" column data frame by pcgroup
    new_input <- new_input[order(new_input$pcgroup, -new_input$signal_moy), ]
    # Write the data frame to a file named "sorted_table.tsv"
    # Erase the rownames of the table
    rownames(new_input) <- NULL
    # Suppress rows which contains Nas
    new_input <- new_input[complete.cases(new_input$pcgroup), ]
    return(new_input)
}

corr_matrix_del <- function(output_tsv, samples_name, param_correlation, param_cytoscape, corr_method) {
    # statmatrix, a dataframe which contains only the columns "name" and the values of the intensity for all the samples
    first_column <- toString(names(output_tsv)[3])
    statmatrix <- output_tsv[(names(output_tsv) %in% c(first_column, samples_name))]
    statmatrix2 <- output_tsv[(names(output_tsv) %in% c(first_column, "signal_moy", samples_name))]
    n <- statmatrix[, first_column]
    # transpose all but the first column (name)
    statmatrix <- as.data.frame(t(statmatrix[, -1]))
    # Rename the columns of the dataframe "statmatrix"
    colnames(statmatrix) <- n
    # Transform dataframe to matrix before doing the cor function
    corr_transpo <- data.matrix(statmatrix)
    # Create a matrix which contains only the data of the samples without the other columns
    statmatrix <- output_tsv[(names(output_tsv) %in% c(samples_name))]
    # Add the rownames previously
    rownames(statmatrix) <- make.unique(as.vector(output_tsv[, first_column]))
    # Do the cor step, with the transposed statmatrix (metabolites as columns)
    corr_transpo <- cor(t(as.matrix(statmatrix)), method = corr_method)
    # Add the columns "signal_moy", "pcgroup" and "name" to the final correlation matrix output
    new_dataframe <- cbind(output_tsv["signal_moy"], output_tsv["pcgroup"], output_tsv[names(output_tsv)[3]], corr_transpo)
    # Add a new column "suppress" which indicates if the metabolite will be removed from the analysis (because it is correlated with the other metabolite)
    new_dataframe["suppress"] <- ""
    # Take the pcgroup names
    pcgroups <- unique(new_dataframe$pcgroup)
    # Remove NAs from the pcgroups vector
    pcgroups <- pcgroups[which(pcgroups != "NA")]
    # Creates a pcgroups list
    list_data <- vector(mode = "list", length = length(pcgroups))
    # Will do the following steps by pcgroup
    for (pcgroup in pcgroups) {
        # Select a subset data frame for each pcgroup
        subset_data <- new_dataframe[new_dataframe$pcgroup == pcgroup, ]
        # Stores the metabolites names for each pcgroup into a matrix one dimension.
        pc_group_elements <- t(as.matrix(subset_data[3]))
        # uses the matrix "pc_group_elements" containing the metabolites to have the correlation data frame "subset_data2" by pcgroup
        subset_data2 <- subset_data[names(subset_data) %in% c(names(output_tsv)[3], pc_group_elements)]
        # We keep the first metabolites for each pcgroup that correspond to the metabolites with the highest amount of signal
        first_metabolite <- pc_group_elements[1]
        # Remove the first metabolite from the pc_group_elements matrix
        pc_group_elements <- pc_group_elements[, c(seq_len(length(pc_group_elements)))]
        # print (pc_group_elements)
        for (metabolite in pc_group_elements) {
            metabolite_dataframe <- subset_data2[subset_data2[, first_column] == metabolite, ]
            # retrives metabolite index from the vector object "pc_group_elements"
            index_metabolite <- as.numeric(which(pc_group_elements == metabolite, arr.ind = TRUE))
            for (f in 2:index_metabolite - 1) {
                if (metabolite != first_metabolite) {
                    value_intensity <- t(data.matrix(metabolite_dataframe[, pc_group_elements[f]]))
                    value_intensity <- value_intensity[1:1]
                    # If one value is >0.75, it means that the metabolite est correlated with at least one previous metabolite in the table, so the boolean variable is set to true.
                    if (value_intensity > param_correlation) {
                        new_dataframe$suppress[new_dataframe[, first_column] == metabolite] <- "DEL"
                    }
                }
            }
        }
    }
    # The dataframe with the column "suppress" at the begginning
    new_dataframe_suppression <- cbind(new_dataframe[first_column], new_dataframe["suppress"], output_tsv["rt"], new_dataframe["signal_moy"], new_dataframe[, !(colnames(new_dataframe) %in% c(first_column, "suppress", "pcgroup", "signal_moy"))])


    # remove all the rows which have the "suppress" data and keep only the selected metabolites
    selected_metabolite_dataframe <- new_dataframe_suppression[new_dataframe_suppression$suppress != "DEL", ]
    # Keep the metabolites selected in a list
    metabolites_selected_list <- as.vector(t(selected_metabolite_dataframe[1]))
    # Add the other columns to keep
    metabolites_selected_list <- c(metabolites_selected_list, first_column, "rt", "signal_moy")
    # Keep the metabolites selected to keep only the columns
    selected_metabolite_dataframe <- selected_metabolite_dataframe[, colnames(selected_metabolite_dataframe) %in% metabolites_selected_list]

    # Export to siff table format for visualization in cytoscape for the selected metabolites
    siff_table <- melt(selected_metabolite_dataframe[, !colnames(selected_metabolite_dataframe) %in% c("pcgroup", "rt", "signal_moy")])
    # Remove the values equal to 1 (correlation between two metabolite identical)
    siff_table <- siff_table[siff_table$value != 1, ]
    # Keep only the values corresponding to the param_cytoscape
    siff_table <- siff_table[siff_table$value >= param_cytoscape, ]
    # Change the order of the columns
    siff_table <- cbind(siff_table[first_column], siff_table["value"], siff_table["variable"])



    # Join the two datasets to keep only the selected metabolite, with all the information about the intensity value for statistics analysis in the next step of the workflow.
    joined_dataframe <- merge(statmatrix2, selected_metabolite_dataframe[first_column], by = first_column)
    # Order by the signal intensity
    joined_dataframe <- joined_dataframe[order(-joined_dataframe$signal_moy), ]
    # Transposition of the dataframe
    joined_dataframe <- joined_dataframe[, !(colnames(joined_dataframe) %in% c("signal_moy"))]


    # Write the different tables into files
    write.table(selected_metabolite_dataframe, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, file = "correlation_matrix_selected.tsv")
    write.table(siff_table, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = "siff_table.tsv")
    output_tsv_vm <- output_tsv[, which(!(colnames(output_tsv) %in% samples_name))]
    output_tsv_vm <- data.frame(output_tsv_vm[, c(3:ncol(output_tsv_vm), 2, 1)], ori.or = seq_len(nrow(output_tsv_vm)))
    output_tsv_vm <- merge(x = output_tsv_vm, y = new_dataframe[, c(3, which(colnames(new_dataframe) == "suppress"))], by.x = 1, by.y = 1, sort = FALSE)
    output_tsv_vm <- output_tsv_vm[order(output_tsv_vm$ori.or), ][, -c(which(colnames(output_tsv_vm) == "ori.or"))]
    write.table(output_tsv_vm, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, file = "sorted_table.tsv")
}


corr_matrix_user <- function(user_matrix_corr, param_cytoscape, corr_method) {
    # read the input table
    input_tsv <- read.csv(user_matrix_corr, header = T, sep = "\t", dec = ".", check.names = FALSE)
    n <- input_tsv$HD
    # transpose all but the first column (name)
    statmatrix_corr <- as.data.frame(t(input_tsv[, -1]))
    # Rename the columns of the dataframe "statmatrix"
    colnames(statmatrix_corr) <- n
    # Do the cor step, with the transposed statmatrix.
    corr_transpo <- cor(t(as.matrix(statmatrix_corr)), method = corr_method)
    # Write the matrix to a tsv file
    write.table(corr_transpo, sep <- "\t", quote = FALSE, col.names = NA, file = "correlation_matrix.tsv")
    # Export to siff table format for visualization in cytoscape for the selected conditions
    siff_table <- melt(corr_transpo)
    # Remove the values equal to 1 (correlation between two metabolite identical)
    siff_table <- siff_table[siff_table$value != 1, ]
    # Keep only the values corresponding to the param_cytoscape
    siff_table <- siff_table[siff_table$value >= param_cytoscape, ]
    # Change the order of the columns
    siff_table <- cbind(Node1 = siff_table["X1"], interaction = siff_table["value"], Node2 = siff_table["X2"])
    # Write the siff table
    write.table(siff_table, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = "siff_table.tsv")
}
do.call(main_function, list_arguments)


cat(
    "\n--------------------------------------------------------------------",
    "\nInformation about R (version, Operating System, attached or loaded packages):\n\n"
)
sessionInfo()
cat(
    "--------------------------------------------------------------------\n",
    "\nJob ending time:\n", format(Sys.time(), "%a %d %b %Y %X")
)
