#' Generates a count table with a single column per replicate
#' 
#' Creates a count table by collapsing multiple replicates into single columns.
#' This is useful, for example, when displaying line plots of expression
#' profiles across time.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts A count matrix
#' @param conditions  Vector of character condition names in the order of columns.
#'
#' @return Data frame with counts averaged within each condition. 
combine_replicates <- function(counts, conditions) {
    # Convert factor to characters if needed and get unique list
    conditions <- as.character(conditions)
    unique_conditions <- unique(conditions)

    # Create a new data.frame by combining replicates from first condition
    rep_counts <- counts[,conditions == unique_conditions[1]]

    if (is.matrix(rep_counts)) {
        # multiple replicates for condition
        combined <- rowMeans(rep_counts)
    } else {
        # single sample for condition
        combined <- rep_counts
    }

    # Add remaining conditions
    for (condition in unique_conditions[2:length(unique_conditions)]) {
        rep_counts <- counts[,conditions == condition]

        if (is.matrix(rep_counts)) {
            # multiple replicates for condition
            combined <- cbind(combined, rowMeans(rep_counts))
        } else {
            # single sample for condition
            combined <- cbind(combined, rep_counts)
        }
    }

    # Assign row and column names
    # Note that adding gene_id causes all of the other columns in the matrix to
    # become non-numeric. Therefore we will generate a separate matrix for output.
    #colnames(combined) <- paste0(rep('condition.', 3), conditions)
    colnames(combined) <- unique_conditions

    return(combined)
}

#' Counts-per-million (CPM) transformation
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param x Gene by sample count matrix
#'
#' @return Matrix of CPM-normalized counts.
counts_per_million <- function (x) {
    sweep(x, 2, colSums(x), '/') * 1E6
}

#' Computes the average gene expression profile for each co-expression module
#' detected.
#'
#' Given a count table with module assignments for each row, returns the average
#' expression profiles corresponding to each module.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts ExpressionSet Count data
#' @param module_assignments Vector of module assignments in the order of the
#' count rows
#' @param method Function to use to average expression profiles in each
#' module (e.g. mean or median)
#'
#' @return Matrix of module-average expression values across all conditions.
get_module_averages <- function(counts, module_assignments, method=mean) {
    # data frame with expression and module assignments
    module_expr <- as.data.frame(exprs(counts) * 1.0)
    module_expr$module <- module_assignments
    module_expr <- tbl_df(module_expr)

    # average expression for genes in each module
    module_averages <- as.matrix(module_expr %>% 
        group_by(module) %>% 
        select(-module) %>% 
        summarise_each(funs(method)) %>%
        select(-module))

    # modules appear in dataframe in sorted order
    rownames(module_averages) <- sort(unique(module_assignments))

    return(module_averages)
}

#' Get short condition names
#'
#' Takes a list of conditions, and a mapping from those condition names to
#' shorter versions, and returns a list of shortened condition names, possibly
#' as numeric.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param conditions long condition names
#' @param mapping data.frame mapping of long and short condition names
#' @param as_numeric Whether or not to convert the conditions to numeric, if
#' possible.
#'
#' @return Vector of numeric or character sample labels to be used in plots.
get_short_condition_names <- function(conditions, mapping, as_numeric=TRUE) {
    short <- mapping$short[match(conditions, mapping$long)]

    # if stage conditions are included (e.g. "proc"), stick to using strings
    if (isTRUE(as_numeric) && !any(is.na(suppressWarnings(as.numeric(short))))) {
        # drop unused conditions and return numeric condition labels
        short <- as.numeric(as.character(factor(short)))
    }
    return(short)
}

#' Converts a wide count table to long format
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts count table in wide format
#' @param condition_mapping data.frame mapping of long and short condition names
#' @param shorten_names Whether or not to shorten condition names
#'
#' @return Data frame of RNA-Seq counts in long format.
melt_counts <- function(counts, condition_mapping, shorten_names=TRUE) {
    # convert data into a more useful format for plotting
    counts_long <- melt(counts)

    # fix column names
    names(counts_long) <- c('id', 'condition', 'log_cpm')

    # If all conditions include a time component (e.g. amast04, amast24)
    # then we can convert the condition to numeric to improve plotting
    if (shorten_names) {
        counts_long$condition <- get_short_condition_names(counts_long$condition,
                                                          condition_mapping)
    }

    # If non-numeric conditions (e.g. 'procyclic') are included, convert to
    # factor and reorder levels of dataframe for better plotting
    if (!is.numeric(counts_long$condition)) {
        counts_long$condition <- factor(counts_long$condition, 
                                       unique(condition_mapping$short))
    }
    return(counts_long)
}

