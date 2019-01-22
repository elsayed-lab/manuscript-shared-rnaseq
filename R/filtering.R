#' Filters low count reads
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts Count matrix or ExpressionSet
#' @param threshold Minimum number of reads a given sample has to have to
#'        contribute towards a gene passing the low count filter
#' @param min_samples Mininum number of samples that require to have at least
#'        `threshold` reads for a gene to pass the filter.
#' @param verbose Whether or not to print out information about the number of
#'        genes filtered.
#'
#' @return Matrix of filtered counts
filter_lowcount_genes <- function (counts, threshold=2, min_samples=2, verbose=TRUE) {
    num_before <- nrow(counts)
    if (class(counts) == "ExpressionSet") {
        keep <- rowSums(exprs(counts) > threshold) >= min_samples
    }
    else {
        keep <- rowSums(counts > threshold) >= min_samples
    }
    counts <- counts[keep, ]
    if (verbose) {
        print(sprintf("Removing %d low-count genes (%d remaining).", 
            num_before - nrow(counts), nrow(counts)))
    }
    return(counts)
}


#' Removes low variance genes from a count matrix
#'
#' Given a specified quantile, all genes with expression variance less than or
#' equal to that cutoff are filtered out of the count matrix.
#'
#' @param counts Input ExpressionSet
#' @param percentile_cutoff Variance pencentile to use as a cutoff.
#' @param include_plot Whether or not to include a plot showing which genes
#'        were removed.
#'
#' @return Filtered ExpressionSet.
filter_low_variance_genes <- function(counts, percentile_cutoff, include_plot=FALSE) {
    # If no variance filter is applied, don't do anything to counts
    if (percentile_cutoff == 0) {
        return(counts)
    }

    # Compute gene averages and variances
    gene_centers   <- rowMeans(exprs(counts))
    gene_variances <- apply(exprs(counts), 1, var)

    # Variance filtering mask
    variance_cutoff <- as.numeric(quantile(gene_variances, percentile_cutoff))
    low_variance_mask <- gene_variances <= variance_cutoff

    # Apply filtering to count tables at different processing stages
    counts <- counts[!low_variance_mask, ]

    # Show which genes were filtered out
    if (include_plot) {
        df <- data.frame(mean=gene_centers, variance=gene_variances,
                        filtered=low_variance_mask)
        ggplot(df, aes(x=log2(mean + 1), y=log2(variance + 1))) +
            geom_point(size=I(0.75), aes(colour=factor(!filtered))) +
            scale_fill_discrete(name="Passed variance filter") +
            #scale_colour_manual(values=c("gray", "blue")) +
            ggtitle("Filtered genes")
    }

    return(counts)
}

#' Removes multicopy genes
#'
#' Filters out parasite multicopy genes both by looking for genes which are
#' annotated as belonging to a known family (e.g. T. cruzi 'mucins'), and by
#' using an external list of multicopy genes generated based on ORF sequence
#' similarity.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts Count matrix
#' @param gene_info Gene annotation data frame
#' @param gene_description_filter A regular expression string indicating
#' description terms to exclude.
#' @param orf_counts_file A file containing a mapping from gene ids to the
#' number of genes with homologous ORFs found in the genome
#' @param orf_cutoff Minimum number of homologous ORFs for a gene to be
#' considered part of a multicopy family.
#' @param verbose Whether or not to print additional information about the
#' genes removed, ORF input structure, etc.
#'
#' @param Count matrix or ExpressionSet with multicopy genes removed.
filter_multicopy_genes <- function (counts, gene_info, gene_description_filter,
                                    orf_counts_file, orf_cutoff=5,
                                    verbose=FALSE) {
    # load counts file
    orf_counts <- tbl_df(read.delim(orf_counts_file, header=FALSE,
                                   col.names=c('gene', 'count', 'gene2')))

    # example of input table structure
    if (verbose) {
        cat("\n**ORF count file format:**\n")
        kable(head(orf_counts))
        cat("\n")
    }
    multicopy_ids <- (orf_counts %>% filter(count >= orf_cutoff))$gene

    # Just to be safe, let's also remove rows based on description
    multicopy_ids <- union(
        multicopy_ids,
        rownames(counts)[grepl(gene_description_filter,
                               gene_info$description)]
    )
    cat(sprintf("\nFound **%d** multicopy genes.\n", length(multicopy_ids)))

    multicopy_genes <- gene_info %>% filter(gene_id %in% multicopy_ids)

    multicopy_gene_classes <- multicopy_genes %>%
        distinct(description) %>%
        select(description) %>%
        arrange(description)

    if (verbose) {
        cat("Multicopy gene descriptions include:\n\n")
        kable(multicopy_gene_classes)
        cat("\n")
    }

    counts[!rownames(counts) %in% multicopy_ids,]
}

