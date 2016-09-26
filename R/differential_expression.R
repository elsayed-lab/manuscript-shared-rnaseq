#' Performs all pairwise differential expression comparisons possible for a
#' given set of conditions using Limma.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param comparisons list of contrasts to run
#' @param counts Count matrix or voom output
#' @param design Design matrix
#' @param pvalue P-value cutoff for significant genes
#' @param log2fc Minimum log2-FC cutoff
#' @param robust Whether or not to use the eBayes `robust` option
#'
#' @return List of limma topTable results and significance-filtered topTable
#'         results for each input contrast
run_pairwise_limma_analyses <- function(comparisons, counts, design, 
                                        max_pvalue=0.05, min_log2fc=0,
                                        robust=FALSE) {
    # Fit a linear model for each gene
    fit <- lmFit(counts, design=design)

    # Number of genes
    num_genes <- nrow(fit$coefficients)

    # eBayes finds an F-statistic from the set of t-statistics for that gene
    design_colnames <- colnames(design)

    i <- 1
    top_genes <- list()
    sig_genes <- list()
    comparison_names <- c()
    num_de <- c()

    for (conds in comparisons) {
        # Reverse order of conditions so comparison subtracts the earlier
        # timepoint from the later one. This assumes that the comparisons
        # have been passed in in pairs of (early, late).
        contrast_formula <- paste(paste0("condition", rev(conds)),
                                 collapse=' - ')

        # Build contrast matrix
        contrast_mat <- makeContrasts(
            #comparison=as.formula(contrast_formula, env=.env),
            contrasts=contrast_formula,
            levels=design
        )
        contrast_fit <- contrasts.fit(fit, contrast_mat)

        eb <- suppressWarnings(eBayes(contrast_fit, robust=robust))

        # Grab highly ranked genes
        top_genes[[i]] <- topTable(eb, number=num_genes)

        # Convert to dplyr df
        top_genes[[i]] <- tbl_df(cbind(top_genes[[i]],
                                      gene_id=row.names(top_genes[[i]])))

        # Filter by p-value and log2-foldchange
        sig_genes[[i]] <- top_genes[[i]] %>% 
            filter(adj.P.Val < max_pvalue & abs(logFC) >= min_log2fc)

        # Check to see if some significant genes were found
        contrast_name <- sprintf("%s vs. %s", conds[2], conds[1])

        comparison_names <- append(comparison_names, 
                                paste(conds[2], conds[1], sep='_vs_'))
        num_de <- append(num_de, nrow(sig_genes[[i]]))

        # increment list counter and store entry name
        i <- i + 1
    }
    names(top_genes) <- comparison_names
    names(sig_genes) <- comparison_names

    return (list(top_genes=top_genes, sig_genes=sig_genes, num_de=num_de))
}

#' Generates an MA plot for a differential expression analysis result
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @seealso \code{\link{run_pairwise_limma_analyses}}
#'
#' @param counts Count table used in differential expression analysis
#' @param top_genes topTable result from limma run
#' @param sig_genes Significance filtered topTable result from limma run
#'
#' @return None
maplot <- function(counts, top_genes, sig_genes) {
    # MA plot
    subtitle= paste("Number of significant DE genes:", nrow(sig_genes), "/",
                    nrow(counts))
    plot(rowMeans(exprs(counts)[top_genes$gene_id,]),
        top_genes$logFC, pch=16, cex=0.5,
        col="darkgrey", main="MA Plot",
        ylab="log FC", xlab="Average Expression", sub=subtitle)

    # Plot any significant DE genes in red
    if (nrow(sig_genes) > 0) {
        points(
            rowMeans(exprs(counts[featureNames(counts)
                                        %in% sig_genes$gene_id,])),
            sig_genes$logFC, col="red", cex=0.5)
    }
    abline(h=c(-1,0,1), col="red")
}
