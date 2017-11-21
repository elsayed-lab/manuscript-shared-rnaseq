#' Parses enrichment analysis results and returns a vector of items that had
#' at least one enriched annotation.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param results A dataframe containing one or more results from goseq.
#'
#' @return Returns a vector containing the names of all co-expression modules
#'         for which at least one annotation was found to be over- or
#'         under-represented.
get_enriched_modules <- function(results) {
    total_enriched <- 0

    # create a vector to keep track of which modules had enrichment
    enriched_colors <- c()

    for (result_name in names(results)) {
        result <- results[[result_name]]

        # determine number of over- and under-represented categories
        over_rep  <- result[result$over_represented_pvalue_adj < 0.05,]
        under_rep <- result[result$under_represented_pvalue_adj < 0.05,]

        if ((nrow(over_rep) + nrow(under_rep)) > 0) {
            enriched_colors <- append(enriched_colors, result_name)
        }

    }
    return(enriched_colors)
}

#' Returns a an empty enrichment test result
#'
#' @return A placeholder enrichment result dataframe
get_enrichment_placeholder <- function () {
    data.frame(category=numeric(0),                                           
               over_represented_pvalue=numeric(0),                       
               under_represented_pvalue=numeric(0),                      
               num_in_subset=numeric(0),                                 
               num_total=numeric(0),                                     
               over_represented_pvalue_adj=numeric(0),                   
               under_represented_pvalue_adj=numeric(0))                  
}

#' Writes a table summarizing enrichment for a set of genes.
#' 
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param results A dataframe containing one or more results from goseq.
#' @param module_output_dir Directory to store enrichment results in
#' @param enrichment_type Filename suffix to use to indicate type of enrichment
#'        performed.
#' @param annotation_mapping An optinal dataframe mapping from annotations
#'        to other useful information such as descriptions.
#' 
#' @return None
output_module_enrichment_results <- function(results, module_output_dir, 
                                             enrichment_type='go',
                                             annotation_mapping=NULL) {
    for (mod_color in names(results)) {
        # Create dataframe to store output rows
        output_df <- data.frame()

        x <- results[[mod_color]]

        # add placeholder entries for annotations not defined in mapping
        if (!is.null(annotation_mapping) && sum(!x$category %in% annotation_mapping$category) > 0) {
            missing <- x$category[!x$category %in% annotation_mapping$category]
            placeholders <- matrix(NA, 
                                   nrow=length(missing), 
                                   ncol=ncol(annotation_mapping))
            colnames(placeholders) <- colnames(annotation_mapping)
            placeholders[,'category'] <- missing
            
            annotation_mapping <- rbind(annotation_mapping, placeholders)
        }

        # determine number of over- and under-represented categories
        over_rep  <- x[x$over_represented_pvalue_adj < 0.05,]
        under_rep <- x[x$under_represented_pvalue_adj < 0.05,]

        if (nrow(over_rep) > 0) {
            out <- over_rep[,c('category', 'over_represented_pvalue_adj',
                                'num_in_subset', 'num_total')]
            colnames(out)[2] <- 'pvalue_adj'

            # add additional fields and print
            if (!is.null(annotation_mapping)) {
                out <- merge(out, annotation_mapping, by='category')
            }
            out <- cbind(out, direction='over-represented')
            output_df <- rbind(output_df, out)
        }

        if (nrow(under_rep) > 0) {
            out <- under_rep[,c('category', 'under_represented_pvalue_adj', 
                                'num_in_subset', 'num_total')]
            colnames(out)[2] <- 'pvalue_adj'

            # add additional fields and print
            if (!is.null(annotation_mapping)) {
                out <- merge(out, annotation_mapping, by='category')
            }
            out <- cbind(out, direction='under-represented')
            output_df <- rbind(output_df, out)
        }

        # write output
        if (nrow(output_df) > 0) {
            filename <- sprintf("%s_%s_enrichment.txt", mod_color, 
                                enrichment_type)
            write.table(output_df,
                        file.path(module_output_dir, filename),
                        row.names=FALSE, quote=FALSE, sep="\t")
        }
    }
}

#' Parses results from goseq enrichment analysis 
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param results A dataframe containing one or more results from goseq.
#' @param subset_sizes A dataframe mapping module/DE result names to number 
#'        of genes in that module or DE result.
#' @param annotation_name Type of annotation bring printed
#' @param annotation_mapping An optinal dataframe mapping from annotations
#'        to other useful information such as descriptions.
#' @param gene_mapping An optional dataframe mapping from gene IDs to
#'        annotations.
#' @param output_dir Directory where table output results will be saved to.
#' @param enrichment_type Filename suffix to use to indicate type of enrichment
#'        performed.
#' @param exclude_unclustered Whether or not the uncluster genes (grey) should
#'        be excluded from the results.
#' @param include_gene_lists If true, a list of genes associated with the
#'        enriched annotations will also be printed; requires that
#'        `gene_mapping` also be specified.
#'
#' @return None
print_enrichment_results <- function(results, subset_sizes, 
                                     annotation_name='annotations', 
                                     annotation_mapping=NULL,
                                     gene_mapping=NULL,
                                     output_dir=NULL,
                                     enrichment_type="",
                                     exclude_unclustered=TRUE,
                                     include_gene_lists=FALSE,
                                     str_max_width=Inf) {
    # counters
    total_enriched <- 0
    pvalues <- c()

    # For GO enrichment results, make GO: categories hyperlinks
    if(enrichment_type == 'go') {
        if (!is.null(gene_mapping)) {
            gene_mapping$category <- sprintf("<a href='http://amigo.geneontology.org/amigo/term/%s'>%s</a>",
                                             gene_mapping$category, gene_mapping$category)
        }
        if (!is.null(annotation_mapping)) {
            annotation_mapping$category <- sprintf("<a href='http://amigo.geneontology.org/amigo/term/%s'>%s</a>",
                                                   annotation_mapping$category, annotation_mapping$category)
        }
    }

    for (result_name in names(results)) {
        result <- results[[result_name]]

        # Skip unclustered genes if enabled
        if (exclude_unclustered && result_name == 'grey') {
            next
        }

        # Skip entries with no enrichment
        if (nrow(result) == 0) {
            next
        }

        # add additional fields
        if (!is.null(annotation_mapping)) {
            result <- merge(result, annotation_mapping, by='category')
        }

        # For GO enrichment results, make GO: categories hyperlinks
        if(enrichment_type == 'go') {
            result$category <- sprintf("<a href='http://amigo.geneontology.org/amigo/term/%s'>%s</a>",
                                       result$category, result$category)
        }

        # determine number of over- and under-represented categories
        over_rep  <- result[result$over_represented_pvalue_adj < 0.05,]
        under_rep <- result[result$under_represented_pvalue_adj < 0.05,]

        # keep track of the number of enriched annotations along with their
        # average adjusted p-value.
        total_enriched <- total_enriched + nrow(over_rep) + nrow(under_rep)

        # module size
        num_genes <- subset_sizes$num_genes[subset_sizes$module_id == result_name]

        # Print results
        if (nrow(over_rep) > 0 || nrow(under_rep) > 0) {
            # With download links
            if (!is.null(output_dir)) {
                # Table containing enriched annotations for module
                enrichment_filename <- sprintf("%s_%s_enrichment.txt", result_name, enrichment_type)
                enrichment_output <- file.path(output_dir, enrichment_filename)
            
                # Table containing information for genes in the module
				genes_output <- file.path(output_dir, sprintf("%s_genes.txt", result_name))

                # Output header with links to tables
                cat(sprintf("\n#### [%s](%s) enrichment (%d [genes](%s))\n", 
                            result_name, enrichment_output, num_genes,
                            genes_output))
            } else {
                # Output header without links to tables
                cat(sprintf("\n#### %s enrichment (%d genes)\n", result_name, num_genes))
            }
        }

        # Over-represented terms
        if (nrow(over_rep) > 0) {
            # fields to display
            out <- over_rep %>% select(-over_represented_pvalue,
                                       -under_represented_pvalue,
                                       -under_represented_pvalue_adj)

            cat(sprintf("\n**Over-represented %s:**\n", annotation_name))

            # print
            print(xkable(out %>% dplyr::rename(adj_pval=over_represented_pvalue_adj), 
                         str_max_width=str_max_width, row.names=FALSE))
            cat('\n')

            # Add adjusted pvalues to vector for averaging purposes
            pvalues <- append(pvalues, over_rep$over_represented_pvalue_adj)

            # Print specific genes responsible, if requested
            if (include_gene_lists) {
                cat("\n**Genes responsible for enrichment:**\n")             
                gene_list <- gene_mapping %>% filter(category %in% out$category & 
                                                     color==result_name)
                gene_list <- gene_list[!duplicated(gene_list),]
                                         
                print(xkable(gene_list, str_max_width=str_max_width, row.names=FALSE))
                cat('\n')
            }
        }

        # Under-represented terms
        if (nrow(under_rep) > 0) {
            # fields to display
            out <- under_rep %>% select(-under_represented_pvalue,
                                        -over_represented_pvalue,
                                        -over_represented_pvalue_adj)

            cat(sprintf("\n**Under-represented %s:**\n", annotation_name))

            print(xkable(out %>% dplyr::rename(adj_pval=under_represented_pvalue_adj),
                         str_max_width=str_max_width, row.names=FALSE))
            cat('\n')

            # Add adjusted pvalues to vector for averaging purposes
            pvalues <- append(pvalues, under_rep$under_represented_pvalue_adj)

            # Print specific genes responsible, if requested
            if (include_gene_lists) {
                cat("\n**Genes responsible for enrichment:**\n")             
                gene_list <- gene_mapping %>% filter(category %in% out$category & 
                                                     color==result_name)
                gene_list <- gene_list[!duplicated(gene_list),]
                                         
                print(xkable(gene_list, str_max_width=str_max_width, row.names=FALSE))
                cat('\n')
            }
        }
    }

    # Summary
    cat(sprintf("\n#### Total enriched %s\n", annotation_name))
    
    if (total_enriched > 0) {
        cat(sprintf("\nTotal: %d (Total -log10(adj.pval) = %f)\n", total_enriched,
                    sum(-log10(pmax(pvalues, 1E-10)))))
    } else {
        cat("\nTotal: 0\n")
    }
}

#' Compute summary statistics for enrichment analysis.
#'
#' Returns some basic statistics about the enrichment results across network
#' modules including the total number of enriched categories observed, the
#' number of unique enriched categories, the average p-value for enriched
#' categories, and the total number of modules with one or more enriched
#' categories.
#'
#' Note that this function currently only includes over-represented categories.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param enrichment_result List returned from call to test_gene_enrichment.
#'
#' @return List of enrichment analysis statistics
summarize_enrichment_result <- function(enrichment_result) {
    # Determine average p-value for most significantly enriched terms
    categories <- c()
    pvalues <- c()

    # keep track of the total number of modules with over-represented GO terms
    num_enriched_modules <- 0

    for (module_color in names(enrichment_result)) {
        # Skip grey module (unclustered genes)
        if (module_color == 'grey') {
            next
        }
        # Enriched categories for a single module
        module <- enrichment_result[[module_color]]
        
        # Keep only over-represented terms
        over_rep <- module[module$over_represented_pvalue_adj < 0.05,]

        categories <- append(categories, over_rep$category)

        # mean adj. p-val, and sum(-log10p(adj. p-vals)) for values <= 0.05;
        # clip at 1E-10 to avoid -Inf values.
        pvalues <- append(pvalues, pmax(over_rep$over_represented_pvalue_adj, 10^-10))

        # increment enriched module counter
        if (nrow(over_rep) > 0) {
            num_enriched_modules <- num_enriched_modules + 1
        }
    }

    mean_pval <- round(mean(pvalues), 4)
    pval_score <- sum(-log10(pvalues))

    # return result as a list
    return(list(
        num_enriched_modules = num_enriched_modules,
        total_categories     = length(pvalues),
        unique_categories    = length(unique(categories)),
        mean_pval            = mean_pval,
        pval_score           = pval_score
    ))
}


#' Function to check for GO enrichment for Human/Mouse
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param gene_subset Subset of gene ids to check for enrichment
#' @param all_genes All gene ids
#' @param gene_annotation_mapping Maping from gene ids to annotations (GO
#'                                terms, KEGG pathways, etc.)
#' @param gene_lengths Dataframe mapping gene ids to transcript lengths
#'
#' @return Dataframe containing enrichment results from a single GOSeq run.
test_host_go_enrichment <- function(gene_subset, all_genes,
                                    genome='hg19', id='ensGene') {

    # Subset membership status
    genes <- as.numeric(all_genes %in% gene_subset)
    names(genes) <- all_genes

    pwf <- nullp(genes, genome, id, plot.fit=FALSE)
    result <- goseq(pwf, genome, id)

    colnames(result) <- c('category', 'over_represented_pvalue',
                            'under_represented_pvalue', 'num_in_subset',
                            'num_total', 'term', 'ontology')

    result <- cbind(
        result,
        over_represented_pvalue_adj=p.adjust(result$over_represented_pvalue, method="BH"),
        under_represented_pvalue_adj=p.adjust(result$under_represented_pvalue, method="BH")
    )

    # Drop non-significant results
    result <- result[
        result$under_represented_pvalue_adj < 0.05 |
        result$over_represented_pvalue_adj  < 0.05,]
    
    return(result)
}


#' Function to check for annotation enrichment in gene subset
#'
#' NOTE: GOSeq does not support dplyr tbl_df instances for the gene2cat field
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param gene_subset Subset of gene ids to check for enrichment
#' @param all_genes All gene ids
#' @param gene_annotation_mapping Maping from gene ids to annotations (GO
#'                                terms, KEGG pathways, etc.)
#' @param gene_lengths Dataframe mapping gene ids to transcript lengths
#'
#' @return Dataframe containing enrichment results from a single GOSeq run.
test_gene_enrichment <- function(gene_subset, all_genes,
                                 gene_annotation_mapping, gene_lengths) {
                                
    # Subset membership status
    member_status <- as.numeric(all_genes %in% gene_subset)
    names(member_status) <- all_genes

    # Relationship between gene length and subset membership
    pwf <- nullp(member_status, bias.data=gene_lengths, plot.fit=FALSE)

    # Check to make sure annotation mapping is in the expected format
    if (all(!c('gene', 'category') %in% colnames(gene_annotation_mapping))) {
        error("Gene-annotation mapping does not have expected fields")
    }

    # Check to see if we are performing GO enrichment, or some other
    # type of annotation enrichment
    go_enrichment <- grepl('GO:', gene_annotation_mapping$category[1])

    # GO enrichment
    if (go_enrichment) {
        # Note: test.cats parameter is ignored when a gene/term mapping is
        # specified manually
        result <- goseq(pwf, gene2cat=gene_annotation_mapping)

        colnames(result) <- c('category', 'over_represented_pvalue',
                              'under_represented_pvalue', 'num_in_subset',
                              'num_total', 'term', 'ontology')
    } else {
        # Other annotation enrichment
        result <- goseq(pwf, gene2cat=gene_annotation_mapping)

        colnames(result) <- c('category', 'over_represented_pvalue',
                              'under_represented_pvalue', 'num_in_subset',
                              'num_total')
    }

    # Adjust for multiple testing using Benjamini and Hochberg method
    result <- cbind(result,
                   over_represented_pvalue_adj=p.adjust(
                       result$over_represented_pvalue, method="BH"),
                   under_represented_pvalue_adj=p.adjust(
                       result$under_represented_pvalue, method="BH")
                  )

    # Drop non-significant results
    result <- result[
        result$under_represented_pvalue_adj < 0.05 |
        result$over_represented_pvalue_adj  < 0.05,]
    
    return(result)
}

#' Writes a GMT version of a gene/annotation mapping to file.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param mapping Dataframe mapping from gene IDs to gene sets
#' @param outfile string Filepath to save GMT output to
#'
#' @return None
gene_set_mapping_to_gmt <- function(mapping, outfile) {
    # open output file for writing
    fp <- file(outfile, 'w')

    # output rows
    entries <- c()

    # iterate over gene sets
    categories <- unique(mapping$category)

    for (i in seq_along(categories)) {
        category <- categories[i]

        # collapse gene set entries into a single tab-delimited string
        gene_entries <- paste0(mapping$gene[mapping$category == category],
                               collapse='\t') 

        # create row in GMT file
        entries <- c(entries, paste0(c(category, 'na', gene_entries), collapse='\t'))
    }

    # save output
    writeLines(entries, fp)

    close(fp)
}

