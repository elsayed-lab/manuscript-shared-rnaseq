#g su' Create a dataframe of module statistics
#'
#' Constructs a summary dataframe containing some useful module statistics such
#' as the number of enriched GO terms in the module, the average correlation of
#' gene expression profiles within the module, etc.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param result Data frame of co-expression analysis results.
#' @param similarity_matrix Co-expression similarity matrix
#' @param wgcna_input Original counts used as input to WGCNA
#' @param main_contrast String specifying primary DE contrast of interest
#' @param ... List of all relevant enrichment results
#'
#' @return Data frame containing module statistics such as median count,
#'              average total expression, etc.
create_module_stats_df <- function(result, similarity_matrix, wgcna_input,
                                   main_contrast, ...) {
    # Module statistics data frame
    # Using dplyr standard-evaluation mode to allow for dynamic variables
    # See http://cran.r-project.org/web/packages/dplyr/vignettes/nse.html
    if (any(grepl('de_', colnames(result)))) {
        # including differential expression results
        module_stats <- result %>% 
                    group_by(color) %>% 
                    summarise_('n()',
                            sprintf('sum(de_%s)', main_contrast),
                            sprintf('round(sum(de_%s) / n(), 2)', main_contrast),
                            mean_expr='round(mean(expr_mean), 2)', 
                            mean_gene_var='round(mean(expr_variance),2)',
                            mean_tx_len='round(mean(transcript_length, na.rm=TRUE), 2)'
                    )
        colnames(module_stats) <- c('color', 'num_genes', 'num_de', 'ratio_de',
                                    'mean_expr', 'mean_gene_var', 'mean_tx_len')
    } else {
        # including differential expression results
        module_stats <- result %>% 
                    group_by(color) %>% 
                    summarise_('n()',
                            mean_expr='round(mean(expr_mean), 2)', 
                            mean_gene_var='round(mean(expr_variance),2)',
                            mean_tx_len='round(mean(transcript_length), 2)'
                    )
        colnames(module_stats) <- c('color', 'num_genes', 'mean_expr',
                                    'mean_gene_var', 'mean_tx_len')
    }
    # Average gene correlation within the module
    min_correlation     <- c()
    max_correlation     <- c()
    mean_correlation    <- c()
    median_correlation  <- c()
        
    for (color in module_stats$color) {
        genes <- result[result$color == color,]$gene_id
        cor_matrix <- similarity_matrix[
            rownames(similarity_matrix) %in% genes,
            colnames(similarity_matrix) %in% genes
        ]

        min_correlation  <- append(min_correlation, round(min(cor_matrix), 2))
        max_correlation  <- append(max_correlation, round(max(cor_matrix), 2))
        mean_correlation <- append(mean_correlation, round(mean(cor_matrix), 2))
        median_correlation <- append(median_correlation, round(median(cor_matrix), 2))
    }

    # Intra-module variance
    intra_module_var <- c()

    for (color in module_stats$color) {
        genes <- result[result$color == color,]$gene_id
        expr <- wgcna_input[rownames(wgcna_input) %in% genes,]
        intra_module_var <- append(intra_module_var, round(mean(var(t(expr))), 2))
    }

    # Create result data frame
    df <- cbind(module_stats, min_correlation, max_correlation,
                    mean_correlation, median_correlation, intra_module_var)

    # Add GO/KEGG/etc. enrichment info
    for (enrichment_type in names(...)) {
        df <- cbind(df, as.vector(sapply(...[[enrichment_type]], nrow)))
        cname <- sprintf("num_enriched_%s", enrichment_type)
        colnames(df) <-  c(colnames(df)[1:ncol(df) - 1], cname)
    }

    # Add total enriched column
    df$num_enriched_total <-  df %>% 
        select(starts_with("num_enriched")) %>% 
        rowSums 

    # Return result, ordered by number of enrichment results
    df %>% arrange(desc(num_enriched_total))
}

#' Outputs separate GraphML files for each co-expression module in the network.
#' 
#' @param adj_mat Co-expression network adjacency matrix
#' @param result Co-expression analysis result
#' @param output_dir Directory to save GraphML output to
#' @param weighted Whether or not to include edge weights in the output
#' @param threshold Edge weight cutoff; when saving a weighted network, only
#'        those edges with weight >= that threshold will be included.
#' @param threshold For weighted networks, if a threshold value between 0 and 
#'        1 is specified, all edges with weights below that value with be 
#'        dropped from the graph.
#' @param max_edge_ratio The maximum number of edges per node in the network to
#'        allow. If the number of edges that would remain for the specified 
#'        threshold exceeds this value, the threshold will be raised to reduce 
#'        the number of edges remaining.
#' @param nodeAttrDataFrame A data frame containing one or more columns 
#'        associated with the vertices in the graph.  The ith row of the data
#'        frame should correspond to the ith entry in the adjacency matrix.
#'
#' @return None
export_module_subnetworks <- function(adj_mat, result,  output_dir="output",
                                      weighted=TRUE, threshold=0.5,
                                      max_edge_ratio=3, nodeAttrDataFrame=NULL) {
    # Create output directory if it doesn't already exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }
    
    # Iterate over network modules and generate a GraphML output file for each
    # one
    for (color in module_stats$color) {
        # get sub-network
        genes <- result[result$color == color,]$gene_id
        module_indices <- rownames(adj_mat) %in% genes

        module_adjacency_matrix <- adj_mat[module_indices, module_indices]
        module_node_attr_dataframe <- nodeAttrDataFrame[module_indices,]

        # filename to use
        outfile <- file.path(output_dir, sprintf("%s.graphml", color))

        g <- export_network_to_graphml(module_adjacency_matrix,
                                       filename=outfile, weighted=weighted,
                                       threshold=threshold,
                                       max_edge_ratio=max_edge_ratio,
                                       nodeAttrDataFrame=module_node_attr_dataframe)
    }
}

#' Converts an adjaceny matrix along with some optional vertex and edge
#' information to a GraphML graph and saves it to disk.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param adj_mat An n-by-n weighted or unweighted adjacency matrix normalized
#' to contain values between 0 and 1.
#' @param filename Name of file to save output to. If file already exists it
#' will be overwritten. (default: network.graphml)
#' @param weighted Whether or not the adjacency matrix should be treated as a 
#' weighted graph. (default: TRUE)
#' @param threshold For weighted networks, if a threshold value between 0 and 
#' 1 is specified, all edges with weights below that value with be dropped from
#' the graph. (default: 0.5)
#' @param max_edge_ratio The maximum number of edges per node in the network to
#' allow. If the number of edges that would remain for the specified threshold
#' exceeds this value, the threshold will be raised to reduce the number of
#' edges remaining. (default: 3)
#' @param nodeAttr A vector with length equal to the number of vertices in the 
#' network, where the ith entry in the vector corresponds to some numeric or 
#' string annotation that should be associated with the ith node in the 
#' adjacency matrix. (default: NULL)
#' @param nodeAttrDataFrame A data frame containing one or more columns 
#' associated with the vertices in the graph.  The ith row of the dataframe 
#' should correspond to the ith entry in the adjacency matrix. (default: NULL)
#' @param edgeAttributes Extra attributes to associate with the graph edges,
#' formatted as a list of matrices of the same dimension and names as the
#' adjacency matrix.
#' @param verbose If true, outputs information about number of nodes and edges
#' included in the output network.
#'
#' Examples
#' --------
#' export_network_to_graphml(sim_matrix, filename='~/network.graphml',
#'                           threshold=0.3, nodeAttrDataFrame=df)
#'
#' @seealso http://www.inside-r.org/packages/cran/WGCNA/docs/exportNetworkToCytoscape
#' @seealso http://graphml.graphdrawing.org/
#'
#' 
#' @return igraph graph object representing the exported graph.
export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
    library('igraph')

    # Determine filename to use
    if (is.null(filename)) {
        filename='network.graphml'
    }

    # TODO 2015/04/09
    # Add option to rescale correlations for each module before applying
    # threshold (this is simpler than the previous approach of trying to
    # determine a different threshold for each module)
    #
    # Still, modules with very low correlations should be given somewhat
    # less priority than those with very high correlations.

    #module_colors <- unique(nodeAttrDataFrame$color)
    #module_genes <- which(nodeAttrDataFrame$color == color)
    #module_adjmat <- adj_mat[module_genes,]
    #num_genes <- length(module_genes)

    # Adjust threshold if needed to limit remaining edges
    max_edges <- max_edge_ratio * nrow(adj_mat)

    edge_to_total_ratio <- max_edges / length(adj_mat)
    edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))

    # Also choose a minimum threshold to make sure that at least some edges
    # are left
    min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))

    threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))

    # Remove edges with weights lower than the cutoff
    adj_mat[abs(adj_mat) < threshold] <- 0

    # Drop any genes with no edges (TODO: Make optional?)
    orphaned <- (colSums(adj_mat) == 0)
    adj_mat <- adj_mat[!orphaned, !orphaned]

    # Also remove annotation entries
    if (!is.null(nodeAttr)) {
        nodeAttr <- nodeAttr[!orphaned]
    }
    if (!is.null(nodeAttrDataFrame)) {
        nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
    }

    # Keep track of non-positive edges and rescale to range 0,1
    is_zero     <- adj_mat == 0
    is_negative <- adj_mat < 0

    adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
    adj_mat[is_zero] <- 0
    adj_mat[is_negative] <- -adj_mat[is_negative]

    if (verbose) {
        message(sprintf("Outputting matrix with %d nodes and %d edges", 
                        nrow(adj_mat), sum(adj_mat > 0)))
    }

    # Create a new graph and add vertices
    # Weighted graph
    if (weighted) {
        g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
    } else {
        adj_mat[adj_mat != 0] <- 1
        g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
    }

    # Add single node annotation from vector
    if (!is.null(nodeAttr)) {
        g <- set.vertex.attribute(g, "attr", value=nodeAttr)
    }

    # Add node one or more node annotations from a data frame
    if (!is.null(nodeAttrDataFrame)) {
        for (colname in colnames(nodeAttrDataFrame)) {
            g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
        }
    }

    # TODO 2014/06/14
    # Add support for setting edge attributes...
    # Convert matrix of edge annotations to a list of edge pairs and values
    #x <- as.data.frame(as.table(diff))
    #E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
    #E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
    edge_correlation_negative <- c()

    # neg_correlations[edge_list]
    edge_list <- get.edgelist(g)

    for (i in 1:nrow(edge_list)) {
        from <- edge_list[i, 1]    
        to   <- edge_list[i, 2]    
    }
    
    # Apply threshold

    # Save graph to a file
    write.graph(g, filename, format='graphml')

    # return igraph
    return(g)
}

#' Determines index of module x in the dendrogram
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param color module color
#' @param module_colors gene module assignments
#' @param gene_tree clustering dendrogram
#'
#' @return Index of specified module in a dendrogram.
get_module_indices <- function(color, module_colors, gene_tree) {
    return(which(module_colors[gene_tree$order] == color))
}

#' Determine the order of modules in the dendrogram
#' 
#' For hybrid clusterings where module assignments are not strictly continuous,
#' this will provide an approximation of the module order in the dendrogram.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param module_colors gene module assignments
#' @param gene_tree clustering dendrogram
#'
#' @return Return module colors in approximately the order they appear in the
#' dendrogram, moving from left to right.
get_module_order <- function(module_colors, gene_tree) {
    unique_colors <- unique(module_colors)

    # get median index for each module color; this should be a good
    # approximation of the modules positions relative to each other
    x <- data.frame()
    for (i in unique_colors) {
        median_index <- median(get_module_indices(i, module_colors, gene_tree))
        x <- rbind(x, data.frame(color=i, median=median_index))
    }

    # sort by median index
    x <- x[order(x$median),]

    # return sorted module colors
    return(x$color)
}

#' Convert labels to colors for plotting
#' 
#' The WGCNA labels2colors is currently limited to generating 436 non-redundant
#' colors; any results with greater than 435 modules (plus gray for unassigned)
#' may not be reliable. This function provides the ability to fall back on a
#' wider range of colors.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param labels Vector of module labels (numbers).
#'
#' @return Color palette to use for co-expression modules.
labels2colors <- function(labels) {
    num_labels <- length(unique(labels))
    # If less than 436, use WGCNA's function to generate colors
    if (num_labels <= 435) {
        return(WGCNA::labels2colors(labels))            
    }
    # Otherwise create a larger color palette
    pal <- c('grey', substring(rainbow(num_labels), 1, 7))
    return (pal[labels + 1])
}

#' Merge modules with similar eigengenes
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param sample_counts Count matrix
#' @param module_colors gene module assignments
#' @param merge_cor Minimum correlation required for two modules to be combined
#' @param verbose Whether or not to print information about the number of
#' modules before and after merging has been performed.
#' 
#' @return List of updated module labels, colors, and eigengenes are merging
#' has been performed.
merge_similar_modules <- function(sample_counts, module_colors, merge_cor, verbose=FALSE) {
    # Call an automatic merging function
    merge_result <- mergeCloseModules(sample_counts, module_colors,
                                     cutHeight=1 - merge_cor,
                                     verbose=ifelse(isTRUE(verbose), 3, 0))

    # The merged module colors and eigengenes
    merged_colors     <- merge_result$colors
    merged_eigengenes <- merge_result$newMEs

    num_modules_before <- length(unique(module_colors))
    print(sprintf("Number of modules before merging: %d", num_modules_before))

    # total number of modules remaining after merge
    num_modules_after <- length(unique(merged_colors))
    print(sprintf("Number of modules after merging: %d", num_modules_after))

    # Replace previous module assignments
    # Note: need to include -1 to ensure comparable color interpolation.
    # (Question: will this work if original tree did not include any 
    #  unassigned genes?)
    colormap_old <- labels2colors(0:(length(unique(module_labels))- 1))
    merged_labels_gapped <- match(merged_colors, colormap_old) - 1

    # Reassign colors
    merged_labels <- as.numeric(as.factor(merged_labels_gapped)) - 1
    merged_colors <- labels2colors(merged_labels)

    # Numbers corresponding to new eigengenes using the old scale
    merged_eigengene_labels_gapped <- match(sub('ME', '', names(merged_eigengenes)),
                                           colormap_old) - 1
    names(merged_eigengenes) <- paste0('ME', labels2colors(as.numeric(
                                as.factor(merged_eigengene_labels_gapped)) -1))

    return(list(
        module_labels=merged_labels,
        module_colors=merged_colors,
        module_eigengenes=merged_eigengenes)
    )
}


#' Plots an eigengene dendrogram
#' 
#' @param module_eigengenes Co-expression module eigengenes
#' @param module_colors gene module assignments
#' @param merge_cor If less than one, modules with correlation greater than or
#' qual to this value will be combined before plotting.
#'
#' @return None
plot_eigengene_dendrogram <- function (module_eigengenes, module_colors, merge_cor=1) {
    # Calculate dissimilarity of module eigengenes
    module_eigengene_dissim <- 1 - cor(module_eigengenes)

    # Cluster module eigengenes
    module_eigengene_tree <- flashClust(as.dist(module_eigengene_dissim),
                                    method="average")
    dendro <- as.dendrogram(module_eigengene_tree)

    # function to get color labels
    # http://gastonsanchez.com/blog/how-to/2012/10/03/Dendrograms.html
    # http://stat.ethz.ch/R-manual/R-patched/library/stats/html/dendrapply.html
    color_labels <- function(n) {
        if (is.leaf(n)) {
            a <- attributes(n)
            #labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
            color <- substring(a$label, 3)
            attr(n, "nodePar") <- c(a$nodePar, lab.col=color)
        }
        return(n)
    }

    # using dendrapply
    eigengene_dendro <- dendrapply(dendro, color_labels)

    # Module merge correlation
    combine_threshold <- 1 - merge_cor

    # Plot the result
    if (length(module_eigengene_tree$order) > 2) {
        plot(module_eigengene_tree, main="Clustering of module eigengenes",
            xlab="", sub="", cex=0.55)
        if (combine_threshold > 0) {
            # Plot the cut line into the dendrogram
            abline(h=combine_threshold, col="red")
        }
        # Alternate way of visualizing dendrogram
        plot(eigengene_dendro, main="Clustering of module eigengenes",
            xlab="", sub="", cex=0.55)
    }
}
