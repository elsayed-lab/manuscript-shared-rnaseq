#' Loads parasite gene annotations
#'
#' Creates a data frame (tibble) of gene information for a given set of gene
#' ID's using the OrganismDbi interface.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param fields Columns to include in the output
#'
#' @return Tibble of parasite gene annotations including: gene_id, chromosome,
#' description, strand, and type.
load_parasite_annotations <- function(orgdb, gene_ids, keytype='ENSEMBL',
                                      fields=c('CHR', 'GENEDESCRIPTION', 'TXSTRAND', 'TYPE')) {

    # Work-around 2017/02/16
    # With recent versions of GenomicFeatures, TXTYPE is used in place of TYPE
    if ('TYPE' %in% columns(orgdb)) {
        type_var <- 'TYPE'
    } else {
        type_var <- 'TXTYPE'
        fields <- sub('TYPE', 'TXTYPE', fields)
    }

    # Gene info
    # Note querying by "GENEID" will exclude noncoding RNAs
    gene_info <- AnnotationDbi::select(orgdb,
                                       keys=gene_ids,
                                       keytype=keytype,
                                       columns=fields)

    # Convert to tbl_df and reorganize
    tbl_df(gene_info) %>%
        select(
            gene_id     = .data[[keytype]],
            chromosome  = CHR,
            description = GENEDESCRIPTION,
            strand      = TXSTRAND,
            type        = .data[[type_var]]
        )
}

#' Load host gene annotations
#'
#' Creates a data frame (tibble) of gene information for a given set of gene
#' ID's using the OrganismDbi interface.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param fields Columns to include in the output
#' @param biomart_dataset Name of the biomaRt dataset to query for gene type
#' @param biomart_host URI for the biomaRt host to connect to
#'
#' @return Tibble of parasite gene annotations including: gene_id, chromosome,
#' description, strand, and type.
load_host_annotations <- function(orgdb, gene_ids, keytype='ENSEMBL',
                                  fields=c('TXCHROM', 'GENENAME', 'TXSTRAND'),
                                  biomart_dataset='hsapiens_gene_ensembl',
                                  biomart_host='www.ensembl.org') {
    # Gene info
    # Note querying by "GENEID" will exclude noncoding RNAs
    gene_info <- AnnotationDbi::select(orgdb, keys=gene_ids, keytype=keytype,
                                       columns=fields)

    # Remove any entries which are either deprecated (e.g. ENSG00000006074) or
    # missing transcript information (e.g. ENSG00000050327).
    #gene_info <- gene_info[!is.na(gene_info$TXSTART),]

    # Convert to tbl_df and reorganize
    gene_info <- tbl_df(gene_info) %>%
        select(
            gene_id     = .data[[keytype]],
            chromosome  = TXCHROM,
            description = GENENAME,
            strand      = TXSTRAND
        )

    # Get gene biotype from biomaRt
    ensembl_mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                            dataset=biomart_dataset, host=biomart_host)
    biomart <- useDataset(biomart_dataset, mart=ensembl_mart)

    biomart_genes <- getBM(attributes=c("ensembl_gene_id", "gene_biotype"), 
                        mart=biomart)
    gene_info$type <- biomart_genes$gene_biotype[match(gene_info$gene_id,
                                                biomart_genes$ensembl_gene_id)]

    return(gene_info)
}

#' Retrieves GO terms associated with a set of genes
#' 
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#' @param keytype Identifier type
#' @param include_ancestors If true, ancestor annotations will be included in
#' mapping.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @return Tibble containing a mapping from gene ID's to GO terms, along with
#' some additional information about each GO term.
load_go_terms <- function(orgdb, gene_ids, keytype='ENSEMBL', include_ancestors=TRUE) {
    # Each gene may be associated with multiple GO terms (1:many)
    go_terms <- AnnotationDbi::select(orgdb, keys=gene_ids, keytype=keytype,
                                      columns=c("GO"))

    # Remove EVIDENCE/ONTOLOGY fields if they are present to normalize handling
    # of GO annotations (ONTOLOGY will be added back in later on); Host queries
    # usually included these fields while parasite databases do not.
    go_terms <- go_terms[,c(keytype, 'GO')]

    # Drop genes with no associated GO terms
    go_terms <- go_terms[complete.cases(go_terms),]

    # Deduplicate (some terms have multiple evidence types)
    go_terms <- go_terms[!duplicated(go_terms),]

    # Get ancestors for each term
    if (include_ancestors) {
        # Requires GO.db to be loaded
        bp <- as.list(GOBPANCESTOR)
        mf <- as.list(GOMFANCESTOR)
        cc <- as.list(GOCCANCESTOR)

        for (i in 1:nrow(go_terms)) {
            gene_id   <- go_terms[i, keytype]
            go_id <- go_terms[i, 'GO']

            # determine which ontology term is associated with and retrieve
            # ancestors accordingly. 
            ancestors = NULL

            # this can take a while...
            if (i %% 1000 == 0) {
                message(sprintf("Retrieving ancestor terms for entry %d/%d", i, nrow(go_terms)))
            }

            if (go_id %in% names(bp)) {
                ancestors <- bp[[go_id]]
            } else if (go_id %in% names(mf)) {
                ancestors <- mf[[go_id]]
            } else if (go_id %in% names(cc)) {
                ancestors <- cc[[go_id]]
            }

            # In some cases, there may be GO terms which have not yet been
            # included in GO.db and therefore will be skipped.
            if (!is.null(ancestors)) {
                rows <- cbind(gene_id, ancestors)
                colnames(rows) <- c(keytype, 'GO')
                go_terms <- rbind(go_terms, rows)
            }
        }
    }

    # Drop genes with no associated GO terms
    go_terms <- go_terms[complete.cases(go_terms),]

    # Deduplicate ancestor terms and drop "all" entries
    go_terms <- go_terms[!duplicated(go_terms),]
    go_terms <- go_terms[!go_terms$GO == 'all',]

    # OrgDb query results may already include an ONTOLOGY field
    fields <- c('TERM', 'GOID')

    if (!'ONTOLOGY' %in% colnames(go_terms)) {
        fields <- append(fields, 'ONTOLOGY')
    }

    go_term_names <- AnnotationDbi::select(GO.db, keys=unique(go_terms$GO),
                                           columns=fields)

    go_terms <- merge(go_terms, go_term_names, by.x='GO', by.y='GOID')

    # remove EVIDENCE field if it exists
    go_terms <- go_terms[,colnames(go_terms) != 'EVIDENCE']

    # Remove redundant annotations which differ only in source/evidence
    go_terms <- unique(tbl_df(go_terms) %>% na.omit())

    return(go_terms)
}

#' Createss a gene/KEGG mapping data frame
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#'
#' @return Data frame mapping gene identifiers to KEGG pathways.
load_kegg_mapping <- function(orgdb, gene_ids, keytype='ENSEMBL') {
    kegg_mapping <- tbl_df(AnnotationDbi::select(orgdb, keys=gene_ids,
                                                 keytype=keytype,
                                                 columns=c('KEGG_PATH'))) %>% na.omit()
    colnames(kegg_mapping) <- c('gene', 'category')

    # goseq does not support tbl_df instances
    return(as.data.frame(kegg_mapping))
}

#' Creates a KEGG pathway/description mapping data frame
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param orgdb An OrganismDb instance
#' @param gene_ids Identifiers of the genes to retrieve annotations for
#'
#' @return Tibble containing names and descriptive information for each KEGG
#' pathway known for an organism.
load_kegg_pathways <- function(orgdb, gene_ids, keytype='ENSEMBL') {
    kegg_pathways <- suppressWarnings(
        tbl_df(AnnotationDbi::select(orgdb, keys=gene_ids,
                            keytype=keytype,
                            columns=c('KEGG_PATH', 'KEGG_NAME', 'KEGG_CLASS', 'KEGG_DESCRIPTION')))
    )
    kegg_pathways <- kegg_pathways %>%
        na.omit() %>%
        select(KEGG_PATH, KEGG_NAME, KEGG_CLASS, KEGG_DESCRIPTION)

    colnames(kegg_pathways) <- c('pathway', 'name', 'class', 'description')

    return(kegg_pathways)
}

#' Maps KEGG identifiers to ENSEMBL, etc. gene identifiers
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' Takes a list of KEGG gene identifiers and returns a list of ENSEMBL or other
#' identifiers corresponding to those genes. 
#'
#' @param kegg_ids List of KEGG identifiers to be mapped
#' @param keytype Annotation type to map to (e.g. ENSEMBL, FlyBase, etc.)
#'
#' @return Vector of mapped gene identifiers for the input set of KEGG gene ID's.
convert_kegg_ids <- function(kegg_ids, keytype='ENSEMBL') {
    # query gene ids 10 at a time (max allowed)
    result <- c()

    # Query each set of genes, three at a time
    for (x in split(kegg_ids, ceiling(seq_along(kegg_ids) / 3))) {
        print(x)
        query <- keggGet(x)
        for (item in query) {
            dblinks <- item$DBLINKS
            gene_id <- dblinks[grepl(keytype, dblinks, ignore.case=TRUE)]
    
            if (length(gene_id) > 0) {
                result <- append(result, unlist(strsplit(gene_id, ': '))[2])

                # Warn when there are one to many mappings
                if (length(gene_id) > 1) {
                    warning(sprintf("One to many KEGG mapping for gene %s", x))
                }
            }
        }
    }
    return(result)
}

#' Generates a GENE/KEGG mapping 
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param pathways Vector of KEGG pathway IDs returned from call to keggLink()
#'        e.g. "path:mmu05134"
#' @param org_abbreviation KEGG identifier for the species of interest (e.g.
#'        "hsa" for Homo sapiens)
#' @param keytype Annotation type to map to (e.g. ENSEMBL, FlyBase, etc.)
#'
#' @return Data frame mapping genes ID's to KEGG pathways
generate_gene_kegg_mapping <- function(pathways, org_abbreviation, 
                                       keytype='ENSEMBL', verbose=FALSE) {
    # data frame to store kegg gene mapping 
    kegg_mapping <- data.frame()

    # Iterate over pathways and query genes for each one
    for (pathway in pathways) {
        if (verbose) {
            message(sprintf("Processing for pathway %s", pathway))
        }
        
        # Get genes in pathway
        kegg_ids <- as.character(keggLink(org_abbreviation, pathway))

        # Convert gene IDS
        if (keytype != 'WBGene') {
            gene_ids <- convert_kegg_ids(kegg_ids, keytype)
        } else {
            gene_ids <- substring(kegg_ids, 10)
        }

        kegg_mapping <- unique(rbind(kegg_mapping,
                               data.frame(gene=gene_ids, pathway=pathway)))
    }

    return(kegg_mapping)
}

#' Generates a KEGG PATHWAY / description mapping
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param pathways Vector of KEGG pathway identifiers
#'
#' @return Data frame containing KEGG pathway identifiers along with pathway
#' class and descriptions information.
generate_kegg_pathway_mapping <- function(pathways, verbose=FALSE) {
    # result data frame
    kegg_pathways <- data.frame()

    for (pathway in pathways) {
        if (verbose) {
            message(sprintf("Processing for pathway %s", pathway))
        }

        # Get pathway info
        meta <- keggGet(pathway)[[1]]
        pathway_desc  <- ifelse(is.null(meta$DESCRIPTION), '', meta$DESCRIPTION)
        pathway_class <- ifelse(is.null(meta$CLASS), '', meta$CLASS)
        kegg_pathways <- rbind(kegg_pathways, 
                              data.frame(pathway=pathway,
                                         name=meta$PATHWAY_MAP,
                                         class=pathway_class,
                                         description=pathway_desc))
    }
    return(kegg_pathways)
}

#
# Parses GMT file at a specified location and returns a gene/annotation mapping dataframe
#
parse_gmt <- function(infile, keytype) {
  # determine maximum number columns (empty column in gmt file may be skipped, so adding
  # one to be safe)
  max_cols <- max(count.fields(infile, sep='\t'), na.rm = TRUE) + 1

  cnames <- c('category', 'source', paste0('gene_', 1:(max_cols - 2)))

  # read in table, filling empty cells with NA's
  gmt <- read.delim(infile, sep = '\t', header = FALSE, col.names = cnames,
                    fill = TRUE, stringsAsFactors = FALSE)

  # drop any empty columns after the first two
  na_cols <- which(apply(gmt, 2, function(x) { all(is.na(x)) }))
  na_cols <- na_cols[!na_cols %in% 1:2]

  if (length(na_cols) > 0) {
    gmt <- gmt[, -na_cols]
  }

  # Check first column of genes to see if weights are specified for each entry
  #  (e.g. <SYMBOL>,1.0)
  if (all(grepl(',', gmt$gene_1))) {
    # for now, we will simply discard the weights and keep the annotation - term mapping 
    for (i in 3:ncol(gmt)) {
        message(i)
        gmt[, i] <- unlist(lapply(strsplit(gmt[, i], ','), '[', 1))
    } 
  }

  # remove any annotation with no associated gene entries (found in some gmt files..)
  mask <- apply(gmt, 1, function(x) { sum(!is.na(x)) }) > 2
  gmt <- gmt[mask, ]
  
  # convert to an n x 2 (annotation, gene) mapping
  res <- do.call('rbind', apply(gmt, 1, function(entry) {
    # gene id's starting in third column
    gids <- entry[3:length(entry)]
    gids <- gids[!is.na(gids)]

    data.frame(category = entry[1], gene = gids, row.names=NULL)
  }))

  # convert entrez gene ids to numer
  if (startsWith(keytype, 'entrez')) {
    res$gene <- as.numeric(as.character(res$gene))
  }

  # convert annotation column to character
  res$category <- as.character(res$category)

  # return dataframe
  res
}
