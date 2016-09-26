#' Measures correlation between each principle component (PC) and each data
#' feature of interest (e.g. condition, batch, etc.)
#'
#' This is based on code adapted from cbcbSEQ
#' (https://github.com/kokrah/cbcbSEQ/) originally written by Kwame Okrah.
#'
#' @author Kwame Okrah
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param dat Matrix of input data
#' @param features Data frame of one of more features associated with each
#' column in the input data matrix.
#'
###############################################################################
pca_feature_correlations <- function(dat, features) {
    dat <- as.matrix(dat)
    features <- as.data.frame(features)

    # Drop any covariates with only a single level
    single_level <- apply(features, 2, function(x) {length(table(x))}) == 1

    features <- features[,!single_level]

    # SVD
    s <- corpcor:::fast.svd(dat - rowMeans(dat))
    rownames(s$v) <- colnames(dat)

    # Create output dataframe
    pc_var <- round((s$d^2) / sum(s$d^2) * 100, 2)

    result <- data.frame(
        pc_var=pc_var,
        pc_var_cum=cumsum(pc_var) 
    )

    # Measure feature correlations
    for (i in 1:ncol(features)) {
        feature_cor <- function(y) {
            round(summary(lm(y~features[,i]))$r.squared*100, 2)
        }
        result <- cbind(result, apply(s$v, 2, feature_cor))
    }
    colnames(result) = c('pc_var', 'pc_var_cum', colnames(features))

    return(result)
}
