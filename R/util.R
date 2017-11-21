#' Generate a color palette using ggplot defaults
#'  http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#'
#' @param n number of colors to include in palette
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length=n + 1)
    hcl(h=hues, l=65, c=100)[1:n]
}

#' Outputs an HTML datatable with all rows or static pandoc-friendly table with
#' a limited set of rows depending on the current output target.
#'
#' This function requires that the document requires that the knitr script is
#' being generated using the rmarkdown `render` function.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param dat Data to display
#' @param nrows Maximum number of rows to show per page for HTML output
#' @param caption Caption to use for table
#' @param digits Number of digits to disply for numeric
#' @param str_max_width Maximum width for string columns
#' @param output_format Desired output format [html|latex]
#' @param datatable_theme Theme to use when rendering the table with datatable
#'
#' @return None
xkable <- function (dat, nrows=10, caption=NULL, digits=getOption("digits"), 
                    str_max_width=Inf, output_format=NA,
                    datatable_style='bootstrap', ...) {
    # Conver to data.frame (necessary for handling mixed datatypes properly)
    dat <- as.data.frame(dat)

    # Trim strings if desired
    if (str_max_width < Inf) {
        str_cols <- sapply(dat, is.character)

        # Iterate over character columns
        for (cname in colnames(dat)[str_cols]) {
            # For all entries that exceed the maximum width, trim
            exceeds_length <- nchar(dat[[cname]]) > str_max_width
            exceeds_length[is.na(exceeds_length)] <- FALSE
            dat[[cname]][exceeds_length] <- paste0(strtrim(dat[[cname]][exceeds_length], str_max_width - 3), '...')
        }
    }

    # Fix types for numeric columns
    numeric_cols <- sapply(dat, Hmisc::all.is.numeric)

    if (nrow(dat) > 0) {
        if (sum(numeric_cols) > 1)  {
            dat[,numeric_cols] <- round(data.matrix(dat[,numeric_cols]), digits)
        } else {
            dat[,numeric_cols] <- round(as.numeric(dat[,numeric_cols]), digits)
        }
    }

    # Determine output format
    if (is.na(output_format)) {
        if (is.null(opts_knit$get("rmarkdown.pandoc.to"))) {
            # Default to latex output (kable)
            output_format <- 'latex'
        } else {
            output_format <- opts_knit$get("rmarkdown.pandoc.to") 
        }
    }

    # HTML output
    if (nrow(dat) > nrows && output_format == 'html') {
        options(DT.options=list(pageLength=nrows))
        # 2017/04/24: printr package not working as expected, reverting to
        # calling datatable directly
        #kable(dat, 'html', caption=caption)

        # 2017/08/04 datatable not rendering to HTML and instead opening in
        # browser
        #datatable(dat, caption=caption)
        #renderDataTable(dat, caption=caption)
        kable(dat, caption=caption, ...)
    } else {
        # Static output
        kable(dat, caption=caption, ...)
    }
}
