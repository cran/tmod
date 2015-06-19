#' Gene expression module data
#'
#' Gene expression module data
#'
#' This data set contains the module annotations (tmod$MODULES), gene to
#' module (tmod$GENES2MODULES) and module to gene (tmod$MODULES2GENES) mappings
#' and a gene list (tmod$GENES).
#'
#' tmod$MODULES and tmod$GENES are data frames, while tmod$MODULES2GENES and
#' tmod$GENES2MODULES are lists with, respectively, module and gene
#' identifiers as names.
#' @examples
#' # list of first 10 modules
#' data(tmod)
#' tmod
#' tmod$MODULES[1:10, ]
#' tmod[1:10]
#' @name tmod-data
NULL

#' @name tmod
#' @rdname tmod-data
NULL

#' @name tmod-class
#' @rdname tmod-data
NULL
