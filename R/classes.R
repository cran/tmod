## we follow limma here

#' An S4 class that stores modules.
#' @rdname tmod-data
#' @importFrom methods setClass setMethod
#' @export
setClass( "tmod",
  representation("list")
)

## how to display a tmod object

#' Shows the tmod object
#'
#' @name show
#' @param object a tmod object
#' @aliases show,tmod-method
#' @rdname extract-methods
#' @docType methods
setMethod( "show", "tmod",
  function(object) {
    .catf( "An object of class \"%s\"\n", class(object) )
    .catf( "\t%d modules, %d genes\n", 
      nrow(object$MODULES),
      nrow(object$GENES) )
  })


## allow easy subsetting of tmod objects

#' Extracts parts of a tmod object
#'
#' @name [
#' @param x a tmod object
#' @param i indices specifying elements to extract or replace
#' @aliases [,tmod-method
#' @rdname extract-methods
#' @docType methods

setMethod("[", signature(x="tmod", i="ANY"),

  function(x, i) {

  object  <- x
  modules <- i
  modules <- object$MODULES[ modules, "ID" ]
  #modules <- modules[ modules %in% object$MODULES$ID ]

  MODULES       <- object$MODULES[modules,,drop=FALSE]
  MODULES2GENES <- object$MODULES2GENES[modules]

  glist <- unique(unlist(MODULES2GENES))

  GENES <- object$GENES[ object$GENES$ID %in% glist,,drop=FALSE]

  GENES2MODULES <- NULL
  if(!is.null(object$GENES2MODULES)) {
    GENES2MODULES <- object$GENES2MODULES[object$GENES$ID]
  }

  new( "tmod", list( 
    MODULES=MODULES,
    GENES=GENES,
    MODULES2GENES=MODULES2GENES,
    GENES2MODULES=GENES2MODULES ))

})
