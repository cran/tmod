## checking a new tmod object
check_tmod <- function(object) {

  if(is.null(object$MODULES))       return("Required list member MODULES missing")
  if(is.null(object$MODULES2GENES)) return("Required list member MODULES2GENES missing")

  if(!is(object$MODULES, "data.frame")) return("MODULES must be a data frame")
  if(!is(object$MODULES2GENES, "list")) return("MODULES2GENES must be a list")

  if(is.null(object$MODULES$ID))    return("MODULES must have columns ID and Title")
  if(is.null(object$MODULES$Title)) return("MODULES must have columns ID and Title")

  if(!all(object$MODULES$ID %in% names(object$MODULES2GENES)))
    return("All MODULES$ID must be found in names of MODULES2GENES")

  TRUE

}

#' @param modules A data frame with at least columns ID and Title
#' @param modules2genes A list with module IDs as names. Each member of the list is a character vector with IDs of genes contained in that module
#' @param genes2modules A list with gene IDs as names. Each member of the list is a character vector with IDs of modules in 
#'        which that gene is contained. This object will be created automatically if the provided parameter is NULL
#' @param genes A data frame with meta-information on genes. Must contain the column ID. If NULL, then a data frame with only one column (ID) will be created automatically.
#' @rdname tmod-class
#' @export
makeTmod <- function(modules, modules2genes, genes2modules=NULL, genes=NULL) {


  if(!is(modules, "data.frame")) stop("MODULES must be a data frame")
  if(!is(modules2genes, "list")) stop("MODULES2GENES must be a list")

  if(is.null(modules$ID))    stop("MODULES must have column ID")
  if(is.null(modules$Title)) stop("MODULES must have column Title")

  if(is.null(genes)) {
    genes <- unique(unlist(modules2genes))
    genes <- data.frame(ID=genes, stringsAsFactors=FALSE)
  } else {
    if(is.null(genes$ID)) {
      stop("genes must have an ID column")
    }
  }

  genes$ID     <- as.character(genes$ID)
  modules$ID   <- as.character(modules$ID)
  rownames(genes)   <- genes$ID
  rownames(modules) <- modules$ID

  if(!all(modules$ID %in% names(modules2genes)))
    stop("All MODULES$ID must be found in names of MODULES2GENES")


  if(is.null(genes2modules)) {
    genes2modules <- sapply(genes$ID, function(g) {
      sel <- sapply(modules2genes, function(m) g %in% m)
      names(modules2genes)[sel]
    }, simplify=F) 
  } else {
    if(!all(genes$ID %in% names(genes2modules)))
      stop("All genes$ID must be found in names of genes2modules")
  }


  mset <- new("tmod", list(MODULES=modules, MODULES2GENES=modules2genes, GENES2MODULES=genes2modules, GENES=genes))


}


#' S4 class for tmod
#'
#' S4 class for tmod
#'
#' An object of class tmod contains the module annotations (tmod$MODULES), gene to
#' module (tmod$GENES2MODULES) and module to gene (tmod$MODULES2GENES) mappings
#' and a gene list (tmod$GENES).
#'
#' tmod$MODULES and tmod$GENES are data frames, while tmod$MODULES2GENES and
#' tmod$GENES2MODULES are lists with, respectively, module and gene
#' identifiers as names. The data frames MODULES and GENES must contain the
#' column "ID", and the data frame MODULES must contain additionally the
#' column "Title".
#' 
#' Objects of class tmod should be constructed 
#' by calling the function makeTmod(). This function check the validity and
#' consistency of the provided objects, and, if necessary, automatically
#' creates the members GENES and GENES2MODULES.
#'
#' See the package vignette for more on constructing custom module sets.
#'
#' @rdname tmod-class
#' @importFrom methods setClass setMethod loadMethod is new
#' @seealso tmod-data
#' @examples
#' # A minimal example
#' m <- data.frame(ID=letters[1:3], Title=LETTERS[1:3])
#' m2g <- list(a=c("g1", "g2"), b=c("g3", "g4"), c=c("g1", "g2", "g4"))
#' mymset <- makeTmod(modules=m, modules2genes=m2g)
#' @export
setClass( "tmod",
  representation("list"),
  validity=check_tmod
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
