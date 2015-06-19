## This function implements a friendly palette to be used in graphics
.mypalette <- function (n = NULL, transparent = "99") {
    pal <- "E69F00 56B4E9 009E73 F0E442 0072B2 D55E00 CC79A7 669900 660099 996600 990066 666633 666600 aa3366 5B4E85 FF6A5C ADAEA3 A0A376 FF8040 A2D6DA DA9CA5"

    pal <- unlist(strsplit(pal, " "))
    pal <- paste("#", pal, transparent, sep = "")
    if (!is.null(n)) {
        if (n > length(pal)) {
            pal <- rep(pal, ceiling(n/length(pal)))
        } else {
            pal <- pal[1:n]
        }
    }
    return(pal)
}

#' Select genes belonging to a module from a data frame
#'
#' Select genes belonging to a module from a data frame
#' 
#' showModule filters a data frame such that only genes from a module are
#' shown.
#' @param df a data frame
#' @param genes a character vector with gene IDs
#' @param module a single character value, ID of the module to be shown
#' @param mset Module set to use; see "tmodUtest" for details
#' @export
showModule <- function(df, genes, module, mset="all") {
  mset <- .getmodules2(NULL, mset)

  if(!is(df, "data.frame")) stop( "df must be a data frame" )
  if(!module %in% mset$MODULES$ID) stop("no such module")

  sel <- genes %in% mset$MODULES2GENES[[module]]

  if(sum(sel) == 0) warning("No genes belonging to that module found")

  df[ genes %in% mset$MODULES2GENES[[module]], ]
}


#' A combined beeswarm / boxplot
#'
#' A combined beeswarm / boxplot
#'
#' This is just a simple wrapper around the beeswarm() and boxplot()
#' commands.
#' @param data a vector of numeric values to be plotted
#' @param group factor describing the groups
#' @param main title of the plot
#' @param pch character to plot the points
#' @param xlab,ylab x and y axis labels
#' @param las see par()
#' @param pwcol colors of the points (see beeswarm)
#' @param ... any additional parameters to be passed to the beeswarm command
#' @examples
#' data(Egambia)
#' E <- as.matrix(Egambia[,-c(1:3)])
#' showGene(E["20799",], rep(c("CTRL", "TB"), each=15))
#' @importFrom beeswarm beeswarm
#' @export
showGene <- function( data, group, main= "", pch= 19, 
                         xlab= "", ylab= "log2 expression", las= 2, pwcol= NULL, ... ) {
  group  <- factor( group )
  pal    <- .mypalette( n= length( unique( group ) ) )

  if(! is.null(pwcol) & length( pwcol ) == 1 ) pwcol <- rep( pwcol, length( group ) )
  if(  is.null(pwcol)) pwcol= pal[ group ]

  
  beeswarm( data ~ group, 
    pch= pch, xlab= xlab, ylab= ylab, main= main, las= las, 
    pwcol= pwcol, bty="n",
    ... )

  boxplot( data ~ group, col= "#ffffff00", add= T, yaxt= "n", xaxt= "n", main= "", outline= FALSE, frame= FALSE )

  return( invisible( list( groups= levels( group ), col= pal ) ) )
}


#' Get genes belonging to a module
#'
#' Get genes belonging to a module
#'
#' Create a data frame mapping each module to a comma separated list of
#' genes. If genelist is provided, then only genes in that list will be
#' shown. An optional column, "fg" informs which genes are in the "foreground"
#' data set.
#' @return data frame containing module to gene mapping
#' @param modules module IDs 
#' @param genelist list of genes 
#' @param mset module set to use
#' @param fg genes which are in the foreground set
#' @export
getGenes <- function(modules, genelist=NULL, fg=NULL, mset="LI") {
  mset <- .getmodules2(modules, mset)

  if(!is.null(genelist)) 
    mset$MODULES2GENES <- lapply(mset$MODULES2GENES, function(x) x[ x %in% genelist ])

  ret <- data.frame(ID=mset$MODULES$ID)
  rownames(ret) <- ret$ID
  ret$N <- sapply(mset$MODULES2GENES, length)
  ret$Genes <- sapply(mset$MODULES2GENES, function(x) paste(x, collapse=","))

  if(!is.null(fg)) {
    ret$fg <- sapply(mset$MODULES2GENES, function(x) paste(x[x %in% fg], collapse=","))
  }
  ret
}

