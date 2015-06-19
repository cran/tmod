.alphacol <- function( x, suff="99" ) {
  x <- col2rgb( x )
  x <- apply( x, 2, function( xx ) paste0( "#", paste( as.hexmode( xx ), collapse= "" ), suff )  )
  x
}

#' Create a visualisation of enrichment
#'
#' Create a visualisation of enrichment
#' 
#' This functions creates a barplot visualizing the enrichment of a
#' module in the foreground (fg) set as compared to the background (bg) set.
#' It is the counterpart 
#' @param fg the foreground set of genes
#' @param bg the background set of genes (gene universe)
#' @param m module for which the plot should be created
#' @param mset Which module set to use (see tmodUtest for details)
#' @param ... additional parameters to be passed to the plotting function
#' @seealso \code{\link{tmod-package}}, \code{\link{evidencePlot}}
#' @examples 
#' set.seed(123)
#' data(tmod)
#' bg <- tmod$GENES$ID
#' fg <- sample( c(tmod$MODULES2GENES[["LI.M127"]], bg[1:1000]))
#' hgEnrichmentPlot(fg, bg, "LI.M127")
#' @export
hgEnrichmentPlot <- function( fg, bg, m, mset="all", ... ) {

  mset <- .getmodules2(NULL, mset)

  if(is.null(mset$MODULES2GENES[[m]])) stop("No such module")

  fg <- unique(fg)
  bg <- unique(c(fg, bg))
  b <- sum(fg %in% mset$MODULES2GENES[[m]])
  n <- length(fg)
  B <- sum(bg %in% mset$MODULES2GENES[[m]])
  N <- length(bg)

  mm <- matrix( c( b/n, 1-b/n, B/N, 1-B/N ), byrow= F, nrow= 2 )

  barplot( mm,
  xlim=c(0, 5), legend.text=c("In module", "Other"), bty="n", names.arg= c(
  "Foreground", "Background" ), col= c( "#E69F0099",  "#56B4E999" ), ... )

    
}



#' Create an evidence plot for a module
#'
#' Create an evidence plot for a module
#'
#' This function creates an evidence plot for a module, based on an
#' ordered list of genes. The plot shows the receiving operator
#' characteristic (ROC) curve and a rug below, which indicates the distribution of the
#' module genes in the sorted list.
#' @param l sorted list of HGNC gene identifiers
#' @param m character vector of modules for which the plot should be created
#' @param mset Which module set to use (see tmodUtest for details)
#' @param scaled if TRUE, the cumulative sums will be divided by the total sum (default)
#' @param filter if TRUE, genes not defined in the module set will be removed
#' @param add if TRUE, the plot will be added to the existing plot
#' @param legend position of the legend. If NULL, no legend will be drawn
#' @seealso \code{\link{tmod-package}}, \code{\link{hgEnrichmentPlot}}
#' @examples 
#' # artificially enriched list of genes
#' set.seed(123)
#' data(tmod)
#' bg <- tmod$GENES$ID
#' fg <- sample( c(tmod$MODULES2GENES[["LI.M127"]], bg[1:1000]))
#' l <- unique(c(fg, bg))
#' evidencePlot(l, "LI.M127")
#' evidencePlot(l, filter=tmod$GENES$ID, "LI.M127")
#' @export
evidencePlot <- function( l, m, mset="all", scaled= TRUE, filter= FALSE, add= FALSE, legend= "topleft" ) {

  m <- as.character( m )
  mset <- .getmodules2(NULL, mset)

  if( ! all( m %in% names( mset$MODULES2GENES ) ) ) stop( "No such module" )
  if( !is.null(filter) ) l <- l[ l %in% mset$GENES$ID ]

  n <- length( l )

  x <- sapply( m, function( mm ) l %in% mset$MODULES2GENES[[mm]] )

  if( scaled ) {
    xcs <- apply( x, 2, function( xx ) cumsum( xx ) / sum( xx ) )
  } else {
    xcs <- apply( x, 2, cumsum )
  }

  r <- range( xcs )
  rd <- 0.2 * (r[2]-r[1])
  r[1] <- r[1] - rd

  if( ! add ) {
    plot( NULL, xlim= c( 1, n ), ylim= r, bty="n", xlab="List of genes", ylab=m, yaxt="n" )
    rect( 0, r[1], n, 0, col= "#eeeeee", border= NA )
    axis(side=2, at= axisTicks( range(xcs), log=FALSE ) )
  }

  for( i in 1:length( m ) ) { lines( 1:n, xcs[,i], col= i ) }

  for( i in 1:length( m ) ) {
    w <- which( x[,i] )
    segments( w, r[1], w, 0, col= i )
  }

  segments( 0, 0, n, r[2], col= "grey" )

  if( ! is.null( legend ) ) {
    legend( legend, m, col= 1:length( m ), lty= 1, bty= "n" )
  }


}
