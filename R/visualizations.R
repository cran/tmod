
## adds an alpha chanel to a color
.alphacol <- function(x, suff="99") {
  x <- col2rgb(x)
  x <- apply(x, 2, function(xx) paste0("#", paste(as.hexmode(xx), collapse= ""), suff) )
  x
}


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
#' fg <- sample(c(tmod$MODULES2GENES[["LI.M127"]], bg[1:1000]))
#' hgEnrichmentPlot(fg, bg, "LI.M127")
#' @export
hgEnrichmentPlot <- function(fg, bg, m, mset="all", ...) {

  mset <- .getmodules2(NULL, mset)

  if(is.null(mset$MODULES2GENES[[m]])) stop("No such module")

  fg <- unique(fg)
  bg <- unique(c(fg, bg))
  b <- sum(fg %in% mset$MODULES2GENES[[m]])
  n <- length(fg)
  B <- sum(bg %in% mset$MODULES2GENES[[m]])
  N <- length(bg)

  mm <- matrix(c(b/n, 1-b/n, B/N, 1-B/N), byrow= F, nrow= 2)

  barplot(mm,
  xlim=c(0, 5), legend.text=c("In module", "Other"), bty="n", names.arg= c(
  "Foreground", "Background"), col= c("#E69F0099",  "#56B4E999"), ...)

    
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
#' @param rug if TRUE, draw a rug-plot beneath the ROC curve
#' @param roc if TRUE, draw a ROC curve above the rug-plot
#' @param filter if TRUE, genes not defined in the module set will be removed
#' @param unique if TRUE, duplicates will be removed
#' @param add if TRUE, the plot will be added to the existing plot
#' @param gene.labels a named character vector with gene labels to be shown (only if rug is plotted)
#' @param gl.cex Text cex (magnification) for gene labels
#' @param style "roc" for receiver-operator characteristic curve (default), and "gsea" for GSEA-style (Kaplan-Meier like plot)
#' @param col a character vector color to be used
#' @param col.rug a character value specifying the color of the rug
#' @param lwd line width (see par())
#' @param lty line type (see par())
#' @param rug.size fraction of the plot that should show the rug. If rug.size is 0, rug is not drawn. If rug.size is 1, ROC curve is not drawn.
#' @param legend position of the legend. If NULL, no legend will be drawn
#' @param ... Further parameters passed to the plotting function
#' @seealso \code{\link{tmod-package}}, \code{\link{hgEnrichmentPlot}}
#' @examples 
#' # artificially enriched list of genes
#' set.seed(123)
#' data(tmod)
#' bg <- tmod$GENES$ID
#' fg <- sample(c(tmod$MODULES2GENES[["LI.M127"]], bg[1:1000]))
#' l <- unique(c(fg, bg))
#' evidencePlot(l, "LI.M127")
#' evidencePlot(l, filter=TRUE, "LI.M127")
#' @export
evidencePlot <- function(l, m, mset="all", scaled= TRUE, rug=TRUE, roc=TRUE,
  filter= FALSE, unique=TRUE, add= FALSE, col="black", 
  col.rug="#eeeeee",
  gene.labels=NULL, gl.cex=1,
  style="roc",
  lwd=1, lty=1, rug.size=0.2,
  legend=NULL, ...) {

  if(rug.size == 0) rug <- FALSE
  if(rug.size == 1) roc <- FALSE
  if(!rug && !roc) stop("Both rug and roc are FALSE, nothing to plot")
  if(!roc) rug.size <- 1
  if(!rug) rug.size <- 0

  style <- match.arg(style, c("roc", "gsea"))

  # standardize the modules
  m <- as.character(m)
  mset <- .getmodules2(NULL, mset)

  if(! all(m %in% names(mset$MODULES2GENES))) stop("No such module")
  if(filter) l <- l[ l %in% mset$GENES$ID ]
  if(unique) l <- unique(l)

  n <- length(l)
  Nm <- length(m)

  # check graphical parameters and fill up if necessary
  col <- col[(((0:(Nm-1)) ) %% length(col))+1]
  lty <- lty[(((0:(Nm-1)) ) %% length(lty))+1]
  lwd <- lwd[(((0:(Nm-1)) ) %% length(lwd))+1]

  # check gene labels
  if(!is.null(gene.labels)) {
    gene.labels <- unique(as.character(gene.labels))
    
    if(is.null(names(gene.labels))) {
      names(gene.labels) <- gene.labels
    }

  }

  gene.labels <- gene.labels[ names(gene.labels) %in% l ]
  gene.labels <- gene.labels[ order(match(names(gene.labels), l)) ]

  # for each module to draw, find out which genes from l are found in that
  # module
  x <- sapply(m, function(mm) l %in% mset$MODULES2GENES[[mm]])



  # cumulative sum or scaled cumulative sum
  if(roc) {
    if(scaled) {
      cfunc <- function(xx) cumsum(xx) / sum(xx)
      r <- c(0, 1)
      ylab="Fraction of genes in module"
    } else {
      cfunc <- function(xx) cumsum(xx)
      r <- range(xcs)
      ylab <- "Number of genes in module"
    }

    cfunc2 <- switch(style,
      roc=cfunc,
      gsea=function(xx) { cfunc(xx) - cfunc(!xx) }
      )

    xcs <- apply(x, 2, cfunc2)

    # additional space for the rug
    r[1] <- - rug.size * (r[2]-r[1]) / (1-rug.size)
  } else {
    ylab <- ""
    r <- c(-1, 0)
  }

  if(!add) {
    plot(NULL, xlim= c(1, n), ylim= r, bty="n", xlab="List of genes", ylab=ylab, yaxt="n", ...)
    if(roc) axis(side=2, at= axisTicks(range(xcs), log=FALSE))
  }

  if(roc) {
    # plot the actual ROC curves
    for(i in 1:Nm) { lines(1:n, xcs[,i], col= col[i], lty=lty[i], lwd=lwd[i]) }

    # draw the diagonal
    segments(0, 0, n, r[2], col= "grey")
  }

  # plot the rug(s)
  if(rug) { 
    step <- r[1] / Nm

    # draw the rug background
    rect(0, (1:Nm) * step - step * 0.8, n, (1:Nm) * step, col= "#eeeeee", border= NA) 

    # draw each individual rug
    for(i in 1:Nm) {
      w <- which(x[,i])
      segments(w, step * i, w, step * (i-0.8), col= col[i], lty=lty[i], lwd=lwd[i])
    }

    # add gene labels
    if(!is.null(gene.labels)) {
      lw <- max(strheight(gene.labels, units="i", cex=gl.cex))
      lw <- grconvertX( c(0, lw), from="i", to="u")
      lw <- lw[2] - lw[1]
      lpos <- lw * 1.5 * (1:length(gene.labels))
      rpos <- match(names(gene.labels), l)

      for(i in 1:length(lpos)) {
        if(rpos[i] + lw > lpos[i]) {
          sel <- i:length(lpos)
          lpos[sel] <- lpos[sel] + (rpos[i] - lpos[i] + lw)
        }
      }

      text(lpos, 0.1, gene.labels, srt=90, pos=4)

      segments(rpos, step * 0.2, lpos + lw/2, 0.08)
    }
  }

  # add a legend
  if(!is.null(legend)) {
    legend(legend, m, col=col, lty=lty, lwd=lwd, bty= "n")
  }
}

