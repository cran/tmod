## ----setup, cache=FALSE, echo=FALSE, results="hide", include=FALSE-------
library(knitr)
opts_chunk$set(cache=FALSE, autodep=FALSE, tidy=FALSE)
opts_knit$set(width=75)

## ----twoA, size="Huge"---------------------------------------------------
library(limma)
library(tmod)
data(Egambia)
design <- cbind(Intercept=rep(1, 30), TB=rep(c(0,1), each= 15))
E <- as.matrix(Egambia[,-c(1:3)])
fit <- eBayes( lmFit(E, design))
tt <- topTable(fit, coef=2, number=Inf, 
  genelist=Egambia[,1:3] )

head(tt, 10)

## ----basic2, fig.width=4, fig.height=4-----------------------------------
group <- rep( c("CTRL", "TB"), each=15)
showGene(E["20799",], group,
  main=Egambia["20799", "GENE_SYMBOL"])

## ----mod1----------------------------------------------------------------
fg <- tt$GENE_SYMBOL[tt$adj.P.Val < 0.05 & abs( tt$logFC ) > 1]
res <- tmodHGtest(fg=fg, bg=tt$GENE_SYMBOL)
res

## ----mod2----------------------------------------------------------------
l    <- tt$GENE_SYMBOL
res2 <- tmodUtest(l)
head(res2)
nrow(res2)

## ----fourB---------------------------------------------------------------
l    <- tt$GENE_SYMBOL
res2 <- tmodCERNOtest(l)
head( res2 )

## ----five,fig=TRUE,fig.width=7,fig.height=4------------------------------
evidencePlot(l, "LI.M75")

## ----limma1--------------------------------------------------------------
res.l <- tmodLimmaTest(fit, Egambia$GENE_SYMBOL)
length(res.l)
names(res.l)
head(res.l$TB)

## ----limma2,fig=TRUE,fig.width=8,fig.height=4----------------------------
plotCI <- function(x, ci.l, ci.r, title="") {
  n <- length(x)
  plot(x, 
    ylab="logFC", xlab="Index", 
    pch=19, ylim=c( min(x-ci.l), max(x+ci.r)),
    main=title)
  segments(1:n, ci.l, 1:n, ci.r, lwd=5, col="#33333333")
}

par(mfrow=c(1,3))

x <- tmodLimmaTopTable(fit, coef="TB")
print(head(x))
x <- x[ x$logFC.TB > 0, ] # only to simplify the output!
x2 <- x[ order(abs(x$logFC.TB), decreasing=T),][1:50,]
plotCI(x2$logFC.TB, x2$ciL.TB, x2$ciR.TB, "logFC")

x2 <- x[ order(x$qval.TB),][1:50,]
plotCI(x2$logFC.TB, x2$ciL.TB, x2$ciR.TB, "q-value")

x2 <- x[ order(x$msd.TB, decreasing=T),][1:50,]
plotCI(x2$logFC.TB, x2$ciL.TB, x2$ciR.TB, "MSD")

## ----limma3--------------------------------------------------------------
x <- tmodLimmaTopTable(fit, coef="TB", genelist=Egambia[,1:3])
x.lfc  <- x[ order(abs(x$logFC.TB), decreasing=T),]
x.qval <- x[ order(x$qval.TB),]
x.msd  <- x[ order(x$msd.TB, decreasing=T),]

head(tmodCERNOtest(x.lfc$GENE_SYMBOL))
head(tmodCERNOtest(x.qval$GENE_SYMBOL))
head(tmodCERNOtest(x.msd$GENE_SYMBOL))

## ----pplot0--------------------------------------------------------------
head(tmodSummary(res.l), 5)

## ----pplot1,fig=TRUE,fig.width=6,fig.height=8----------------------------
tmodPanelPlot(res.l, text.cex=0.8)

## ----pplot2--------------------------------------------------------------
pie <- tmodLimmaDecideTests(fit, genes=Egambia$GENE_SYMBOL)
head(pie$TB[ order( pie$TB[,"Up"], decreasing=T), ])
data(tmod)
tmod$MODULES["DC.M3.4",]

## ----pplot3,fig=TRUE,fig.width=5,fig.height=6----------------------------
tmodPanelPlot(res.l, pie=pie, text.cex=0.8)

## ----pplot4,fig=TRUE,fig.width=6,fig.height=6----------------------------
tmodPanelPlot(res.l, 
  pie=pie, pie.style="rug", 
  grid="between")

## ----pplot5--------------------------------------------------------------
tt.I <- 
  topTable(fit, coef="Intercept", number=Inf, sort.by="n")
tt.TB <- topTable(fit, coef="TB", number=Inf, sort.by="n")
pie2 <- tmodDecideTests(Egambia$GENE_SYMBOL,
  lfc=cbind(tt.I$logFC, tt.TB$logFC),
  pval=cbind(tt.I$adj.P.Val, tt.TB$adj.P.Val))
identical(pie[[1]], pie2[[1]])

## ----fiveB---------------------------------------------------------------
l    <- tt$GENE_SYMBOL
res2 <- tmodUtest(l, mset="all")
head( res2 )

## ----six,fig=TRUE,fig.width=8,fig.height=5-------------------------------
library(pca3d)
mypal <- c("#E69F00", "#56B4E9")
pca <- prcomp(t(E), scale.=TRUE)
par(mfrow=c(1, 2))
l<-pca2d(pca, group=group, 
  palette=mypal)
cols <- as.character(l$colors)
legend("topleft", as.character(l$groups),
       pch=l$shapes,
       col=cols, bty="n")
l<-pca2d(pca, group=group, components=3:4,
  palette=mypal)
legend("topleft", as.character(l$groups),
       pch=l$shapes,
       col=cols, bty="n")

## ----seven---------------------------------------------------------------
o <- order(abs(pca$rotation[,4]), decreasing=TRUE)
l <- Egambia$GENE_SYMBOL[o]
res <- tmodUtest(l)
head(res)

## ----pcasum--------------------------------------------------------------
# Calculate enrichment for each component
gs   <- Egambia$GENE_SYMBOL
# function calculating the enrichment of a PC
gn.f <- function(r) {
    tmodCERNOtest(gs[order(abs(r), decreasing=T)],
                qval=0.01)
}
x <- apply(pca$rotation, 2, gn.f)
tmodSummary(x, filter.empty=TRUE)[1:5,]

## ----pcsum2,fig=TRUE,fig.width=10,fig.height=8,results="hide"------------
tmodPanelPlot(x)

## ----pcsum3,fig=TRUE,fig.width=10,fig.height=8---------------------------
qfnc <- function(r) quantile(r, 0.75)
qqs <- apply(pca$rotation[,1:10], 2, qfnc)
pie <- tmodDecideTests(gs, lfc=pca$rotation[,1:10], lfc.thr=qqs)
tmodPanelPlot(x[1:10], pie=pie, 
  pie.style="rug", grid="between")

## ----eight,fig=TRUE,fig.width=6,fig.height=5,results="hide"--------------
library(tagcloud)
w <- -log10(res$P.Value)
c <- smoothPalette(res$AUC, min=0.5)
tags <- strmultline(res$Title)
tagcloud(tags, weights=w, col=c)

## ----perm1,fig=FALSE-----------------------------------------------------
permset <- function(data, design) {
  require(limma)
  data <- data[, sample(1:ncol(data)) ]
  fit  <- eBayes(lmFit(data, design))
  tt   <- topTable(fit, coef=2, number=Inf, sort.by="n")
  order(tt$P.Value)
}

## ----perm2,fig=FALSE-----------------------------------------------------
# same design as before
design <- cbind(Intercept=rep(1, 30), 
  TB=rep(c(0,1), each= 15))
E      <- as.matrix(Egambia[,-c(1:3)])
N      <- 250  # small number for the sake of example
set.seed(54321)
perms  <- sapply(1:N, function(x) permset(E, design))
pauc   <- tmodAUC(Egambia$GENE_SYMBOL, perms)
dim(perms)

## ----perm3,fig=FALSE-----------------------------------------------------
fit <- eBayes(lmFit(E, design))
tt  <- topTable(fit, coef=2, number=Inf, 
  genelist=Egambia[,1:3])
res <- tmodCERNOtest(tt$GENE_SYMBOL, qval=Inf, order.by="n")
all(res$ID == rownames(perms))
fnsum <- function(m) sum(pauc[m,] >= res[m,"AUC"])
sums <- sapply(res$ID, fnsum)
res$perm.P.Val <- sums / N
res$perm.P.Val.adj <- p.adjust(res$perm.P.Val)
res <- res[order(res$AUC, decreasing=T),]
head(res[order(res$perm.P.Val),
  c("ID", "Title", "AUC", "adj.P.Val", "perm.P.Val.adj") ])

## ----ten-----------------------------------------------------------------
data(tmod)
res <- tmodUtest(tt$GENE_SYMBOL, qval=Inf)
gstest <- function(x) {
  sel <- tt$GENE_SYMBOL %in% tmod$MODULES2GENES[[x]]
  geneSetTest(sel, tt$logFC)
}
gst <- sapply(res$ID, gstest)

## ----eleven,fig=TRUE,fig.width=6,fig.height=6----------------------------
plot(res$P.Value, gst, 
  log="xy", pch=19,
  col="#33333366",
  xlab="P Values from tmod", 
  ylab="P Values from geneSetTest")
abline(0,1)
abline(h=0.01, col="grey")
abline(v=0.01, col="grey")

## ----exmset--------------------------------------------------------------
mymset <- makeTmod(
  modules=data.frame(ID=c("A", "B"),
             Title=c("A title", 
                      "B title")),
  modules2genes=list(
    A=c("G1", "G2"),
    B=c("G3", "G4"))
)
mymset

## ----msigdbparse1,eval=FALSE---------------------------------------------
#  msig <- tmodImportMSigDB("msigdb_v5.0.xml")
#  msig

## ----msigdb6,eval=FALSE--------------------------------------------------
#  res <- tmodCERNOtest(tt$GENE_SYMBOL, mset=msig )
#  head(res)

## ----msigdb7,eval=FALSE--------------------------------------------------
#  sel <- msig$MODULES$Category == "H"
#  tmodCERNOtest(tt$GENE_SYMBOL, mset=msig[sel] )

## ----msigdb1,eval=FALSE--------------------------------------------------
#  library(XML)
#  foo  <- xmlParse( "/home/january/Projects/R/pulemodule/vignette/msigdb_v5.0.xml" )
#  foo2 <- xmlToList(foo)

## ----msigdb2,eval=FALSE--------------------------------------------------
#  path1 <- foo2[[1]]
#  class(path1)

## ----msigdb2b,eval=FALSE-------------------------------------------------
#  names(path1)

## ----msigdb3,eval=FALSE--------------------------------------------------
#  orgs <- sapply(foo2, function(x) x["ORGANISM"])
#  unique(orgs)
#  
#  foo3 <- foo2[ orgs == "Homo sapiens" ]
#  foo3 <- foo3[ ! sapply(foo3, is.null) ]

## ----msigdb4,eval=FALSE--------------------------------------------------
#  modules <- t(sapply(foo3,
#    function(x)
#      x[ c("SYSTEMATIC_NAME", "STANDARD_NAME", "CATEGORY_CODE", "SUBCATEGORY_CODE") ]))
#  colnames(modules) <- c( "ID", "Title", "Category", "Subcategory" )
#  modules <- data.frame(modules, stringsAsFactors=FALSE)

## ----msigdb5,eval=FALSE--------------------------------------------------
#  m2g <- lapply(foo3,
#    function(x) strsplit( x["MEMBERS_SYMBOLIZED"], "," )[[1]])
#  names(m2g) <- modules$ID
#  
#  msig <- makeTmod(modules=modules, modules2genes=m2g)

## ----wpex3---------------------------------------------------------------
pathways <- data.frame( ID=files$ID, 
  Title=files$Title,
  stringsAsFactors=FALSE )
pathways$N <- sapply(p2GENES, length)
pathways$URL <- 
  paste0("http://www.wikipathways.org/index.php/Pathway:",
  pathways$ID )
sel <- pathways$N > 4
pathways <- pathways[ sel, ]
rownames(pathways) <- pathways$ID

## ----wpex4---------------------------------------------------------------
Hspaths <- makeTmod(modules=pathways, modules2genes=p2GENES)

## ----wpex5---------------------------------------------------------------
tmodCERNOtest(tt$GENE_SYMBOL, mset=Hspaths)

## ----wpex6---------------------------------------------------------------
sel <- grep( "Interferon", 
  Hspaths$MODULES$Title, ignore.case=T )
tmodCERNOtest(tt$GENE_SYMBOL, mset=Hspaths, 
  modules=Hspaths$MODULES$ID[sel]) 

