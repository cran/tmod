## ----setup, echo=FALSE, results="hide", include=FALSE-------------------------
library(knitr)
opts_chunk$set(cache=TRUE, autodep=TRUE, tidy=FALSE, fig.width=5, warning=FALSE, fig.height=5, width=60)
opts_knit$set(width=60)
is.internet <- FALSE

## ----eval=TRUE,echo=FALSE-----------------------------------------------------
library(ggplot2)
theme_set(theme_minimal())

## ----twoA,eval=TRUE-----------------------------------------------------------
library(tmod)
data(Egambia)
E <- as.matrix(Egambia[,-c(1:3)])

## ----limma,eval=FALSE---------------------------------------------------------
#  library(limma)
#  design <- cbind(Intercept=rep(1, 30), TB=rep(c(0,1), each= 15))
#  fit <- eBayes(lmFit(E, design))
#  tt <- topTable(fit, coef=2, number=Inf,
#    genelist=Egambia[,1:3])

## ----twoA.2-------------------------------------------------------------------
data(EgambiaResults)
tt <- EgambiaResults

## ----twoA2,echo=FALSE, results="asis"-----------------------------------------
library(pander)
options(digits=3)
tmp <- tt[,c("GENE_SYMBOL", "GENE_NAME", "logFC", "adj.P.Val")]
rownames(tmp) <- NULL
pandoc.table(head(tmp), split.tables=Inf, justify="llrr")

## ----basic2, fig.width=4, fig.height=4----------------------------------------
group <- rep( c("CTRL", "TB"), each=15)
showGene(E["20799",], group,
  main=Egambia["20799", "GENE_SYMBOL"])

## ----fourB--------------------------------------------------------------------
l    <- tt$GENE_SYMBOL
resC <- tmodCERNOtest(l)
head(resC, 15)

## ----panelplot00,fig.width=14,fig.height=7------------------------------------
library(cowplot)

g1 <- ggPanelplot(list(Gambia=resC))

## calculate the number of significant genes
## per module
sgenes <- tmodDecideTests(g=tt$GENE_SYMBOL,
  lfc=tt$logFC,
  pval=tt$adj.P.Val)
names(sgenes) <- "Gambia"
g2 <- ggPanelplot(list(Gambia=resC), sgenes = sgenes)
plot_grid(g1, g2, labels=c("A", "B"))

## ----mod1-------------------------------------------------
fg <- tt$GENE_SYMBOL[tt$adj.P.Val < 0.05 & abs( tt$logFC ) > 1]
resHG <- tmodHGtest(fg=fg, bg=tt$GENE_SYMBOL)
options(width=60)
resHG

## ----mod2-------------------------------------------------
l    <- tt$GENE_SYMBOL
resU <- tmodUtest(l)
head(resU)
nrow(resU)

## ----mod2cerno--------------------------------------------
l    <- tt$GENE_SYMBOL
resCERNO <- tmodCERNOtest(l)
head(resCERNO)
nrow(resCERNO)

## ----plage------------------------------------------------
tmodPLAGEtest(Egambia$GENE_SYMBOL, Egambia[,-c(1:3)], group=group)

## ----heatmap0---------------------------------------------
m <- "LI.M75"
## getModuleMembers returns a list – you can choose to 
## select multiple modules
genes <- getModuleMembers(m)[[1]]
sel <- Egambia$GENE_SYMBOL %in% genes
x <- data.matrix(Egambia)[sel, -c(1:3)] # expression matrix

## ----heatmap1,fig.width=12,fig.height=7,echo=FALSE--------
x0 <- t(scale(t(x)))

cols <- colorRampPalette(c("blue", "white", "red"))(17) 
library(gplots)
heatmap.2(x0, trace="n", 
  labRow=Egambia$GENE_SYMBOL[sel], 
  scale="n", 
  dendrogram="r", Colv=F, 
  col=cols,
  breaks=seq(-2.5, 2.5, length.out=18))

## ----heatmap2,fig.width=10,fig.height=7,echo=FALSE--------
## per group means and standard deviations
ms <- apply(x, 1, function(xx) tapply(xx, group, mean))
sd <- apply(x, 1, function(xx) tapply(xx, group, sd))

library(plotwidgets)
plot(NULL, bty="n", ylim=range(x), xlim=c(0.5, 2.5), 
  xaxt="n", xlab="group", ylab="Expression")
axis(1, at=c(1,2), labels=unique(group))
errbar <- function(x, y0, y1, w=0.1, ...) {
  w <- w/2
  segments(x, y0, x, y1, ...)
  segments(x-w, y0, x+w, y0, ...)
  segments(x-w, y1, x+w, y1, ...)
}

cols <- plotPals("default")
n <- ncol(ms)
segments(1, ms[1,], 2, ms[2,], col=cols, lwd=2)
points(rep(1, n), ms[1,], pch=19, col=cols)
points(rep(2, n), ms[2,], pch=19, col=cols)
errbar(1, ms[1,]-sd[1,], ms[1,]+sd[1,], col=cols)
errbar(2, ms[2,]-sd[2,], ms[2,]+sd[2,], col=cols)

## ----eigengene,fig.width=8,fig.height=5-------------------
par(mfrow=c(1,2))
eig <- eigengene(Egambia[,-c(1:3)], Egambia$GENE_SYMBOL)
showGene(eig["LI.M75", ], group, 
  ylab="Eigengene",
  main="antiviral Interferon signature")
showGene(eig["LI.M16", ], group, 
  ylab="Eigengene",
  main="TLR and inflammatory signaling")

## ----five,fig=TRUE,fig.width=7,fig.height=4---------------
l    <- tt$GENE_SYMBOL
theme_set(theme_minimal())
ggEvidencePlot(l, "LI.M75") 

## ----fig.width=12-----------------------------------------
library(purrr)
sel <- c("LI.M67", "LI.M37.0")
plots <- map(sel, ~ ggEvidencePlot(l, .x, gene.labels=FALSE))
plot_grid(plotlist = plots, labels=sel)

## ---------------------------------------------------------
foo <- tmodCERNOtest(l) %>% dplyr::filter(ID %in% sel) 
foo %>% knitr::kable()

## ----fourC,size="tiny"------------------------------------
resAll <- list(CERNO=resC, U=resU, HG=resHG)
#head(tmodSummary(resAll))

## ----fourC2,results="asis",echo=FALSE,size="tiny"---------
tmp <- tmodSummary(resAll)
rownames(tmp) <- NULL
pandoc.table(head(tmp), split.tables=Inf, justify="llrrrrrr")

## ----panelplots, fig.width=8, fig.height=6----------------
resAll$HG$AUC <- log10(resAll$HG$E) - 0.5
ggPanelplot(resAll)

## ----panelplots2, fig.width=8, fig.height=5---------------
ggPanelplot(resAll, q_thr=1e-3)

## ----pie, fig.width=10,fig.height=6-----------------------
degs <- tmodDecideTests(g=tt$GENE_SYMBOL, lfc=tt$logFC,
                       pval=tt$adj.P.Val)[[1]]
degs <- list(CERNO=degs, HG=degs, U=degs)
ggPanelplot(resAll, sgenes = degs)

## ----ankrd22----------------------------------------------
x <- E[ match("ANKRD22", Egambia$GENE_SYMBOL), ]
cors <- t(cor(x, t(E)))
ord <- order(abs(cors), decreasing=TRUE)
head(tmodCERNOtest(Egambia$GENE_SYMBOL[ ord ]))

## ----m75--------------------------------------------------
g <- getGenes("LI.M75", genes=Egambia$GENE_SYMBOL, 
              as.list=TRUE)
x <- E[ match(g[[1]], Egambia$GENE_SYMBOL), ]

## calculating the "eigengene" (PC1)
pca <- prcomp(t(x), scale.=T)
eigen <- pca$x[,1]
cors <- t(cor(eigen, t(E)))

## order all genes by the correlation between the gene and the PC1
ord <- order(abs(cors), decreasing=TRUE)
head(tmodCERNOtest(Egambia$GENE_SYMBOL[ ord ]))

## ----six,fig=TRUE,fig.width=8,fig.height=5----------------
mypal <- c("#E69F00", "#56B4E9")
pca <- prcomp(t(E), scale.=TRUE)

col <- mypal[ factor(group) ]
par(mfrow=c(1, 2))
l<-pcaplot(pca, group=group, col=col)
 
legend("topleft", as.character(l$groups),
       pch=l$pch,
       col=l$colors, bty="n")
l<-pcaplot(pca, group=group, col=col, components=3:4)
legend("topleft", as.character(l$groups),
       pch=l$pch,
       col=l$colors, bty="n")

## ----seven------------------------------------------------
o <- order(abs(pca$rotation[,4]), decreasing=TRUE)
l <- Egambia$GENE_SYMBOL[o]
res <- tmodUtest(l)
head(res)

## ----pcasum-----------------------------------------------
# Calculate enrichment for each component
gs   <- Egambia$GENE_SYMBOL
# function calculating the enrichment of a PC
gn.f <- function(r) {
    tmodCERNOtest(gs[order(abs(r), decreasing=T)],
                qval=0.01)
}
x <- apply(pca$rotation, 2, gn.f)
tmodSummary(x, filter.empty=TRUE)[1:5,]

## ----pcsum2,fig=TRUE,fig.width=8,fig.height=5,results="hide"----
ggPanelplot(x)

## ----pcsum3,fig=TRUE,fig.width=8,fig.height=5-------------
qfnc <- function(r) quantile(r, 0.75)
qqs <- apply(pca$rotation[,1:10], 2, qfnc)
gloadings <- tmodDecideTests(gs, lfc=pca$rotation[,1:10], lfc.thr=qqs)
ggPanelplot(x[1:10], sgenes = gloadings) 

## ----fiveB------------------------------------------------
l    <- tt$GENE_SYMBOL
res2 <- tmodUtest(l, mset="LI")
head( res2 )

## ----tmodobj----------------------------------------------
data(tmod)
tmod

## ----tmodobj2---------------------------------------------
length(tmod)
sel <- grep("Interferon", tmod$MODULES$Title, ignore.case=TRUE)
ifn <- tmod[sel]
ifn
length(ifn)

## ----exmset-----------------------------------------------
mymset <- makeTmodGS(
  gs=data.frame(ID=c("A", "B"),
             Title=c("A title", 
                      "B title")),
  gs2gene=list(
    A=c("G1", "G2"),
    B=c("G3", "G4"))
)
mymset

## ----eval=FALSE-------------------------------------------
#  library(msigdbr)
#  msig <- msigdbr()
#  msig <- makeTmodFromDataFrame(df=msig,
#    feature_col="gene_symbol",
#    module_col="gs_id", title_col="gs_name",
#    extra_module_cols=c("gs_cat", "gs_subcat", "gs_url",
#                        "gs_exact_source", "gs_description"))

## ----msigdbparse1, eval=FALSE-----------------------------
#  msig <- tmodImportMSigDB("msigdb_v2022.1.Hs.xml")

## ----msigdb6,eval=FALSE-----------------------------------
#  res <- tmodCERNOtest(tt$GENE_SYMBOL, mset=msig)

## ----msigdb7,eval=FALSE-----------------------------------
#  sel <- msig$gs$gs_cat == "H"
#  tmodCERNOtest(tt$GENE_SYMBOL, mset=msig[sel] )

## ----biomart,eval=FALSE-----------------------------------
#  library(biomaRt)
#  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#  bm <- getBM(filters="hgnc_symbol",
#              values = Egambia$GENE_SYMBOL,
#              attributes = c( "hgnc_symbol", "entrezgene", "reactome", "go_id", "name_1006", "go_linkage_type"),
#              mart=mart)

## ----biomart2,eval=FALSE----------------------------------
#  m2g_r <- with(bm[ bm$reactome != "", ], split(hgnc_symbol, reactome))
#  m2g_g <- with(bm[ bm$go_id != "", ], split(hgnc_symbol, go_id))
#  
#  ll <- lengths(m2g_r)
#  m2g_r <- m2g_r[ ll >= 5 & ll <= 250 ]
#  ll <- lengths(m2g_g)
#  m2g_g <- m2g_g[ ll >= 5 & ll <= 250 ]
#  
#  m_r <- data.frame(ID=names(m2g_r), Title=names(m2g_r))
#  m_g <- data.frame(ID=names(m2g_g),
#    Title=bm$name_1006[ match(names(m2g_g), bm$go_id)])
#  
#  ensemblR  <- makeTmod(modules=m_r, modules2genes=m2g_r)
#  ensemblGO <- makeTmod(modules=m_g, modules2genes=m2g_g)
#  
#  ## these objects are no longer necessary
#  rm(bm, m_g, m_r, m2g_r, m2g_g)

## ----orghs,eval=FALSE-------------------------------------
#  library(org.Hs.eg.db)
#  mtab <- toTable(org.Hs.egGO)

## ----orghs2,eval=FALSE------------------------------------
#  mtab <- mtab[ mtab$Ontology == "BP", ]
#  m2g <- split(mtab$gene_id, mtab$go_id)
#  ## remove the rather large object
#  rm(mtab)
#  ll <- lengths(m2g)
#  m2g <- m2g[ ll >= 10 & ll <= 100 ]
#  length(m2g)

## ----orghs3,eval=FALSE------------------------------------
#  library(GO.db)
#  gt <- toTable(GOTERM)
#  m <- data.frame(ID=names(m2g))
#  m$Title <- gt$Term[ match(m$ID, gt$go_id) ]
#  
#  goset <- makeTmod(modules=m, modules2genes=m2g)
#  rm(gt, m2g, m)

## ----kegg, eval=FALSE-------------------------------------
#  library(KEGGREST)
#  pathways <- keggLink("pathway", "hsa")
#  
#  ## get pathway Names in addition to IDs
#  paths    <- sapply(unique(pathways), function(p) keggGet(p)[[1]]$NAME)
#  m <- data.frame(ID=unique(pathways), Title=paths)
#  
#  ## m2g is the mapping from modules (pathways) to genes
#  m2g <- split(names(pathways), pathways)
#  
#  ## kegg object can now be used with tmod
#  kegg <- makeTmod(modules=m, modules2genes=m2g)

## ----kegg2, eval=FALSE------------------------------------
#  eg <- paste0("hsa:", tt$EG)
#  tmodCERNOtest(eg, mset="kegg")

## ----msigdb1,eval=FALSE-----------------------------------
#  library(XML)
#  foo  <- xmlParse( "msigdb_v2022.1.Hs.xml" )
#  foo2 <- xmlToList(foo)

## ----msigdb2,eval=FALSE-----------------------------------
#  path1 <- foo2[[1]]
#  class(path1)

## ----msigdb2b,eval=FALSE----------------------------------
#  names(path1)

## ----msigdb3,eval=FALSE-----------------------------------
#  orgs <- sapply(foo2, function(x) x["ORGANISM"])
#  unique(orgs)
#  
#  foo3 <- foo2[ orgs == "Homo sapiens" ]
#  foo3 <- foo3[ ! sapply(foo3, is.null) ]

## ----msigdb4,eval=FALSE-----------------------------------
#  modules <- t(sapply(foo3,
#    function(x)
#      x[ c("SYSTEMATIC_NAME", "STANDARD_NAME", "CATEGORY_CODE", "SUBCATEGORY_CODE") ]))
#  colnames(modules) <- c( "ID", "Title", "Category", "Subcategory" )
#  modules <- data.frame(modules, stringsAsFactors=FALSE)
#  nrow(modules)

## ----msigdb5,eval=FALSE-----------------------------------
#  m2g <- lapply(foo3,
#    function(x) strsplit( x["MEMBERS_SYMBOLIZED"], "," )[[1]])
#  names(m2g) <- modules$ID
#  
#  mymsig <- makeTmod(modules=modules, modules2genes=m2g)
#  mymsig

