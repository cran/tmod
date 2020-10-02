## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  collapse = TRUE,
  comment = "#>"
)
library(dplyr)

## ----setup--------------------------------------------------------------------
library(tmod)

## ----preparation1-------------------------------------------------------------
data(vaccination)
vaccination %>% 
  select(GeneName, Description, logFC.F.D1, qval.F.D1) %>%
  arrange(qval.F.D1)

## ----tmod_results1, eval=FALSE------------------------------------------------
#  days <- paste0("D", c(1:5, 7))
#  res <- lapply(days, function(d) {
#    q <- paste0("qval.F.", d)
#    l <- vaccination$GeneName[ order(vaccination[[q]]) ]
#    tmodCERNOtest(l, mset="LI", qval=.1)
#  })
#  
#  lfcs <- vaccination %>% select(starts_with("logFC"))
#  pvals <- vaccination %>% select(starts_with("qval"))
#  
#  pie <- tmodDecideTests(vaccination$GeneName, lfc=lfcs, pval=pvals)
#  names(pie) <- names(res) <- days

