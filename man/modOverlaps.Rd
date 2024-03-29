% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modstats.R
\name{modOverlaps}
\alias{modOverlaps}
\title{Calculate overlaps of the modules}
\usage{
modOverlaps(modules = NULL, mset = NULL, stat = "jaccard")
}
\arguments{
\item{modules}{either a character vector with module IDs from mset, or a list which
contains the module members}

\item{mset}{Which module set to use. Either a character vector ("LI", "DC" or "all", default: all) or an object of class tmod (see "Custom module definitions" below)}

\item{stat}{Type of statistics to return. 
"jaccard": Jaccard index (default);
"number": number of common genes;
"soerensen": Soerensen-Dice coefficient;
"overlap": Szymkiewicz-Simpson coefficient.}
}
\description{
Calculate overlaps of the modules
}
\details{
For a set of modules (aka gene sets) determine the similarity between
these. You can run modOverlaps either on a character vector of module
IDs or on a list. In the first case, the module elements are taken from
`mset`, and if that is NULL, from the default tmod module set. 

Alternatively, you can provide a list in which each element is a
character vector. In this case, the names of the list are the module IDs,
and the character vectors contain the associated elements.

The different statistics available are:
 * "number": total number of common genes (size of the overlap)
 * "jaccard": Jaccard index, i.e. \eqn{\frac{|A \cap B|}{|A \cup B|}}
   (number of common elements divided by the total number of unique elements);
 * "soerensen": Soerensen-Dice coefficient, defined as \eqn{\frac{2 \cdot |A \cap B|}{|A| + |B|}} – number of common genes in relation to the total number of elements (divided by two, such that the maximum is 1)
   (number of common elements divided by the average size of both gene sets)
 * "overlap": Szymkiewicz-Simpson coefficient, defined as \eqn{\frac{|A \cap B|}{\min(|A|, |B|)}} – this is the number of common genes scaled by the size of the smaller of the two gene sets
   (number of common elements divided by the size of the smaller gene set)
}
