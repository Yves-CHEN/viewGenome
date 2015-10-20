\name{viewGenome-package}
\alias{viewGenome-package}
\alias{viewGenome}
\alias{myplot}
\alias{plotPloidy}
\docType{package}
\title{
View human genome in cytobands
}
% \details{
% \tabular{ll}{
% Package: \tab viewGenome\cr
% Type: \tab Package\cr
% Version: \tab 1.0\cr
% Date: \tab 2015-07-18\cr
% License: \tab GPL-2
% }
% }
\description{
    This is designed to view human genome with a chromsomal ideogram (HG19).
    The plot function has a "chr" parameter,
    which specify the specific chromosome to plot. \cr
    e.g.  \link{https://www.pinterest.com/pin/475411304393375351/} \cr
    The plots includes    
    \itemize{
     \item{ the regular 2D plot with a cytoband in the x-axis,    }
     \item{ A digital SKY-Plot. This provides a view of the genome, with copy number alterations.}
    }
}

\usage{

\method{myplot}{default} (x,y, chr, chrWidth = 0.3, col = 1, graybg = F, checks = c(1,10000000),...)
\method{plotPloidy}{default} (chr, CN, col=2, lwd=0.4, chrWidth=0.6, maxCN=12, ...)
\method{genomeInSubplots}{default} (x, y, grouping)
}
\arguments{
    \item{x}{ A numeric vector for the location. }
    \item{y}{ A numeric vector for copy numbers (or other informations, e.g. depth). }
    \item{chrWidth}{ Indicate the width of the cytobands plotted at bottom. }
    \item{col}{ A single numeric number of color code. }
    \item{...}{ These parameters will be passed to plot function. }
    \item{CN}{CN is a data.frame, with 4 colums of startPos, endPos, copy number, minor allele copy number. }
    \item{grouping}{grouping should have length as x. 
        It specifies which of the x belonging to which group. }
    \item{pch}{ specify the shape of the brush for painting. }
    \item{type}{ type of the dots. e.g. "l" for line, "p" for point. see par}
}

\author{
    Wenhan Chen
}
%\references{
%}
\keyword{d.SKY, chromosomal ideograms, cytobands}
%\seealso{
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
    library(viewGenome)
    runtest()
}
