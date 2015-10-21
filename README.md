# viewGenome
viewGenome is a Rpackage with a collection of plots designed to view human genome by per chromosomes with cytobands, or view in d.SKY. This is designed to view human genome with a chromsomal ideogram (HG19). The plot function has a "chr" parameter, which specify the specific chromosome to plot
[view in pinit](https://www.pinterest.com/pin/475411304393375351/).
The plots includes    
* the regular 2D plot with a cytoband in the x-axis,    
* A digital SKY-Plot. This provides a view of the genome, with copy number alterations.

## Requirements
R >= (3.1.2), quantsmooth
## Install
R CMD INSTALL ./viewGenome
## Usage
    library(viewGenome)
    runtest()   # This generates a test.pdf file if the installation works through.
