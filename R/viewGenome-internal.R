#library(quantsmooth)



getChrCol  <-  function (chr)
{

    data("skyColTheme")
    toCol  <- function(colCode) { rgb(colCode[1], colCode[2], colCode[3], maxColorValue = 255)}
    colForChr = toCol(colTab[as.character(chr), ])
    colForChr

}



#################################################################################
####  This adds a chromosomal ideogram to default plot, by specifiying which chromsomes
####          chr  : the id of an autosome, in between 1 and 22
####          chrWidth : specifiy the width of the ideogram.
####          
####  This function can also generate a ggplot2-alike gray checked background
####          graybg = T :  specify the gray background. Default value is F.
####          checks:  Only works when specify graybg = T. This specifiy how to 
####                   divide the gray backgroundi into squares. By default, 
####                   x axis (chromsomse) is by every 1M, y axis (Copy number) is by 1.
#################################################################################
myplot  <- function(x,y, chr, chrWidth = 0.3,
                    cex.lab=2, font.lab=2, type="p", cex=1.5,  cex.axis=1.6, cex.main=3, col = -1,
                    lend = 1, ylim = c(0, 10), xlim = NULL, xaxt = "n", graybg = F,
                    checks=c(1, 10000000), ...)
{
    grayc = rgb(0.9, 0.9, 0.9)  # light gray
    par(oma=c(0,cex.lab,0,0))
    xlim = c(0, as.integer(lengthChromosome(chr,"hg19")))
    plot(x, y, col =  1,  ylim = ylim, xlim = xlim, xaxt = xaxt, type = "n", ...)

    if(graybg == T)
    {
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = grayc)
        yscale = checks[2]  # every 10M
        abline(v= seq(1, max(xlim), yscale),
             lwd = 3, col ="white") 
        xscale = checks[1]  
        abline(h=seq(0, max(ylim), xscale), lwd = 3,  col = "white")
    }
    if(length(col) == 1 && col == -1)
    {
        col = 1
    }
    points(x, y,  col = col,  ylim = ylim, xlim = xlim, xaxt = "n", type = type, ...)
    theLabel = paste(as.character(as.integer(seq(min(xlim), max(xlim),1e7)/1e6)), "M", sep="")
    axis(1, at=seq(min(xlim), max(xlim),1e7), cex.axis=cex.axis,cex.lab=cex.lab, labels=theLabel)
    chrWidth = (max(ylim) - min(ylim))/11 * chrWidth
    paintCytobands(chr,pos = c(0,min(ylim)), units = "hg19", orientation="h", bleach=0.1,
                   legend=T, width=chrWidth, cex.leg=1.4)
}



#################################################################################
####  This generates the SKY alike, digital SKY.
####          chr  : the id of an autosome, in between 1 and 22
####          lwd  : specifiy the width of the copy number bands.
####          chrWidth : specifiy the width of the ideogram.
####          maxCN : specify biggest copy number value possible.
####          
#################################################################################

plotPloidy <- function(chr, CN, col=2, lwd=0.4, pch=15, chrWidth=0.6, maxCN=12)
{
    drawCN <- function(data, col=2, lwd=10, pch=15)
    {
        start =data[1] * -1
        end  =data[2]  * -1
        cn     =data[3]
        mCN     =data[4]
        if(cn > 0)
        {
            for(i in c(1:cn))
            {
                res = 100
                for(j in 1:res)
                {
                    lines(c(i+lwd/res * j, i+lwd/res * j), c(start,end), col=col, lwd=1,pch=20, lend=1)
                }
            }
        }
        if(cn == 0)
        {
            lwd = 1.2
            res = 60
            for(j in 1:res)
            {
                lines(c(-1.0+lwd/res * j, -1.0+lwd/res * j), c(start,end), col=3, lwd=1, pch=20)
            }
        }
        if(mCN == 0 && cn != 0)
        {
            lwd = 0.2
            res = 30
            for(j in 1:res)
            {
                lines(c(-0.8+lwd/res * j, -0.8+lwd/res * j), c(start,end), col=2, lwd=1, pch=20)
            }
        }
    }

    par(bg = 'Black')
    plot(c(0-chrWidth,maxCN),c(0,lengthChromosome(chr,"hg19"))*(-1),
         type="n",  xaxt="n",yaxt="n",xlab="",ylab="", col.main="white" )
    legend("topright", inset=.03, legend=paste("chr", chr, sep=""), pt.cex =2.5, cex=3, text.col="red")
    paintCytobands(chr, units = "hg19",orientation="v",bleach=0.5, legend=F, width=chrWidth)
    sel = CN[,3]  == 0
    apply(CN[!sel,, drop=F], 1, drawCN, col=col, lwd=lwd, pch=pch)
    if(sum(sel) > 0 & nrow(CN[sel,,drop = F]) > 0)
        apply(CN[sel,,drop = F], 1, drawCN, col=col, lwd=lwd, pch=pch)

}








subplots <- function(nrows, ncols)
{
    par(mfrow=c(nrows, ncols),  
        omi=c(0.05,0.8,0.08,0.00), 
        plt=c(0.10,0.95,0.05,0.65))                 
}




##########################################################################
#  1. Grouping is used as the title for each of the subplots.
#     The number of different groups defines the number of subplots.
#  2. By default all the subplots are laid in one line unless the number 
#     of rows is set to a different number.
##########################################################################
genomeInSubplots <- function(x, y, grouping, chroms = unique(grouping),
                             xaxt = "n", yaxt = "n", ylim=c(0,1), ylab = "", rows = 1)
{
    
    # ‘omi’ A vector of the form ‘c(bottom, left, top, right)’ giving
    #       the size of the outer margins in inches

    # ‘plt’ A vector of the form ‘c(x1, x2, y1, y2)’ giving the
    #      coordinates of the plot region as fractions of the current
    #      figure region.

 #   old.par = par(mfrow=c(rows, floor(length(chroms)/rows)),  
 #                 omi=c(0.05,0.8,0.08,0.00), 
 #                 plt=c(0.10,0.95,0.05,0.85))

    count = 1
    for (chrID in chroms)
    {
        sel = (grouping == chrID)
        if(count == 1)
        {
            #plot(x[sel], y[sel], main = chrID, ylab = "tmp", xlab = "test")
            plot(x[sel], y[sel], main = chrID, xaxt = xaxt, ylim = ylim, ylab = ylab)
            title(outer = T, ylab = ylab)
            #axis(2, col="dodgerblue", col.ticks="green", col.axis="orange", cex.axis=2)
        } else
        {
            plot(x[sel], y[sel], main = chrID, xaxt = xaxt, yaxt=yaxt, ylab = "", ylim = ylim)
        }

        chrWidth = 0.3
    chrWidth = (max(ylim) - min(ylim))/11 * chrWidth
        paintCytobands(1,pos = c(0, min(ylim)), units = "hg19", orientation="h", bleach=0.1, legend=T, width=chrWidth, cex.leg=0.9)
        count = count +1
    }
}

runtest  <-  function()
{

    xx = seq(1, 50000000, 1000000)
    yy = rep(5, length(xx))
    CN = cbind(xx,xx + 100000, yy, yy -1)
    chr = 4
    colForChr = getChrCol(chr)
    pdf(file = "test.pdf", width = 16)
    myplot(xx, yy, chr=chr, col = colForChr)
    myplot(xx, yy, chr=chr, col = colForChr, pch = 4)
    myplot(xx, yy, chr=chr, col = colForChr, pch = 4, type = "o")
    myplot(xx, yy, chr=chr, col = colForChr, pch = 4, type = "l")
    myplot(xx, yy, chr=chr, col = colForChr, pch = 4, type = "l", chrWidth = 0.4)
    myplot(xx, yy, chr=chr, col = colForChr, pch = 4, type = "l", chrWidth = 0.4, ylim = c(-3, 7))

    #  Test: chrWidth, chr color
    par(mfrow = c(1, 2))
    plotPloidy(CN, chr=4, col = getChrCol (4))
    plotPloidy(CN, chr=5, col = getChrCol (5), chrWidth = 1.2)

    #  Test: LOH
    par(mfrow = c(1, 2))
    CN = cbind(xx,xx + 100000, yy, yy - yy )
    plotPloidy(CN, chr=4, col = getChrCol (4))
    plotPloidy(CN, chr=5, col = getChrCol (5))

    #  Test: Double deletions
    par(mfrow = c(1, 2))
    CN = cbind(xx,xx + 100000, yy - yy, yy - yy )
    plotPloidy(CN, chr=4, col = getChrCol (4))
    plotPloidy(CN, chr=5, col = getChrCol (5))

    dev.off()

}

