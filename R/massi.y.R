
massi.y <- function(exprs, y.probes){
  
  exprs$ID <- rownames(exprs) # set probe as ID
  y.probes$ID <- row.names(y.probes) # set probe as ID
  y.values <- as.data.frame(merge(exprs, y.probes, by="ID")) # extract matched probes from expression matrix using ID
  cal.cv <- function(x) ( 100*sd(x)/mean(x) ) # define function to calculate CV
  y.values$CV <- apply(y.values[, -which(names(y.values) == "ID")], MARGIN=1, FUN=cal.cv) # calculate CV for each probe
  quantiles <- quantile(y.values$CV) # calculate quantiles for probe CV
  ## generate probe CV plot
  cv.max = (round(1.1*(max(y.values$CV)))+1)
  probe.names <- y.values$ID
  y.plot <- function(){
    barplot.default(height=y.values$CV, ylab="Probe CV (%)", 
                    xlab="", names.arg=probe.names, col="blue",
                    ylim=c(0,cv.max), las=2, xpd=F, cex.names=0.5)
    title(xlab="Y chromosome probes", line=4)
    abline(h=quantiles[2], col="red", lty=2)
    abline(h=quantiles[3], col="red", lty=2)
    abline(h=quantiles[4], col="red", lty=2)
    # draw an axis on the right, with text and ticks 
    axis(4, line=-0.5, at=c(quantiles[1], quantiles[2], quantiles[3],quantiles[4]),
         labels=c(1,2,3,4), tick=T, las=2)
    mtext("Threshold", side=4, line=1, las=3, cex.lab=1, at=quantiles[3], )
  }
  
  ## save y.plot as .pdf
  pdf(file="massi.yplot.pdf", paper="a4r")
  y.plot()
  dev.off()
## display y.plot
  dev.new()
  y.plot()
}
