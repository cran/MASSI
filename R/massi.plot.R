library(gplots)
library(fpc)

massi.plot <- function(y.values=y.values, y.subset.values=y.subset.values,
                       massi.output=massi.output) {
  
  ## Null variable to pass check
  
  
  ## generate heatmap of Y chromosome probe subset
  ord <- order(rowSums(abs(y.subset.values)),decreasing=T)
  pdf(file="MASSI.figures.pdf", onefile=T, paper="a4r", useDingbats=F)
  heatmap.2(x=as.matrix(y.subset.values[ord,]), keysize=2, cexRow=0.7,
            key=T, trace="none", dendrogram="row", col=redgreen(75), scale="row")
  
  ## generate mean and SD expression barplot
  ## plot values
  massi.output.sort <- massi.output[order(massi.output$kmedoids.sex),] # sort data by sex
  probe.means <- massi.output.sort$mean.probe.value # samples probe mean values
  probe.sd <- massi.output.sort$sample.sd # sample probe sd values
  sample.names <- massi.output.sort$sample.ID # set x-axis names
  plot.top <- ceiling(max(probe.means+probe.sd*1.1)) # set y-axis upper limit
  plot.bottom <- floor(min(probe.means-probe.sd*1.1)) # set y-axis lower limit
  sample.sex <- massi.output.sort$kmedoids.sex # set the factor for bar color
  # create the plot
  barCenters <- barplot(probe.means, xpd=F, names.arg=massi.output$ID, cex.names=0.7,
                        ylab="Chr.Y probe expression (Mean +/- SD)",
                        xlab="",
                        col=c("red", "green")[as.factor(sample.sex)],
                        las=2, ylim=c(plot.bottom,plot.top))
  segments(barCenters, probe.means-probe.sd, # add the sd bars
           barCenters, probe.means+probe.sd, lwd=0.8)
  legend(x=1, y=plot.top, fill=c("red", "green"), title="MASSI sex", ## add legend to plot
         legend=c("female", "male"), cex=1, )
  
  ## generate PC plot of clusters
  y.kmedoids <- get("y.kmedoids")
  clusplot(t(y.subset.values), y.kmedoids$clustering, color=TRUE, shade=FALSE, main="",cex.txt=0.5,
           labels=2, lines=0)

  
  dev.off()
}
