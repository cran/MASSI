### required libraries
library(fpc)
library(gplots)

### MASSI espresso function
massi.espresso <- function(exprs, y.probes, threshold=3){
  
  #Check that the input data is in the correct data.frame class
  class.exprs <- class(exprs)
  class.y.probes <- class(y.probes)
  if(class.exprs != "data.frame") stop("input data must be in data.frame class")
  if(class.y.probes != "data.frame") stop("input data must be in data.frame class")
  
  exprs$ID <- rownames(exprs) # set probe as ID
  y.probes$ID <- row.names(y.probes) # set probe as ID
  y.values <- as.data.frame(merge(exprs, y.probes, by="ID")) # extract matched probes from expression matrix using ID
  cal.cv <- function(x) ( 100*sd(x)/mean(x) ) # define function to calculate CV
  y.values$CV <- apply(y.values[, -which(names(y.values) == "ID")], MARGIN=1, FUN=cal.cv) # calculate CV for each probe
  quantiles <- quantile(y.values$CV) # calculate quantiles for probe CV
  cv.threshold <- quantiles[threshold] # set threshold (4=75%, 3=50% ,2=25%, 1=all)
  cv.threshold <<- quantiles[threshold] 
  cv.cutoff <- function(x) ((x)>=cv.threshold) # define function to select probes above threshold
  y.values$above.threshold<- cv.cutoff(x=y.values$CV) # identify probes above threshold
  y.values <<- y.values
  
  y.subset.values <- as.data.frame(subset(y.values[, -which(names(y.values) == "CV")],
                                          subset=y.values$above.threshold == TRUE)) #extract probes above threshold
  
  # change data frame into format for cluster analysis
  y.subset.values$above.threshold <- NULL
  row.names(y.subset.values) <- NULL
  row.names(y.subset.values) <- y.subset.values$ID
  y.subset.values$ID <- NULL
  y.subset.values <<- y.subset.values
  
  ## generate probe CV plot
  cv.max = (round(1.1*(max(y.values$CV)))+1)
  probe.names <- y.values$ID
  
  barplot.default(height=y.values$CV, ylab="Probe CV (%)", 
                  xlab="", names.arg=probe.names, col="blue",
                  ylim=c(0,cv.max), las=2, xpd=F, cex.names=0.5)
  title(xlab="Y chromosome probes", line=4)
  abline(h=cv.threshold, col="red",)
  
  ## perfom clustering using k-medoids
  y.subset.values.t <- data.frame(t(y.subset.values)) #transpose the probe subset values matrix
  max.clusters <- as.numeric(nrow(y.subset.values.t))-1 # set the maximum number of clusters to n-1
  clust.num <- pamk(data=y.subset.values.t, krange=c(2,max.clusters)) ### Find the optimum number of clusters 2:n-1
  pam.k = as.numeric(clust.num$nc) ## set number of clusters as numeric
  ifelse(test=pam.k==2, yes="Good, Optimum number of clusters is 2, proceed to massi.cluster", no="Bad, Optimum number of clusters is >2, try reducing the Y chromosome probe list by setting a higher threshold")
  
  y.subset.values.t <- data.frame(t(y.subset.values)) #transpose the probe subset values matrix
  y.kmedoids <- pam(x=y.subset.values.t, k=2) ### perform clustering of pam.k clusters
  y.kmedoids <<- pam(x=y.subset.values.t, k=2) 
  sample.sex <- data.frame(y.kmedoids$clustering) ## get samples and clusters into data.frame
  
  ## now the samples are classified, the group with the highest mean is assigned as male
  sample.means <- data.frame(colMeans(y.subset.values)) # calculate the mean of probe values for each sample
  sample.sex <- merge(sample.sex, sample.means, by="row.names") # merge the clusters and mean expression by sample ID
  sample.sex$ID <- sample.sex$Row.names
  sample.sex$Row.names <- NULL
  
  ## calculate the mean for each cluster of the NRY probe mean
  cluster.1 <- subset(sample.sex, subset=sample.sex$y.kmedoids.clustering==1) #subset cluster 1
  cluster.1.mean <- mean(cluster.1$colMeans.y.subset.values.) # calculate cluster 1 mean
  
  cluster.2 <- subset(sample.sex, subset=sample.sex$y.kmedoids.clustering==2) # subset cluster 2
  cluster.2.mean <- mean(cluster.2$colMeans.y.subset.values.) # calculate cluster 2 mean
  
  # assign cluster with the highest mean as male, and lowest mean as female
  c1.sex <- ifelse(cluster.1.mean>cluster.2.mean, yes=as.character("male"), no=as.character("female"))
  c2.sex <- ifelse(cluster.1.mean<cluster.2.mean, yes=as.character("male"), no=as.character("female"))  
  
  # create a column for sex and substitute the cluster id for "Male" or "Female"
  sample.sex$kmedoids.sex <- as.character(sample.sex$y.kmedoids.clustering) # create the column for sex
  sample.sex$kmedoids.sex[sample.sex$kmedoids.sex == "1"] <- c1.sex
  sample.sex$kmedoids.sex[sample.sex$kmedoids.sex == "2"] <- c2.sex
  
  sample.sd <- data.frame(sapply(y.subset.values, FUN=sd)) ## add sample sd to output
  
  sample.sd$ID <- row.names(sample.sd)
  #massi.output <- merge(x=massi.output, y=sample.sd, by="sample.ID")
  #massi.output <<- massi.output ## output of clustering results
  #write.table(x=massi.output, file="MASSI.results.txt", quote=F, sep="\t", row.names=F)
  
  massi.results <- data.frame(sample.sex$ID) # Add sample ID
  massi.results$ID <- massi.results$sample.sex.ID
  massi.results$mean.probe.value <- sample.sex$colMeans.y.subset.values. # add mean probe values
  massi.results<- merge(massi.results, sample.sd, by="ID") # add sample/probe sd
  massi.results$sample.sd <- massi.results$sapply.y.subset.values..FUN...sd.
  massi.results$sapply.y.subset.values..FUN...sd. <- NULL
  
  massi.results <- merge(massi.results, sample.sex, by="ID") # add clustering results
  massi.results$sample.sex.ID <- NULL # remove redundant fields
  massi.results$y.kmedoids.clustering <- NULL
  massi.results$colMeans.y.subset.values. <- NULL
  
  massi.output <- massi.results
  massi.output <<- massi.results
  
  write.table(massi.results, file="massi.results.txt", quote=F, sep="\t", row.names=F)
  print(massi.results)
  
  my.hclust <- function(d) hclust(d, method="ward")
  
  pdf(file="MASSI.figures.pdf", onefile=T, paper="a4r", useDingbats=F)
  heatmap.2(x=as.matrix(y.subset.values), hclustfun=my.hclust, keysize=2, cexRow=0.7,
            key=T, trace="none", dendrogram="column", col=redgreen(75), scale="row")
  
  
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
                        xlab="Samples",
                        col=c("red", "green")[as.factor(sample.sex)],
                        las=2, ylim=c(plot.bottom,plot.top))
  
  segments(barCenters, probe.means-probe.sd, # add the sd bars
           barCenters, probe.means+probe.sd, lwd=0.8)
  
  
  legend(x=1, y=plot.top, fill=c("red", "green"), title="MASSI sex", ## add legend to plot
         legend=c("female", "male"), cex=1, )
  
  
  clusplot(t(y.subset.values), y.kmedoids$clustering, color=TRUE, shade=FALSE, main="",cex.txt=0.5,
           labels=2, lines=0)
  dev.off()  
}
