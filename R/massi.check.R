library("fpc")

massi.check <- function(exprs, y.probes, threshold=3) {
    
  exprs$ID <- rownames(exprs) # set probe as ID
  y.probes$ID <- row.names(y.probes) # set probe as ID
  y.values <- as.data.frame(merge(exprs, y.probes, by="ID")) # extract matched probes from expression matrix using ID
  cal.cv <- function(x) ( 100*sd(x)/mean(x) ) # define function to calculate CV
  y.values$CV <- apply(y.values[, -which(names(y.values) == "ID")], MARGIN=1, FUN=cal.cv) # calculate CV for each probe
  quantiles <- quantile(y.values$CV) # calculate quantiles for probe CV
  cv.threshold <- quantiles[threshold] # set threshold (4=75%, 3=50% ,2=25%, 1=all)
  cv.threshold <<- quantiles[threshold] # set threshold (4=75%, 3=50% ,2=25%, 1=all)
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
}

