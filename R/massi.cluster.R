require("fpc")

massi.cluster <- function(y.subset.values){
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
  
  # add sample average z-score
  z.score <- scale(t(y.subset.values))
  z.score.mean <- data.frame(rowMeans(z.score))
  z.score.mean$ID <- row.names(z.score.mean)
  massi.results <- merge(massi.results, z.score.mean, by="ID")
  
  massi.results <- merge(massi.results, sample.sex, by="ID") # add clustering results
  massi.results$sample.sex.ID <- NULL # remove redundant fields
  massi.results$y.kmedoids.clustering <- NULL
  massi.results$colMeans.y.subset.values. <- NULL
  
  massi.output <- massi.results
  massi.output <<- massi.results
  
  write.table(massi.results, file="massi.results.txt", quote=F, sep="\t", row.names=F)
  print(massi.results)
}


