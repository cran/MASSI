### R code from vignette source 'MASSI_vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: load the example data
###################################################
library(MASSI)
data(MASSI.test.dataset)


###################################################
### code chunk number 2: load the test probe list
###################################################
data(MASSI.test.probes)


###################################################
### code chunk number 3: load data from ExpressionSet (eval = FALSE)
###################################################
## norm.log.values <- exprs(object=ExpressionSet)


###################################################
### code chunk number 4: extract Y from test dataset
###################################################
massi.y(exprs=MASSI.test.dataset, y.probes=MASSI.test.probes)


###################################################
### code chunk number 5: fig1too
###################################################
exprs <- MASSI.test.dataset
y.probes <- MASSI.test.probes
exprs$ID <- rownames(exprs)
    y.probes$ID <- row.names(y.probes)
    y.values <- as.data.frame(merge(exprs, y.probes, by = "ID"))
    cal.cv <- function(x) (100 * sd(x)/mean(x))
    y.values$CV <- apply(y.values[, -which(names(y.values) == 
        "ID")], MARGIN = 1, FUN = cal.cv)
    quantiles <- quantile(y.values$CV)
    cv.max = (round(1.1 * (max(y.values$CV))) + 1)
    probe.names <- y.values$ID
    y.plot <- function() {
        barplot.default(height = y.values$CV, ylab = "Probe CV (%)", 
            xlab = "", names.arg = probe.names, col = "blue", 
            ylim = c(0, cv.max), las = 2, xpd = F, cex.names = 0.5)
        title(xlab = "Y chromosome probes", line = 4)
        abline(h = quantiles[2], col = "red", lty = 2)
        abline(h = quantiles[3], col = "red", lty = 2)
        abline(h = quantiles[4], col = "red", lty = 2)
        axis(4, line = -0.5, at = c(quantiles[1], quantiles[2], 
            quantiles[3], quantiles[4]), labels = c(1, 2, 3, 
            4), tick = T, las = 2)
        mtext("Threshold", side = 4, line = 1, las = 3, cex.lab = 1, 
            at = quantiles[3], )
    }



###################################################
### code chunk number 6: fig1
###################################################
y.plot()


###################################################
### code chunk number 7: load the included probe lists
###################################################
 data(y.probes)


###################################################
### code chunk number 8: run massi.check
###################################################
massi.check(exprs=MASSI.test.dataset, y.probes=MASSI.test.probes,
            threshold=3)


###################################################
### code chunk number 9: MASSI_vignette.Rnw:121-122
###################################################
results <- massi.cluster(y.subset.values=y.subset.values)


###################################################
### code chunk number 10: MASSI_vignette.Rnw:124-125
###################################################
head(results)


###################################################
### code chunk number 11: MASSI_vignette.Rnw:134-136
###################################################
massi.plot(y.values=y.values,
           y.subset.values=y.subset.values, massi.output=massi.output)


###################################################
### code chunk number 12: fig2too
###################################################
massi.output.sort <- massi.output[order(massi.output$kmedoids.sex),]
    probe.means <- massi.output.sort$mean.probe.value
    probe.sd <- massi.output.sort$sample.sd
    sample.names <- massi.output.sort$sample.ID
    plot.top <- ceiling(max(probe.means + probe.sd * 1.1))
    plot.bottom <- floor(min(probe.means - probe.sd * 1.1))
    sample.sex <- massi.output.sort$kmedoids.sex


###################################################
### code chunk number 13: fig2
###################################################
barCenters <- barplot(probe.means, xpd = F, names.arg = massi.output$ID, 
        cex.names = 0.7, ylab = "Chr.Y probe expression (Mean +/- SD)", 
        xlab = "Samples", col = c("red", "green")[as.factor(sample.sex)], 
        las = 2, ylim = c(plot.bottom, plot.top))
segments(barCenters, probe.means - probe.sd, barCenters, 
        probe.means + probe.sd, lwd = 0.8)
    legend(x = 1, y = plot.top, fill = c("red", "green"), title = "MASSI sex", 
        legend = c("female", "male"), cex = 1, )


###################################################
### code chunk number 14: fig3too
###################################################
ord <- order(rowSums(abs(y.subset.values)), decreasing = T)


###################################################
### code chunk number 15: fig3
###################################################
require(gplots)
heatmap.2(x = as.matrix(y.subset.values[ord, ]), keysize = 2, 
    cexRow = 0.7, key = T, trace = "none", dendrogram = "row", 
    col = redgreen(75), scale = "row")


###################################################
### code chunk number 16: fig4too
###################################################
y.kmedoids <- get("y.kmedoids")


###################################################
### code chunk number 17: fig4
###################################################
clusplot(t(y.subset.values), y.kmedoids$clustering, color = TRUE, 
        shade = FALSE, main = "", cex.txt = 0.5, labels = 2, 
        lines = 0)


###################################################
### code chunk number 18: massi.espresso function
###################################################
data(MASSI.test.dataset)
data(MASSI.test.probes)
massi.espresso(exprs=MASSI.test.dataset, y.probes=MASSI.test.probes)


