library(BiSeq)
args <- commandArgs(TRUE)
#inputDir          = args[1]
#compareGroup      = args[2]
#compareGroupName  = args[3]

inputDir          = "/home/meng_yang/project/201704/meth_diff_analysis/BiSeq/cov_bed100_test"
compareGroup      = "LC002-220,LC002-228,LC002-117,LC002-180--LC002-094,LC002-098,LC002-104,LC002-112"
compareGroupName  = "IA--BEN"
sampleName  = unlist(strsplit(compareGroup,"--|,"))
compareGroupNames = unlist(strsplit(compareGroup,"--"))
sampleName1 = unlist(strsplit(compareGroupNames[1],","))
sampleName2 = unlist(strsplit(compareGroupNames[2],","))
groupName1 =  unlist(strsplit(compareGroupName,"--"))[1]
groupName2 =  unlist(strsplit(compareGroupName,"--"))[2]

repNum1 = length(sampleName1)
repNum2 = length(sampleName2)

outdir = paste("out",compareGroupName,sep="/")
dir.create(outdir,recursive = TRUE)
outname = paste(outdir,compareGroupName,sep="/")

###read Bismark data
ppp_colData <-  DataFrame(group=factor(c(rep(groupName1,repNum1),rep(groupName2,repNum2))),row.names=sampleName)
ppp_file <- list.files(path=inputDir, pattern = "cov", full.names = T)
bs <- readBismark(ppp_file,colData = ppp_colData)

###The BSrel class,bs.rel is a container for 'relative' methylation level of BS data
bs.rel <- rawToRel(bs)

#******Quality  control**********
qc <- covStatistics(bs)
write.table(qc, file=paste(outname,".Sample_wise_coverage.stat.xls",sep=""),row.names = TRUE,quote = FALSE,sep="\t")
pngfile= paste(outname,".Sample_wise_coverage_distributions.png",sep="")
png(pngfile)
covBoxplots(bs, col = "cornflowerblue", las = 2)
dev.off()
#******Quality control end*******
######Detection of DMRs within groups of samples######
###a.Definition of CpG clusters.
###Then bs.clust.unlim is a again BSraw object but restricted to CpG sites within CpG clusters
bs.clust.unlim <- clusterSites(object = bs,
                                 groups=colData(bs)$group,
                                 perc.samples = 1/4,
                                 min.sites = 5,
                                 max.dist = 100,
                                 mc.cores = 6)

###b.Smooth methylation data
###Then predictedMeth is a BSrel object with smoothed relative methylation levels for each CpG site within CpG clusters.
ind.cov <- totalReads(bs.clust.unlim) > 0
quant <- quantile(totalReads(bs.clust.unlim)[ind.cov], 1)
print (quant)
bs.clust.lim <- limitCov(bs.clust.unlim, maxCov=quant)
predictedMeth <- predictMeth(object = bs.clust.lim,mc.cores=6)
#******Smooth methylation data plot*********
a <- totalReads(bs.clust.lim)
write.table(a,file=paste(outname,".bs.clust.lim.xls",sep=""),row.names = TRUE,quote = FALSE,sep="\t")
pngfile = paste(outname,".Sample_wise_coverage_distributions_after_coveraget_limitation.png",sep="")
png(pngfile)
covBoxplots(bs.clust.lim, col = "cornflowerblue", las = 2)
dev.off()

#****error****png("Raw_data_together_with_smoothed_methylation_levels.png")
#ppp_region <- GRanges(seqnames="chr1",ranges=IRanges(start=872335,872616))
#plotMeth(object.raw = bs[,6],
#         object.rel = predictedMeth[,6],
#         region = ppp_region,
#         lwd.lines = 2,
#         col.points = "blue",
#         cex = 1.5)
#dev.off()
#******Smooth methylation data end**********
###c.Model and test group effect
###betaResults is a data.frame containing model and test information for each CpG site
#******plot Smoothed methylation levels in tow groups******************
pngfile = paste(outname,".Smoothed_methylation_levels_in_two_groups.png",sep="")
png(pngfile)
treat <- predictedMeth[, colData(predictedMeth)$group == groupName1]
control <- predictedMeth[, colData(predictedMeth)$group == groupName2]
mean.treat <- rowMeans(methLevel(treat),na.rm=TRUE)
mean.control <- rowMeans(methLevel(control),na.rm=TRUE)
plot(mean.control,
     mean.treat,
     col = "blue",
     xlab = paste("Methylation in ",groupName2,sep=""),
     ylab = paste("Methylation in ",groupName1,sep=""))
dev.off()

#******plot Smoothed methylation levels in tow groups end******************
betaResults <- betaRegression(formula = ~group,
                              link = "probit",
                              object = predictedMeth,
                              type = "BR",
			      mc.cores=6)
###or
###d.Test CpG clusters for differential methylation

predictedMethNull <- predictedMeth[,c(1:min(repNum1,repNum2),(repNum1+1):(repNum1+min(repNum1,repNum2)))]
colData(predictedMethNull)$group.null <- rep(c(1,2), min(repNum1,repNum2))
betaResultsNull <- betaRegression(formula = ~group.null,
                                  link = "probit",
                                  object = predictedMethNull,
                                  type="BR",
				  mc.cores=6)
###or
#estimate the variogram for Z score obtained for the resampled data
vario <- makeVariogram(betaResultsNull)
#plot(vario$variogram$v)
vario.sm <- smoothVariogram(vario, sill = 0.9)
#lines(vario.sm$variogram[,c("h", "v.sm")],col = "red", lwd = 1.5)
#grid()
#We replace the pValsList object (which consists of the test results of the resampled data) by the test results of interest (for group eect):\
vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)
vario.sm$pValsList <- vario.aux$pValsList
#The correlation of the Z scores between two locations in a cluster can now be estimated
locCor <- estLocCor(vario.sm)

for (i in c(0.9,0.5,0.4,0.3,0.2,0.1)){
	for(j in c(1,0.9,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05)){
		tmp = paste(i,j,sep="-")
		tmp0 = paste(tmp,".clusters.trimmed.xls",sep="")
		tmp1 = paste(tmp,".clusters.reject.xls",sep="")
		tmp = paste(tmp,".dmr.xls",sep="")
		clusters.rej <- testClusters(locCor, FDR.cluster = i)
		clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = j)
		DMRs <- findDMRs(clusters.trimmed, max.dist = 100, diff.dir = TRUE)
		write.table(clusters.rej$clusters.reject,file=paste(outname,tmp1,sep="."),row.names = FALSE,quote = FALSE,sep="\t")
		write.table(clusters.trimmed,file=paste(outname,tmp0,sep="."),row.names = FALSE,quote = FALSE,sep="\t")		
		write.table(DMRs,file=paste(outname,tmp,sep="."),row.names = FALSE,quote = FALSE,sep="\t")		
	}
}


#####clusters.rej <- testClusters(locCor, FDR.cluster = 0.2)

###e.Trim significant CpG clusters
#Then clusters.trimmed is a data.frame object containing all differentially methylated CpG sites.
#####clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = 0.1)
#####write.table(clusters.trimmed,"clusters.trimmed.xls",row.names = FALSE,quote = FALSE,sep="\t")
###f.Definition of DMR boundaries
########DMRs <- findDMRs(clusters.trimmed, max.dist = 100, diff.dir = TRUE)

#Definition of DMR boundaries
########write.table(DMRs,"dmr.xls",row.names = FALSE,quote = FALSE,sep="\t
