
# Usage: R CMD BATCH genomic-heatmap.R [link data] [output dir] [naming tag 1] [naming tag 2]
# E.g. R CMD BATCH genomic-heatmap.R links.txt ~/test/ Lac0 Lac1

args<-commandArgs(TRUE)

# Read params and data
f=read.table(file=args[1])
outdir = paste(args[2], "/figs", sep="")
tag1=args[3]
tag2=args[4]
library(lattice)
f2=as.matrix(f)

# Generate genomic heatmap
pdf(file=paste(outdir, "/heat-genomic-", tag1, "-", tag2, ".pdf", sep=""), height=10, width=10)
levelplot(log(as.matrix(f2)+0.01), col.regions=colorRampPalette(c("blue","red","yellow")), scales=list(y=list(draw=FALSE), x=list(draw=FALSE)), xlab="", ylab="")
dev.off()

# Compute spearman transformation
f2=as.matrix(f)
f3=f2
for(i in 1:dim(f)[1]){
for(j in 1:dim(f)[1]){
f2[i,j]=cor(f[,i], f[,j], method="spearman")
}}

# Generate spearman-transformed genomic heatmap
pdf(file=paste(outdir, "/heat-genomic-spearman", tag1, "-", tag2, ".pdf", sep=""), height=10, width=10)
levelplot(as.matrix(f2), col.regions=colorRampPalette(c("blue","red","yellow")), scales=list(y=list(draw=FALSE), x=list(draw=FALSE)), xlab="", ylab="")
dev.off()



