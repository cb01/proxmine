#!/usr/bin/Rscript
#
# reads in a "mini sam" e.g. a SAM cut down to the columns
# genome 1, position, genome 2, position, insert distance (template len)
# first map reads with: 
# bwa mem Hi_Bact.read1.fq.gz Hi_Bact.read2.fq.gz | samtools view -S -b - | samtools sort - hi_bact
# then generate minisam with: 
# samtools view -q 20 hi_bact.bam | cut -f 3,4,7,8,9 > minisam.txt
#
# copyleft 2013 Aaron Darling
#

minisam <- read.table(commandArgs(TRUE)[1])
# only plot things that are probably long distance Hi-C reads
#minisam <- minisam[minisam$V5!=0]

locus <- levels(unique(minisam$V1))
if(locus[1] == "*"){
	locus<-locus[2:length(locus)]
}

lnames <- list(
"gi|116491818|ref|NC_008525.1|"="P. pentosaceus ATCC25745",
"gi|170079663|ref|NC_010473.1|"="E. coli K12 DH10B",
"gi|387823261|ref|NC_012892.2|"="E. coli BL21(DE3)",
"gi|116332681|ref|NC_008497.1|"="L. brevis ATCC367",
"gi|116334867|ref|NC_008498.1|"="L. brevis plasmid 1",
"gi|116334879|ref|NC_008499.1|"="L. brevis plasmid 2",
"gi|83716035|ref|NC_007650.1|"="B. thailandensis E264 chrI",
"gi|83718394|ref|NC_007651.1|"="B. thailandensis E264 chrII"
)

lengths <- list(
"gi|116491818|ref|NC_008525.1|"=1832387,
"gi|170079663|ref|NC_010473.1|"=4686137,
"gi|387823261|ref|NC_012892.2|"=4558947,
"gi|116332681|ref|NC_008497.1|"=2291220,
"gi|116334867|ref|NC_008498.1|"=13413,
"gi|116334879|ref|NC_008499.1|"=35595,
"gi|83716035|ref|NC_007650.1|"=2914771,
"gi|83718394|ref|NC_007651.1|"=3809201
)


# make all the insert distributions
insdists <- list()
densities <- list()
max_reads <- 0
for(i in 1:length(locus)){
	insdists[[i]] <- abs(minisam$V5[minisam$V1 == locus[i]& minisam$V3 == "="])
	insdists[[i]] <- insdists[[i]][insdists[[i]]>1000]
	cur_len <- as.numeric(lengths[locus[i]])
	insdists[[i]][ insdists[[i]] >  cur_len / 2 ] <- cur_len - insdists[[i]][ insdists[[i]] > cur_len / 2 ]
	densities[[i]] <- hist(insdists[[i]], xlim=c(0,2500000), breaks=seq(0,2500000,length.out=100))
	max_reads <- max(max_reads, densities[[i]]$counts)
	# following line eliminates the edge-effect
	densities[[i]]$counts[ densities[[i]]$counts == 0 ] <- NA
}

for(i in 1:length(locus)){
	densities[[i]]$counts[ min(which(is.na(densities[[i]]$counts))) - 1 ] <- NA
}

pdf("ins_dist.pdf",width=8,height=5)
plot(y=densities[[1]]$counts, x=densities[[1]]$mids, main="Hi-C read pair insert distributions",ylim=c(0,max_reads),xlim=c(0,2750000),ylab="Count",xlab="Insert distance in nucleotides",type="l")
for(i in 2:length(locus)){
	lines(y=densities[[i]]$counts, x=densities[[i]]$mids, col=i)
}

legnames <- c(lnames[locus[1]], lnames[locus[4:length(locus)]])

legend("topright",legend=legnames, col=c(1,4,5,6,7,8),lty=rep(1,(length(lnames)-2)))
dev.off()

