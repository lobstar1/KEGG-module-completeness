#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-k", "--kegg"), type="character", action = "store", help = "gene2module (pathway) tsv file"),
  make_option(c("-g", "--gene2ko"), type="character", action = "store", help = "gene2ko (kegg ortholog) tsv file"),
  make_option(c("-s", "--selected"), type="character", action = "store", help = "selected gene lists for enrichment analysis"),
  make_option(c("-o", "--output"), type="character", action = "store", default = "Module", help = "Output prefix name, default Module"),
  make_option(c("-p", "--padj"), type="double", default = 0.05, action = "store", help = "Expected padj value for the DEGs, default 0.05"),
  make_option(c("-P", "--pvalue"), default = FALSE, action = "store_true", help = "Use pvalue rather than padj, default fause"),
  make_option(c("-c", "--count"), type="character", action = "store", default = "NA", help = "Count table to generate DEGs, NA means equal threshold for each genes to represent KO, default NA")
 )

opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE, usage = "usage: %prog [options]"))



library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(ggupset)
Gene2kegg <- read.table(opt$kegg,header=F,sep="\t")
keggAnn <- read.table("/nfs_genome/BioDB/kegg/MO_name.tsv",header=F,sep="\t")

if(opt$count != "NA")
{
	count<-read.table(opt$count, header=T,sep="\t", row.names=1)
	threshold<-log(rowSums(count[,6:dim(count)[2]]),2)%/%1
} else
{
	threshold<-rep(1,dim(Gene2kegg)[1])
	names(threshold)<-Gene2kegg[,1]
}

diff <- read.table(opt$selected,header=T,sep="\t",row.names=1)

if(opt$pvalue)
{
        num = 5
	x =enricher(rownames(diff),TERM2GENE=Gene2kegg[,c(2,1)],TERM2NAME=data.frame(name=keggAnn$V1,def=keggAnn$V2),pvalueCutoff=1, qvalueCutoff=1)
} else
{
        num = 6
	x =enricher(rownames(diff),TERM2GENE=Gene2kegg[,c(2,1)],TERM2NAME=data.frame(name=keggAnn$V1,def=keggAnn$V2),pvalueCutoff=opt$padj, qvalueCutoff=(opt$padj*4))
}
enrichRes <- x[x[,num]<=opt$padj]

num_sig <- dim(enrichRes)[1]

head(enrichRes)
write.table(enrichRes, file = paste0(opt$output,"_mdEnriched.tsv"), quote = FALSE, sep="\t")


x2 <- pairwise_termsim(x)
if(num_sig >= 50)
{
        num_show = 50
} else
{
        num_show = num_sig
}

library(RColorBrewer)
col <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(256))

pdf(paste0(opt$output,"_dotplot.pdf"),width=15,height=30)
dotplot(x,showCategory=num_show)
dev.off()
pdf(paste0(opt$output,"_emetplot.pdf"),width=20,height=20)
emapplot(x2,showCategory=num_show)
dev.off()
pdf(paste0(opt$output,"_cnetplot.pdf"),width=20,height=20)


fc<-diff$log2FoldChange
#try(names(fc)<-rownames(diff))
if(is.null(fc)){
fc <- rep(1, dim(diff)[1])
names(fc)<-rownames(diff)
} else {
names(fc)<-rownames(diff)
switch = 1

}

fc2 <- fc
fc2[fc2 > 5] = 5
fc2[fc2 < -5] = -5

#if !(is.null(fc)){
cnetplot(x2,foldChange=fc2,showCategory=num_show,node_label="category") + scale_color_gradientn(colours=col)
dev.off()
pdf(paste0(opt$output,"_heatplot.pdf"),width=30,height=12)

heatplot(x2,foldChange=fc2,showCategory=num_show) + scale_fill_gradientn(colours=col)
dev.off()


pdf(paste0(opt$output,"_upsetplot.pdf"),width=30,height=12)

upsetplot(x2,n=num_show)
dev.off()

