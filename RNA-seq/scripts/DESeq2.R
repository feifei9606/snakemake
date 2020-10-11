library(DESeq2)
library(optparse)
option_list = list(
  make_option(c("-c", "--counts"), type="character", default="mapped/count/readcount/count.txt", help="counts tabulated file from Feature Counts", metavar="character"),
  make_option(c("-s", "--samplefile"), type="character", default="samples.tsv", help="sample files used to get conditions for DESEq2 model fit", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="Rplot/DESeq2_results_all.csv", help="where to place differential expression files", metavar="character")
) 

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load countmatrix
countdata <- read.table(opt$counts, header=TRUE, skip=1, row.names=1,stringsAsFactors = F, check.names = F)
countdata <- countdata[ ,6:ncol(countdata)]

colnames(countdata) <- gsub("mapped/[sb]am/(.*?).[sb]am$","\\1",colnames(countdata))
countdata <- countdata[-which(rowSums(countdata)<4),]
countdata <- as.matrix(countdata)

# experimental conditions 
samplefile <- read.delim(file = opt$samplefile,header = T,stringsAsFactors = F)
condition <- factor(x = samplefile[,"condition"],levels = unique(samplefile[,"condition"]))

## DESeq2
colData <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~condition)
dds <- DESeq(dds)
res <- results(dds) 
summary(res)
resOrdered <- res[order(res$padj), ] 

write.table(resOrdered, file=opt$outdir,sep = "\t",quote=F)
#deg <- subset(resOrdered, padj <= 0.01 & abs(log2FoldChange) >= 2) 
#write.table(deg, file=opt$outdir+"/DESeq2_results_significant.csv",sep = "\t",quote=F)

## volcano plot
# library(ggplot2)
# volcano_data <- read.table("DESeq2_results_all.csv",sep = "\t")
# volcano_data <- na.omit(volcano_data) 
# significant <- as.factor(abs(volcano_data$log2FoldChange) >=2 & volcano_data$padj <= 0.01) 
# ggplot(volcano_data, aes(x = log2FoldChange, y = - log10(padj))) + geom_point(aes(shape = significant, color = significant)) + xlim(c(-10, 10)) + labs(x = "log2FoldChange", y = "-log10 padj") + scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) + scale_shape_discrete(labels =c ("no", "yes")) + scale_color_discrete(labels = c("no", "yes"))
# options(bitmapType = 'cairo', device = 'png')
# ggsave("Rplot/volcano.pdf",limitsize = FALSE)
