library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])
# Get the sample order so that we can annotate the PC table
sample_order = read.table(args[2], header = F)
colnames(sample_order) = c("sample")

# Read in the covariance matrix made by PCAngsd
cov <- as.matrix(read.table(args[3]))

# Decompose to get the PCs
e<-eigen(cov)

# Set the column and row names of the PCs before writing out
rownames(e$vectors) = sample_order$sample
colnames(e$vectors) = paste0("PC", 1:ncol(e$vectors))

write.table(e$vectors, "PCs.txt")

e.values = e$values/(sum(e$values))*100

write.table(e.values, "PC.values.txt", row.names = F, col.names = F)

pc.plot = ggplot(as.data.frame(e$vectors), aes(x=PC1, y=PC2)) + geom_point() + xlab(paste0("PC1 ", format(round(e.values[1], 2), nsmall = 2), "%")) + ylab(paste0("PC2 ", format(round(e.values[2], 2), nsmall = 2), "%"))
ggsave("PC.plot.png")