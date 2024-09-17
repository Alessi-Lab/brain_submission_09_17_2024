a <- read.table("H:/20230905_150637_Total_Proteome.wide.protein.quant.2.tsv", sep = "\t", header=TRUE,comment.char = "", na.strings = c("#VALUE!"))
sample.annotation <-read.table("H:/SAMPLE.INFO.BASICS_v5.txt", sep="\t", header=TRUE, comment.char = "")
library(QFeatures)
library(stringr)
a <- a[a["PG.NrOfStrippedSequencesIdentified..Experiment.wide."]>1,]
data <- readQFeatures(a, ecol = rep(3:157), name="PG.ProteinGroups")
condition <- c()
sample <- c()
rename.cols <- c()
for (i in colnames(a)[3:157]) {
  x <- unlist(strsplit(i, "\\."))
  sample <- cbind(sample, x[4])
  d <- sample.annotation[sample.annotation["Sample.ID"]==x[4]]
  condition <- cbind(condition, paste(str_replace(d[4], " ", ""), str_replace(d[5], " ","")))
  rename.cols <- cbind(rename.cols, paste(str_replace(d[4], " ", ""), str_replace(d[5], " ",""), x[4], sep = "."))
}

library(data.table)
library(limma)
setnames(a, old=colnames(a)[3:157], new=rename.cols)

data$group <- condition
data$sample <- sample
colData(data)
rowDataNames(data)
d <- selectRowData(data, c("PG.ProteinGroups"))
d <- zeroIsNA(d, i = seq_along(d))
nNA(d, i = seq_along(d))
d <- filterNA(d, i = seq_along(d), pNA = 0.7)
d <- impute(d, method="knn", i="PG.ProteinGroups")
d <- addAssay(d, logTransform(d[[2]]), name="log2")
d <- addAssay(d,
              normalize(d[["log2"]], method = "quantiles.robust"),
              name = "norm")
normalized.df <- cbind(as.data.frame(rowData(d[["norm"]])), as.data.frame(assay(d, "norm")))

setnames(normalized.df, old=colnames(normalized.df)[2:156], new=rename.cols)

write.table(normalized.df, "H:/20230905_150637_Total_Proteome.wide.protein.quant.normalized.0.7.tsv", sep="\t", row.names = FALSE)


design <- model.matrix(~0+d$group)
colnames(design) <- gsub("d\\$group", "", colnames(design))
colnames(design) <- str_replace(colnames(design), " ", ".")
fit <- lmFit(assay(d, "norm"), design)
contrast.matrix <- makeContrasts(L2PD.G2019S-iPD.., levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
write.table(cbind(rowData(d[[4]]), result), "H:/20230905_150637_Total_Proteome.wide.protein.quant.de.0.7.tsv", sep="\t", row.names = FALSE)
write.table(a, "H:/20230905_150637_Total_Proteome.wide.protein.quant.raw.tsv", sep="\t", row.names = FALSE)