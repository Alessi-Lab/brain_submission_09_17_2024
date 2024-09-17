a <- read.table("H:/20230905_150637_Total_Proteome.wide.protein.quant.2.tsv", sep = "\t", header=TRUE,comment.char = "", na.strings = c("#VALUE!"))
sample.annotation <-read.table("H:/SAMPLE.INFO.BASICS_v5.txt", sep="\t", header=TRUE, comment.char = "")
library(QFeatures)
library(stringr)
a <- a[a["PG.NrOfStrippedSequencesIdentified..Experiment.wide."]>1,]
data <- readQFeatures(a, ecol = rep(3:157), name="PG.ProteinGroups")
selectedSamples <- c()
condition <- c()
sample <- c()
rename.cols <- c()
for (i in colnames(a)[3:157]) {
  x <- unlist(strsplit(i, "\\."))
  sample <- cbind(sample, x[4])
  d <- sample.annotation[sample.annotation["Sample.ID"]==x[4]]
  if (str_starts(x[4], "B")) {
    if (d[5] != ".") {
      condition <- cbind(condition, c("mutant"))
      rename.cols <- cbind(rename.cols, paste("mutant", x[4], sep = "."))
    } else {
      condition <- cbind(condition, paste(str_replace(d[4], " ", "")))
      rename.cols <- cbind(rename.cols, paste(str_replace(d[4], " ", ""), x[4], sep = "."))
    }
    selectedSamples <- c(selectedSamples, i)
  }
}
e <- a[c("PG.ProteinGroups", selectedSamples)]

data <- readQFeatures(e, ecol = selectedSamples, name="PG.ProteinGroups")
library(data.table)
setnames(e, old=colnames(e)[2:70], new=rename.cols)

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

library("limma")
design <- model.matrix(~0+d$group)
colnames(design) <- gsub("d\\$group", "", colnames(design))
colnames(design) <- str_replace(colnames(design), " ", ".")
fit <- lmFit(assay(d, "norm"), design)
contrast.matrix <- makeContrasts(mutant-iPD, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
write.table(cbind(rowData(d[[4]]), result), "H:/20230905_150637_Total_Proteome.wide.protein.quant.de.0.7.b.mutant.tsv", sep="\t", row.names = FALSE)
write.table(e, "H:/20230905_150637_Total_Proteome.wide.protein.quant.mutant.control.ipd.raw.tsv", sep="\t", row.names = FALSE)