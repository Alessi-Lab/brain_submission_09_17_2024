library(QFeatures)
library(stringr)
library(data.table)
a <- read.table("H:/20230905_150637_Total_Proteome.wide.protein.quant.2.tsv", sep = "\t", header=TRUE,comment.char = "", na.strings = c("#VALUE!"))
sample.annotation <-read.table("H:/SAMPLE.INFO.BASICS_v5.txt", sep="\t", header=TRUE, comment.char = "")
a <- a[a["PG.NrOfStrippedSequencesIdentified..Experiment.wide."]>1,]
selectedSamples <- c()
condition <- c()
sample <- c()
rename.cols <- c()
for (i in colnames(a)[3:157]) {
  x <- unlist(strsplit(i, "\\."))
  d <- sample.annotation[sample.annotation["Sample.ID"]==x[4]]
  if (str_starts(x[4], "D")||(str_starts(x[4], "B") && (d[5]=="R1441G"))) {
    sample <- cbind(sample, x[4])
    selectedSamples <- c(selectedSamples, i)
    condition <- cbind(condition, paste(str_replace(d[4], " ", ""), str_replace(d[5], " ","")))
    rename.cols <- cbind(rename.cols, paste(str_replace(d[4], " ", ""), str_replace(d[5], " ",""), x[4], sep = "."))
  }
}
e <- a[c(selectedSamples, "PG.ProteinGroups")]

data <- readQFeatures(e, ecol = selectedSamples, name="PG.ProteinGroups")


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

setnames(normalized.df, old=colnames(normalized.df)[2:(length(selectedSamples)+1)], new=rename.cols)

write.table(normalized.df, "H:/20230905_150637_Total_Proteome.wide.protein.quant.d.normalized.0.7.tsv", sep="\t", row.names = FALSE)

library("limma")
design <- model.matrix(~0+d$group)
colnames(design) <- gsub("d\\$group", "", colnames(design))
colnames(design) <- str_replace(colnames(design), " ", ".")
fit <- lmFit(assay(d, "norm"), design)
contrast.matrix <- makeContrasts(L2PD.R1441G-iPD.., levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
write.table(cbind(rowData(d[[4]]), result), "H:/20230905_150637_Total_Proteome.wide.protein.quant.de.0.7.d.tsv", sep="\t", row.names = FALSE)
write.table(a, "H:/20230905_150637_Total_Proteome.wide.protein.quant.raw.s.tsv", sep="\t", row.names = FALSE)