
a <- read.table("H:/20230713_143813_22_phospho_no-norm_Report.output.tsv", sep = "\t", header=TRUE,comment.char = "", na.strings = c("#VALUE!", "NA"))
filtered <- a[, which(colMeans(!is.na(a[1:155])) > 0.3)]
l <- length(colnames(filtered))
to.keep <-  c("PTM_collapse_key", "EG.ModifiedPeptide", "PEP.StrippedSequence", "PG.UniProtIds", "PTM_localization")
for (c in to.keep) {
  filtered[c] <- a[c]
}
sample.annotation <-read.table("H:/SAMPLE.INFO.BASICS_v5.txt", sep="\t", header=TRUE, comment.char = "")
library(QFeatures)
library(stringr)
data <- readQFeatures(filtered, ecol = rep(1:l), name="PTM_collapse_key")
condition <- c()
sample <- c()
rename.cols <- c()
for (i in colnames(filtered)[1:l]) {
  x <- unlist(strsplit(i, "_"))
  sample <- cbind(sample, x[1])
  d <- sample.annotation[sample.annotation["Sample.ID"]==x[1]]
  
  if (d[5] != ".") {
    condition <- cbind(condition, d[5])
  } else {
    condition <- cbind(condition, d[4])
  }
}

for (i in colnames(a)[1:155]) {
  x <- unlist(strsplit(i, "_"))
  sample <- cbind(sample, x[1])
  d <- sample.annotation[sample.annotation["Sample.ID"]==x[1]]
  if (d[5] != ".") {
    rename.cols <- cbind(rename.cols, paste(d[5], x[1], sep = "."))
  } else {
    rename.cols <- cbind(rename.cols, paste(d[4], x[1], sep = "."))
  }
}

library(data.table)
setnames(a, old=colnames(a)[1:155], new=rename.cols)

data$group <- condition
data$sample <- sample
colData(data)
rowDataNames(data)
d <- selectRowData(data, to.keep)
d <- zeroIsNA(d, i = seq_along(d))
nNA(d, i = seq_along(d))
d <- filterNA(d, i = seq_along(d), pNA = 0.3)

filtered <- assay(d[[seq_along(d)]])
d <- impute(d, method="knn", i="PTM_collapse_key")
d <- addAssay(d, logTransform(d[[2]]), name="log2")
d <- addAssay(d,
              normalize(d[["log2"]], method = "quantiles.robust"),
              name = "norm")

library("limma")
design <- model.matrix(~0+d$group)
colnames(design) <- gsub("d\\$group", "", colnames(design))
colnames(design) <- str_replace(colnames(design), " ", ".")
fit <- lmFit(assay(d, "norm"), design)
contrast.matrix <- makeContrasts(G2019S-iPD, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
fin.df <- cbind(rowData(d[[4]]), result)
write.table(fin.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.de.0.3G2019S-iPD.tsv", sep="\t", row.names = FALSE)
raw.df <- a[a$PTM_collapse_key %in% fin.df$PTM_collapse_key,]
raw.df <- subset(raw.df, select = c(colnames(raw.df)[1:155], "PTM_collapse_key"))
write.table(raw.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.raw.carriers.tsv", sep="\t", row.names = FALSE)