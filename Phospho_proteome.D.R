
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
library(data.table)
library(limma)
condition <- c()
sample <- c()
rename.cols <- c()
selected.samples <- c()
for (i in colnames(filtered)[1:l]) {
  x <- unlist(strsplit(i, "_"))
  d <- sample.annotation[sample.annotation["Sample.ID"]==x[1]]
  if (d[7] == "D" || d[5] == "R1441G") {
    condition <- cbind(condition, paste(str_replace(d[4], " ", ""), str_replace(d[5], " ","")))
    sample <- cbind(sample, x[1])
    selected.samples <- c(selected.samples, i)
  }
}

for (i in colnames(filtered)[1:l]) {
  x <- unlist(strsplit(i, "_"))
  sample <- cbind(sample, x[1])
  d <- sample.annotation[sample.annotation["Sample.ID"]==x[1]]
  if (d[7] == "D" || d[5] == "R1441G") {
    rename.cols <- cbind(rename.cols, paste(str_replace(d[4], " ", ""), str_replace(d[5], " ",""), x[1], sep = "."))
  }
}

raw.filtered <- subset(a, select= c(selected.samples,"PTM_collapse_key" ))
setnames(raw.filtered, old=colnames(raw.filtered)[1:(length(colnames(raw.filtered))-1)], new=rename.cols)

data <- readQFeatures(subset(filtered, select = c(selected.samples,"PTM_collapse_key", "EG.ModifiedPeptide", "PEP.StrippedSequence", "PG.UniProtIds", "PTM_localization")), ecol = rep(1:length(selected.samples)), name="PTM_collapse_key")
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

normalized.df <- cbind(as.data.frame(rowData(d[["norm"]])), as.data.frame(assay(d, "norm")))

setnames(normalized.df, old=colnames(normalized.df)[6:length(normalized.df)], new=rename.cols)

write.table(normalized.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.D.normalized.tsv", sep="\t", row.names = FALSE)
design <- model.matrix(~0+d$group)
colnames(design) <- gsub("d\\$group", "", colnames(design))
colnames(design) <- str_replace(colnames(design), " ", ".")
fit <- lmFit(assay(d, "norm"), design)
contrast.matrix <- makeContrasts(L2PD.R1441G-Control.., levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
fin.df <- cbind(rowData(d[[4]]), result)
write.table(fin.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.de.0.3.D.L2PD.R1441G-Control.tsv", sep="\t", row.names = FALSE)
contrast.matrix <- makeContrasts(L2NMC.R1441G-Control.., levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
fin.df <- cbind(rowData(d[[4]]), result)
write.table(fin.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.de.0.3.D.L2NMC.R1441G-Control.tsv", sep="\t", row.names = FALSE)
contrast.matrix <- makeContrasts(L2PD.R1441G-L2NMC.R1441G, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
fin.df <- cbind(rowData(d[[4]]), result)
write.table(fin.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.de.0.3.D.L2PD.R1441G-L2NMC.R1441G.tsv", sep="\t", row.names = FALSE)
contrast.matrix <- makeContrasts(L2PD.R1441G-iPD.., levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
fin.df <- cbind(rowData(d[[4]]), result)
write.table(fin.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.de.0.3.D.L2PD.R1441G-iPD.tsv", sep="\t", row.names = FALSE)
contrast.matrix <- makeContrasts(iPD..-Control.., levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, coef=1, number=Inf, sort.by = "none")
fin.df <- cbind(rowData(d[[4]]), result)
write.table(fin.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.de.0.3.D.iPD-Control.tsv", sep="\t", row.names = FALSE)

raw.df <- raw.filtered[raw.filtered$PTM_collapse_key %in% fin.df$PTM_collapse_key,]
raw.df <- subset(raw.df, select = c(colnames(raw.df)[1:(length(colnames(raw.df))-1)], "PTM_collapse_key"))
write.table(raw.df, "H:/20230713_143813_22_phospho_no-norm_Report.output.D.raw.tsv", sep="\t", row.names = FALSE)