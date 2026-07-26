args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: visualize.R LOGCPM.csv SAMPLES.csv GENES.txt OUTPUT_DIR")

logcpm <- read.csv(args[1], check.names = FALSE)
samples <- read.csv(args[2], stringsAsFactors = FALSE)
genes <- unique(readLines(args[3], warn = FALSE))
if (!all(c("sample_id", "group") %in% names(samples))) stop("sample sheet needs sample_id and group columns")
rownames(logcpm) <- logcpm$gene_id
values <- as.matrix(logcpm[intersect(genes, rownames(logcpm)), samples$sample_id, drop = FALSE])
storage.mode(values) <- "numeric"
values <- values[apply(values, 1, sd) > 0, , drop = FALSE]
if (nrow(values) < 2) stop("fewer than two variable selected genes are available")

dir.create(args[4], recursive = TRUE, showWarnings = FALSE)
colors <- ifelse(samples$group == "case", "#c44e52", "#4c72b0")
names(colors) <- samples$sample_id
scaled <- t(scale(t(values)))

png(file.path(args[4], "top5_heatmap.png"), width = 1000, height = 700, res = 130)
heatmap(scaled, Colv = NA, scale = "none", ColSideColors = colors[colnames(scaled)], margins = c(8, 10), main = "Top internal logFC genes")
legend("topright", legend = c("control", "case"), fill = c("#4c72b0", "#c44e52"), bty = "n")
dev.off()

pca <- prcomp(t(values), scale. = TRUE)
x <- pca$x[, 1]
y <- if (ncol(pca$x) >= 2) pca$x[, 2] else rep(0, length(x))
variance <- summary(pca)$importance[2, ] * 100
png(file.path(args[4], "top5_pca.png"), width = 900, height = 700, res = 130)
plot(x, y, col = colors[names(x)], pch = 19, xlab = sprintf("PC1 (%.1f%%)", variance[1]), ylab = sprintf("PC2 (%.1f%%)", ifelse(length(variance) >= 2, variance[2], 0)), main = "PCA: top internal logFC genes")
text(x, y, labels = names(x), pos = 3, cex = 0.7)
legend("topright", legend = c("control", "case"), col = c("#4c72b0", "#c44e52"), pch = 19, bty = "n")
dev.off()
