args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: edger.R COUNTS.(csv|tsv) SAMPLES.(csv|tsv) OUTPUT_DIR")

read_tabular <- function(path) {
  sep <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\\t"
  read.table(path, header = TRUE, sep = sep, check.names = FALSE, comment.char = "", quote = "")
}

counts <- read_tabular(args[1])
samples <- read_tabular(args[2])
if (!all(c("sample_id", "group") %in% names(samples))) stop("sample sheet needs sample_id and group columns")
if (ncol(counts) < 3) stop("count matrix needs a gene column and at least two sample columns")
if (anyDuplicated(counts[[1]])) stop("gene identifiers in the first column must be unique")
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) stop("matrix must contain raw non-negative integer counts")
if (!setequal(colnames(counts), samples$sample_id)) stop("sample IDs must exactly match count-matrix columns")
samples <- samples[match(colnames(counts), samples$sample_id), ]
if (!setequal(unique(samples$group), c("case", "control"))) stop("group values must be case and control")
group <- factor(samples$group, levels = c("control", "case"))
if (min(table(group)) < 2) stop("each group needs at least two samples")

suppressPackageStartupMessages(library(edgeR))
y <- DGEList(counts = counts, group = group)
y <- y[filterByExpr(y, group = group), , keep.lib.sizes = FALSE]
if (nrow(y) == 0) stop("filterByExpr removed every gene")
y <- normLibSizes(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
test <- glmQLFTest(fit, coef = "groupcase")
result <- topTags(test, n = Inf, sort.by = "none")$table
result <- cbind(gene_id = rownames(result), result)
logcpm <- cbind(gene_id = rownames(y), cpm(y, log = TRUE, prior.count = 2))
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(args[3], "edger_results.csv"), row.names = FALSE)
write.csv(logcpm, file.path(args[3], "logcpm.csv"), row.names = FALSE)
