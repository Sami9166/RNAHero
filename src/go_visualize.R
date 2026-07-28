args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: go_visualize.R GO.csv OUTPUT.png COHORT")

go <- read.csv(args[1], stringsAsFactors = FALSE, check.names = FALSE)
go <- go[order(as.numeric(go$p_value)), , drop = FALSE]
go <- head(go, 20)
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
png(args[2], width = 1600, height = 900, res = 140)
if (!nrow(go)) {
  plot.new()
  text(0.5, 0.5, "No significant GO terms", cex = 1.5)
} else {
  score <- -log10(as.numeric(go$p_value))
  order <- order(score)
  labels <- substr(go$term_name[order], 1, 55)
  colors <- ifelse(go$direction[order] == "up", "#c44e52", "#4c72b0")
  size <- 1 + 3 * as.numeric(go$intersection_size[order]) / max(as.numeric(go$intersection_size))
  par(mar = c(5, 34, 4, 2))
  plot(score[order], seq_along(order), pch = 19, cex = size,
       col = colors, yaxt = "n", ylab = "", xlab = "-log10(GO adjusted p-value)", main = paste("GO enrichment:", args[3]))
  axis(2, at = seq_along(order), labels = labels, las = 1, cex.axis = 0.75)
  legend("bottomright", legend = c("up", "down"), col = c("#c44e52", "#4c72b0"), pch = 19, bty = "n")
}
dev.off()
