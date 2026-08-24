
.require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("缺少R包：", paste(missing, collapse = ", "))
}

.prepare_output_dir <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) {
    stop("输出目录非空，请指定新目录：", path)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

.ascii_input <- function(path) {
  if (!file.exists(path)) stop("找不到输入文件：", path)
  ext <- tools::file_ext(path)
  tmp <- tempfile(fileext = if (nzchar(ext)) paste0(".", ext) else "")
  if (!file.copy(path, tmp, overwrite = TRUE)) stop("无法建立输入文件临时副本：", path)
  tmp
}

.read_fpkm <- function(path) {
  tmp <- .ascii_input(path); on.exit(unlink(tmp), add = TRUE)
  x <- read.delim(tmp, check.names = FALSE, stringsAsFactors = FALSE)
  sample_cols <- c("FI1", "FI2", "FI3", "LI1", "LI2")
  if (!all(sample_cols %in% names(x))) stop("FPKM必须包含：", paste(sample_cols, collapse = ", "))
  gene_col <- names(x)[1]
  if (anyDuplicated(x[[gene_col]])) stop("FPKM基因ID存在重复：", gene_col)
  mat <- as.matrix(x[, sample_cols, drop = FALSE]); storage.mode(mat) <- "numeric"
  rownames(mat) <- x[[gene_col]]
  if (any(!is.finite(mat))) stop("FPKM表达矩阵含非有限值。")
  mat
}

.module_map <- data.frame(
  ModuleColor = c("red", "yellow", "green", "black", "purple", "lightcyan",
                  "turquoise", "cyan", "brown", "magenta", "pink", "tan",
                  "greenyellow", "midnightblue", "blue", "salmon"),
  DisplayModule = paste0("ME", c("one", "two", "three", "four", "five", "six",
                                   "seven", "eight", "nine", "ten", "eleven", "twelve",
                                   "thirteen", "fourteen", "fifteen", "sixteen")),
  stringsAsFactors = FALSE
)

.save_table_pair <- function(x, stem) {
  write.csv(x, paste0(stem, ".csv"), row.names = FALSE, fileEncoding = "UTF-8")
  openxlsx::write.xlsx(x, paste0(stem, ".xlsx"), overwrite = TRUE)
}

.open_figure <- function(path, type, width, height) {
  if (type == "pdf") pdf(path, width = width, height = height, useDingbats = FALSE)
  else tiff(path, width = width, height = height, units = "in", res = 600, compression = "lzw")
}

.bundle_colors <- function(bundle) {
  cols <- bundle$module_colors
  if (is.null(cols)) cols <- WGCNA::labels2colors(bundle$net$colors)
  as.character(cols)
}

.validate_bundle <- function(bundle) {
  required <- c("net", "filtered_expression", "sample_groups", "gene_ids",
                "dendrogram_order", "module_colors", "module_mapping")
  miss <- setdiff(required, names(bundle))
  if (length(miss)) stop("WGCNA对象缺少字段：", paste(miss, collapse = ", "))
  if (length(bundle$gene_ids) != 22577L) stop("入网基因数不是22,577，而是：", length(bundle$gene_ids))
  observed <- sort(unique(.bundle_colors(bundle)))
  observed <- setdiff(observed, "grey")
  if (length(observed) != 16L) stop("非grey模块数不是16，而是：", length(observed))
  if (!setequal(observed, .module_map$ModuleColor)) stop("模块颜色与固定ME映射不一致。")
  invisible(TRUE)
}

.draw_fig1a <- function(bundle, path, type) {
  expr <- bundle$filtered_expression
  groups <- bundle$sample_groups
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  # image()从下向上绘制行，因此先放L、再放F，最终显示为上F、下L。
  group_means <- rbind(
    L = rowMeans(z[, names(groups)[groups == "L"], drop = FALSE]),
    F = rowMeans(z[, names(groups)[groups == "F"], drop = FALSE])
  )
  ord <- bundle$dendrogram_order
  colors <- .bundle_colors(bundle)[ord]
  dend <- bundle$net$dendrograms[[1]]
  pal <- grDevices::colorRampPalette(c("#67A9CF", "#F7F7F7", "#EF8A62"))(101)
  .open_figure(path, type, 11, 6.8); on.exit(dev.off(), add = TRUE)
  layout(matrix(1:3, ncol = 1), heights = c(4.7, 0.35, 1.15))
  par(mar = c(0.5, 8, 2.5, 2))
  plot(dend, labels = FALSE, hang = 0.04, main = "Gene dendrogram and module colors",
       ylab = "Height", xlab = "", sub = "")
  par(mar = c(0, 8, 0, 2))
  module_levels <- unique(colors)
  image(matrix(match(colors, module_levels), ncol = 1), col = module_levels,
        axes = FALSE, xlab = "", ylab = "")
  mtext("Module colors", side = 2, line = 1.8, las = 1, cex = 0.75)
  par(mar = c(2.2, 8, 0, 2))
  image(t(group_means[, ord, drop = FALSE]), col = pal, zlim = c(-2, 2),
        axes = FALSE, xlab = "", ylab = "")
  axis(2, at = c(0, 1), labels = c("L", "F"), las = 1, tick = FALSE)
  # 
  par(fig = c(0.022, 0.039, 0.074, 0.172), new = TRUE, mar = c(0, 0, 0, 0))
  image(1, seq(-2, 2, length.out = length(pal)),
        matrix(seq_along(pal), nrow = 1), col = pal, axes = FALSE,
        xlab = "", ylab = "")
  axis(4, at = c(-2, -1, 0, 1, 2), labels = c(-2, -1, 0, 1, 2),
       las = 1, tick = FALSE, cex.axis = 0.45)
  mtext("Z-score", side = 3, line = 0.2, cex = 0.52)
}

# Stage 01 — Fig.1A
stage_01_fig1a <- function(fpkm_file, output_dir, existing_net_file = NULL,
                           power = 12, min_module_size = 30, merge_cut_height = 0.25) {
  .require_packages(c("WGCNA", "openxlsx"))
  out <- .prepare_output_dir(output_dir)
  fpkm <- .read_fpkm(fpkm_file)
  dat_expr <- as.data.frame(t(fpkm))
  qc <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  if (!qc$allOK) dat_expr <- dat_expr[qc$goodSamples, qc$goodGenes, drop = FALSE]
  groups <- setNames(ifelse(grepl("^FI", rownames(dat_expr)), "F",
                            ifelse(grepl("^LI", rownames(dat_expr)), "L", NA_character_)),
                     rownames(dat_expr))
  if (anyNA(groups) || !identical(sort(unique(groups)), c("F", "L"))) {
    stop("样本名必须以FI或LI开头。")
  }
  if (is.null(existing_net_file)) {
    net <- WGCNA::blockwiseModules(
      dat_expr, power = power, maxBlockSize = 25000, TOMType = "unsigned",
      minModuleSize = min_module_size, reassignThreshold = 0,
      mergeCutHeight = merge_cut_height, numericLabels = TRUE,
      pamRespectsDendro = FALSE, saveTOMs = FALSE, corType = "bicor",
      maxPOutliers = 0.05, verbose = 3
    )
  } else {
    tmp <- .ascii_input(existing_net_file); on.exit(unlink(tmp), add = TRUE)
    net <- readRDS(tmp)
  }
  if (length(net$colors) != ncol(dat_expr)) {
    stop("WGCNA网络对象与过滤后FPKM的基因数不一致。")
  }
  colors <- WGCNA::labels2colors(net$colors)
  names(colors) <- colnames(dat_expr)
  ord <- net$dendrograms[[1]]$order
  bundle <- list(
    net = net,
    filtered_expression = t(as.matrix(dat_expr)),
    sample_groups = groups,
    gene_ids = colnames(dat_expr),
    dendrogram_order = ord,
    dendrogram_gene_ids = colnames(dat_expr)[ord],
    module_colors = colors,
    module_mapping = .module_map,
    parameters = list(power = power, minModuleSize = min_module_size,
                      mergeCutHeight = merge_cut_height)
  )
  .validate_bundle(bundle)
  saveRDS(bundle, file.path(out, "WGCNA_net_object.rds"), compress = "xz")
  gene_modules <- data.frame(GeneID = bundle$gene_ids, ModuleColor = colors,
                             DisplayModule = .module_map$DisplayModule[match(colors, .module_map$ModuleColor)])
  .save_table_pair(gene_modules, file.path(out, "WGCNA_gene_module_assignment"))
  for (type in c("pdf", "tiff")) {
    .draw_fig1a(bundle, file.path(out, paste0("Fig1A.", type)), type)
  }
  invisible(bundle)
}

# Stage 02 — Fig.1B
stage_02_fig1b <- function(wgcna_rds, output_dir) {
  .require_packages(c("WGCNA", "openxlsx"))
  out <- .prepare_output_dir(output_dir)
  tmp <- .ascii_input(wgcna_rds); on.exit(unlink(tmp), add = TRUE)
  bundle <- readRDS(tmp); .validate_bundle(bundle)
  colors <- as.character(WGCNA::labels2colors(bundle$net$colors))
  counts <- as.data.frame(table(colors), stringsAsFactors = FALSE)
  names(counts) <- c("ModuleColor", "GeneCount")
  counts <- merge(.module_map, counts, by = "ModuleColor", all.x = TRUE, sort = FALSE)
  counts <- counts[match(.module_map$ModuleColor, counts$ModuleColor), ]
  if (anyNA(counts$GeneCount) || sum(counts$GeneCount) != 22577L) stop("模块计数校验失败。")
  # 
  counts <- counts[order(-counts$GeneCount), , drop = FALSE]
  .save_table_pair(counts, file.path(out, "Fig1B_module_gene_counts"))
  draw <- function(path, type) {
    .open_figure(path, type, 8, 5.6); on.exit(dev.off(), add = TRUE)
    par(mar = c(8, 4.5, 1.5, 1))
    barplot(counts$GeneCount, names.arg = counts$DisplayModule, col = counts$ModuleColor,
            border = NA, las = 2, ylab = "Number of genes", cex.names = 0.8)
  }
  for (type in c("pdf", "tiff")) draw(file.path(out, paste0("Fig1B.", type)), type)
  invisible(counts)
}

.me_color <- function(x) {
  raw <- sub("^ME", "", x)
  num <- suppressWarnings(as.numeric(raw))
  ifelse(is.na(num), raw, WGCNA::labels2colors(num))
}

.draw_correlation_heatmap <- function(cor_mat, p_mat, labels, path, type) {
  pal <- grDevices::colorRampPalette(c("#3E9BE8", "white", "#FF3B1F"))(101)
  .open_figure(path, type, 5.5, 8.3); on.exit(dev.off(), add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(0.055, 2.5, 0.25))
  module_cols <- .module_map$ModuleColor[match(labels, .module_map$DisplayModule)]
  par(mar = c(4, 0, 1.5, 0))
  image(1, seq_along(labels), matrix(seq_along(labels), nrow = 1),
        col = rev(module_cols), axes = FALSE, xlab = "", ylab = "")
  par(mar = c(4, 6.7, 1.5, 0.3))
  image(seq_len(ncol(cor_mat)), seq_len(nrow(cor_mat)), t(cor_mat[nrow(cor_mat):1, , drop = FALSE]),
        col = pal, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_len(ncol(cor_mat)), labels = colnames(cor_mat), tick = FALSE)
  axis(2, at = seq_len(nrow(cor_mat)), labels = rev(labels), las = 1, tick = FALSE)
  for (i in seq_len(nrow(cor_mat))) for (j in seq_len(ncol(cor_mat))) {
    ii <- nrow(cor_mat) - i + 1
    text(j, ii, sprintf("%.2f\n(%.2g)", cor_mat[i, j], p_mat[i, j]), cex = 0.72)
  }
  box()
  par(mar = c(4, 0.2, 1.5, 2.5))
  image(1, seq(-1, 1, length.out = length(pal)),
        matrix(seq_along(pal), nrow = 1), col = pal, axes = FALSE,
        xlab = "", ylab = "")
  axis(4, at = c(-1, -0.5, 0, 0.5, 1), las = 1, tick = FALSE, cex.axis = 0.72)
}

# Stage 03 — Fig.1C
stage_03_fig1c <- function(wgcna_rds, output_dir) {
  .require_packages(c("WGCNA", "openxlsx"))
  out <- .prepare_output_dir(output_dir)
  tmp <- .ascii_input(wgcna_rds); on.exit(unlink(tmp), add = TRUE)
  bundle <- readRDS(tmp); .validate_bundle(bundle)
  mes <- WGCNA::orderMEs(bundle$net$MEs)
  me_colors <- .me_color(colnames(mes))
  display <- .module_map$DisplayModule[match(me_colors, .module_map$ModuleColor)]
  keep <- !is.na(display)
  mes <- mes[, keep, drop = FALSE]; display <- display[keep]
  ord <- match(.module_map$DisplayModule, display)
  if (anyNA(ord)) stop("网络对象的模块特征向量无法完整映射到MEone–MEsixteen。")
  mes <- mes[, ord, drop = FALSE]
  traits <- cbind(F = as.integer(bundle$sample_groups[rownames(mes)] == "F"),
                  L = as.integer(bundle$sample_groups[rownames(mes)] == "L"))
  cors <- stats::cor(mes, traits, use = "pairwise.complete.obs")
  ps <- WGCNA::corPvalueStudent(cors, nrow(mes))
  tab <- data.frame(DisplayModule = .module_map$DisplayModule,
                    ModuleColor = .module_map$ModuleColor,
                    Cor_F = cors[, "F"], P_F = ps[, "F"],
                    Cor_L = cors[, "L"], P_L = ps[, "L"])
  .save_table_pair(tab, file.path(out, "Fig1C_module_trait_correlations"))
  for (type in c("pdf", "tiff")) {
    .draw_correlation_heatmap(cors, ps, .module_map$DisplayModule,
                              file.path(out, paste0("Fig1C.", type)), type)
  }
  invisible(tab)
}
