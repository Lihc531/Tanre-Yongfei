# Canonical recovered two-cohort CAP preprocessing workflow; paths refactored for public repository.
############################################################
# Figure2_FINAL_STABLE_FROM_TWO_SCRIPTS.R
# GSE103119: GEO matrix / microarray -> limma
# GSE196399: count matrix -> edgeR + voom + limma
############################################################

rm(list = ls())

## =========================
## 0) Packages
## =========================
pkgs_cran <- c("ggplot2", "dplyr", "patchwork", "ggVennDiagram", "umap", "ggrepel", "ggforce", "ggalluvial")
pkgs_bioc <- c("GEOquery", "Biobase", "limma", "edgeR",
               "AnnotationDbi", "illuminaHumanv4.db", "org.Hs.eg.db")

missing_cran <- pkgs_cran[!vapply(pkgs_cran, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- pkgs_bioc[!vapply(pkgs_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) + length(missing_bioc) > 0) {
  stop("Missing R packages. See environment/R_PACKAGES.md. CRAN: ",
       paste(missing_cran, collapse=", "), "; Bioconductor: ",
       paste(missing_bioc, collapse=", "))
}

library(GEOquery)
library(Biobase)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(illuminaHumanv4.db)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggVennDiagram)
library(umap)
library(ggrepel)
library(ggforce)
library(ggalluvial)

## =========================
## 1) Parameters
## =========================

args <- commandArgs(trailingOnly = TRUE)
gse103119_file <- if (length(args) >= 1) args[1] else "data/restricted_inputs/GSE103119_series_matrix.txt.gz"
gse196399_file <- if (length(args) >= 2) args[2] else "data/restricted_inputs/GSE196399_count_matrix.csv.gz"
out_dir <- if (length(args) >= 3) args[3] else "results/generated/transcriptomics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff  <- 0.05
logfc_cutoff <- 1

col_cap     <- "#E64B35"
col_control <- "#4DBBD5"
col_ns      <- "grey75"

theme_pub <- function() {
  theme_classic(base_size = 13) +
    theme(
      axis.line = element_line(color = "black"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

## =========================
## 2) Utility functions
## =========================

make_valid_group <- function(group) {
  group <- factor(group)
  group <- droplevels(group)
  levels(group) <- make.names(levels(group))
  group
}

check_two_group <- function(expr, group, name = "dataset") {
  if (length(group) != ncol(expr)) {
    stop(name, ": group长度与表达矩阵列数不一致。group = ",
         length(group), ", ncol(expr) = ", ncol(expr))
  }
  if (nlevels(group) != 2) {
    stop(name, ": 当前不是两组比较。检测到分组为：",
         paste(levels(group), collapse = ", "))
  }
  message(name, " group table:")
  print(table(group))
}

safe_log2_if_needed <- function(expr) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  
  expr[!is.finite(expr)] <- NA
  
  qx <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
  print(qx)
  
  # 关键修复：
  # 如果表达矩阵含负值，通常说明已经是log2/normalized array，不再log2
  if (qx[1] < 0) {
    message("Expression contains negative values. Treat as already normalized/log2-like. Skip log2.")
  } else if (qx[5] > 100 || (qx[5] - qx[1] > 50)) {
    message("Expression seems not log2-transformed. Applying log2(x + 1).")
    expr <- log2(expr + 1)
  } else {
    message("Expression seems already log2-like. Skip log2 transform.")
  }
  
  expr[!is.finite(expr)] <- NA
  
  miss_rate <- rowMeans(is.na(expr))
  expr <- expr[miss_rate <= 0.2, , drop = FALSE]
  
  row_mean <- rowMeans(expr, na.rm = TRUE)
  idx <- which(is.na(expr), arr.ind = TRUE)
  
  if (nrow(idx) > 0) {
    expr[idx] <- row_mean[idx[, 1]]
  }
  
  expr <- expr[apply(expr, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
  
  return(expr)
}

collapse_by_symbol <- function(res, symbol_col = "SYMBOL") {
  res <- res[!is.na(res[[symbol_col]]) & res[[symbol_col]] != "", ]
  res <- res[order(res$adj.P.Val, -abs(res$logFC)), ]
  res_gene <- res[!duplicated(res[[symbol_col]]), ]
  res_gene$gene <- res_gene[[symbol_col]]
  res_gene
}

## =========================
## 3) Load GSE103119: matrix / microarray
## =========================

load_gse103119_from_original_file <- function(file) {
  
  con <- if (grepl("\\.gz$", file)) gzfile(file, open = "rt") else base::file(file, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con)
  
  ## 1. 解析 sample title
  title_line <- lines[grep("^!Sample_title", lines)]
  title_vec <- strsplit(title_line, "\t")[[1]][-1]
  
  title_vec <- gsub('^"|"$', "", title_vec)
  
  group <- rep(NA_character_, length(title_vec))
  
  group[grepl("-Healthy Control\\s", title_vec, ignore.case = TRUE)] <- "Control"
  group[grepl("-Pneumonia\\s", title_vec, ignore.case = TRUE)] <- "CAP"
  
  message("Group table parsed from original !Sample_title:")
  print(table(group, useNA = "ifany"))
  
  if (any(is.na(group))) {
    print(title_vec[is.na(group)])
    stop("仍有样本无法识别分组，请检查 !Sample_title 格式。")
  }
  
  group <- factor(group, levels = c("Control", "CAP"))
  
  ## 2. 解析表达矩阵
  start <- grep("^!series_matrix_table_begin", lines)
  end   <- grep("^!series_matrix_table_end", lines)
  
  expr_raw <- read.table(
    text = lines[(start + 1):(end - 1)],
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )
  
  rownames(expr_raw) <- expr_raw[, 1]
  expr <- expr_raw[, -1]
  
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  
  ## 3. 检查维度
  message("Expression dim:")
  print(dim(expr))
  
  if (ncol(expr) != length(group)) {
    stop("表达矩阵列数与group长度不一致：ncol(expr) = ",
         ncol(expr), ", length(group) = ", length(group))
  }
  
  ## 4. 检查原始分布
  message("Raw expression quantiles:")
  print(quantile(expr, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
  
  ## 5. 关键修复：不能直接 log2(expr + 1)
  expr <- log2(pmax(expr, 1))
  
  ## 6. 芯片数据推荐做分位数标准化
  expr <- limma::normalizeBetweenArrays(expr, method = "quantile")
  
  ## 7. 去掉异常值/无变化探针
  expr[!is.finite(expr)] <- NA
  
  miss_rate <- rowMeans(is.na(expr))
  expr <- expr[miss_rate <= 0.2, , drop = FALSE]
  
  row_mean <- rowMeans(expr, na.rm = TRUE)
  idx <- which(is.na(expr), arr.ind = TRUE)
  
  if (nrow(idx) > 0) {
    expr[idx] <- row_mean[idx[, 1]]
  }
  
  expr <- expr[apply(expr, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
  
  message("Final expression dim:")
  print(dim(expr))
  
  message("Final group table:")
  print(table(group))
  
  return(list(expr = expr, group = group, title = title_vec))
}

## =========================
## 4) Load GSE196399: count matrix
## =========================

load_gse196399_count <- function(file) {
  
  if (grepl("\\.gz$", file)) {
    counts <- read.csv(gzfile(file), row.names = 1, check.names = FALSE)
  } else {
    counts <- read.csv(file, row.names = 1, check.names = FALSE)
  }
  
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  
  counts[!is.finite(counts)] <- 0
  counts[counts < 0] <- 0
  
  message("GSE196399 count dim:")
  print(dim(counts))
  
  ## 自动分组：NC = Control, SP = CAP
  group <- ifelse(grepl("^NC", colnames(counts), ignore.case = TRUE),
                  "Control",
                  ifelse(grepl("^SP", colnames(counts), ignore.case = TRUE),
                         "CAP", NA))
  
  if (any(is.na(group))) {
    stop("GSE196399: 有样本无法通过列名NC/SP识别分组，请检查colnames(counts)。")
  }
  
  group <- factor(group, levels = c("Control", "CAP"))
  group <- make_valid_group(group)
  
  check_two_group(counts, group, "GSE196399")
  
  list(counts = counts, group = group)
}

## =========================
## 5) DEG functions
## =========================

run_limma_matrix <- function(expr, group, dataset_name = "matrix_dataset") {
  
  group <- make_valid_group(group)
  check_two_group(expr, group, dataset_name)
  
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  fit <- lmFit(expr, design)
  
  contrast_matrix <- makeContrasts(
    contrasts = "CAP-Control",
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  res <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
  res$feature_id <- rownames(res)
  
  res
}

run_voom_count <- function(counts, group, dataset_name = "count_dataset") {
  
  group <- make_valid_group(group)
  check_two_group(counts, group, dataset_name)
  
  dge <- DGEList(counts = counts, group = group)
  
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  dge <- calcNormFactors(dge, method = "TMM")
  
  v <- voom(dge, design, plot = FALSE)
  
  fit <- lmFit(v, design)
  
  contrast_matrix <- makeContrasts(
    contrasts = "CAP-Control",
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  res <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
  res$feature_id <- rownames(res)
  
  list(res = res, voom_expr = v$E)
}

## =========================
## 6) Annotation
## =========================

annotate_gse103119_probe <- function(res) {
  
  ann <- AnnotationDbi::select(
    illuminaHumanv4.db,
    keys = res$feature_id,
    keytype = "PROBEID",
    columns = c("SYMBOL", "GENENAME")
  )
  
  res2 <- merge(res, ann, by.x = "feature_id", by.y = "PROBEID", all.x = TRUE)
  
  gene_level <- collapse_by_symbol(res2, symbol_col = "SYMBOL")
  
  gene_level
}

annotate_gse196399_ensembl <- function(res) {
  
  ensembl <- gsub("\\..*$", "", res$feature_id)
  
  ann <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ensembl,
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENENAME")
  )
  
  map_df <- data.frame(
    feature_id = res$feature_id,
    ENSEMBL = ensembl,
    stringsAsFactors = FALSE
  )
  
  res2 <- merge(map_df, res, by = "feature_id", all.x = TRUE)
  res2 <- merge(res2, ann, by = "ENSEMBL", all.x = TRUE)
  
  gene_level <- collapse_by_symbol(res2, symbol_col = "SYMBOL")
  
  gene_level
}

## =========================
## 7) Plot functions
## =========================

plot_volcano <- function(res, title) {
  
  # 仅对“优美火山图”代码做最低限度适配：
  # logfc -> logFC；pval -> adj.P.Val；gene 保留当前基因列
  deg <- res %>%
    dplyr::filter(
      is.finite(logFC),
      is.finite(adj.P.Val),
      adj.P.Val > 0
    )
  
  logfc_low  <- -logfc_cutoff
  logfc_high <-  logfc_cutoff
  p <- padj_cutoff
  num <- 8
  
  # “优美火山图”原配色
  bullcol <- c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  
  # 提取核心基因
  core_gene_data <- deg %>%
    arrange(adj.P.Val) %>%
    slice_head(n = num) %>%
    mutate(pvalue_log = -log10(adj.P.Val))
  
  # 绘制渐变火山图
  ggplot(deg, aes(logFC, -log10(adj.P.Val))) +
    geom_hline(yintercept = -log10(p), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(logfc_low, logfc_high), linetype = "dashed", color = "black") +
    geom_point(aes(size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
    scale_color_gradientn(colors = bullcol) +
    scale_size_continuous(range = c(1, 3)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "right", legend.justification = c(0, 1)) +
    guides(color = guide_colorbar(title = "-Log10(adj.P.Val)", barheight = 8, barwidth = 1), size = "none") +
    xlab("Log2FC") + ylab("-Log10(adj.P.Val)") +
    labs(title = title) +
    ggrepel::geom_text_repel(
      data = core_gene_data,
      aes(label = gene, color = pvalue_log),
      size = 3.5,
      max.overlaps = 15,
      force = 10,
      box.padding = grid::unit(0.5, "lines"),
      nudge_x = 0.5,
      nudge_y = 1
    )
}

plot_umap_safe <- function(expr, group, title, top_n = 5000) {
  
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  
  expr[!is.finite(expr)] <- NA
  
  keep <- rowSums(is.na(expr)) == 0
  expr <- expr[keep, , drop = FALSE]
  
  vars <- apply(expr, 1, var)
  vars <- sort(vars, decreasing = TRUE)
  
  top_features <- names(vars)[seq_len(min(top_n, length(vars)))]
  expr_use <- expr[top_features, , drop = FALSE]
  
  set.seed(123)
  um <- umap(t(expr_use))
  
  df <- data.frame(
    UMAP1 = um$layout[, 1],
    UMAP2 = um$layout[, 2],
    group = group
  )
  
  ggplot(df, aes(UMAP1, UMAP2, color = group)) +
    geom_point(size = 3, alpha = 0.85) +
    stat_ellipse(level = 0.95, linewidth = 0.6) +
    scale_color_manual(values = c("Control" = col_control,
                                  "CAP" = col_cap)) +
    labs(title = title, x = "UMAP1", y = "UMAP2") +
    theme_pub()
}



############################################################
# B: logFC correlation plot
############################################################

plot_logfc_correlation <- function(df_cor) {
  
  cor_s <- suppressWarnings(
    cor.test(
      df_cor$GSE103119_logFC,
      df_cor$GSE196399_logFC,
      method = "spearman"
    )
  )
  
  rho_label <- paste0(
    "Spearman rho = ",
    round(unname(cor_s$estimate), 3),
    "\nP = ",
    formatC(cor_s$p.value, format = "e", digits = 2)
  )
  
  ggplot(df_cor, aes(GSE103119_logFC, GSE196399_logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey45", linewidth = 0.5) +
    geom_point(alpha = 0.55, color = col_control, size = 1.35) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.8) +
    annotate(
      "text",
      x = min(df_cor$GSE103119_logFC, na.rm = TRUE),
      y = max(df_cor$GSE196399_logFC, na.rm = TRUE),
      label = rho_label,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = "Cross-cohort logFC consistency",
      x = "GSE103119 logFC",
      y = "GSE196399 logFC"
    ) +
    theme_pub()
}

############################################################
# C: custom vertical Venn plot
############################################################

make_lens_df <- function(r = 1, y_top = 0.45, y_bottom = -0.45, n = 200) {
  
  d <- abs(y_top - y_bottom)
  
  if (d >= 2 * r) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  
  y_mid <- (y_top + y_bottom) / 2
  x_int <- sqrt(r^2 - (d / 2)^2)
  
  a_top_left  <- atan2(y_mid - y_top, -x_int)
  a_top_right <- atan2(y_mid - y_top,  x_int)
  
  a_bot_right <- atan2(y_mid - y_bottom,  x_int)
  a_bot_left  <- atan2(y_mid - y_bottom, -x_int)
  
  top_angles <- seq(a_top_left, a_top_right, length.out = n)
  bot_angles <- seq(a_bot_right, a_bot_left, length.out = n)
  
  top_arc <- data.frame(
    x = r * cos(top_angles),
    y = y_top + r * sin(top_angles)
  )
  
  bot_arc <- data.frame(
    x = r * cos(bot_angles),
    y = y_bottom + r * sin(bot_angles)
  )
  
  rbind(top_arc, bot_arc)
}

plot_custom_venn <- function(set1,
                             set2,
                             title,
                             top_fill,
                             bottom_fill,
                             overlap_fill) {
  
  set1 <- unique(set1)
  set2 <- unique(set2)
  
  n1 <- length(set1)
  n2 <- length(set2)
  n_shared <- length(intersect(set1, set2))
  
  pct_smaller <- ifelse(
    min(n1, n2) > 0,
    100 * n_shared / min(n1, n2),
    NA
  )
  
  overlap_label <- paste0(
    n_shared,
    " shared\n",
    round(pct_smaller, 1),
    "% of smaller set"
  )
  
  circle_df <- data.frame(
    x0 = c(0, 0),
    y0 = c(0.45, -0.45),
    r = c(1, 1),
    dataset = c("GSE103119", "GSE196399")
  )
  
  lens_df <- make_lens_df(r = 1, y_top = 0.45, y_bottom = -0.45)
  
  ggplot() +
    ggforce::geom_circle(
      data = circle_df[1, ],
      aes(x0 = x0, y0 = y0, r = r),
      fill = top_fill,
      color = "black",
      linewidth = 0.7,
      alpha = 0.65
    ) +
    ggforce::geom_circle(
      data = circle_df[2, ],
      aes(x0 = x0, y0 = y0, r = r),
      fill = bottom_fill,
      color = "black",
      linewidth = 0.7,
      alpha = 0.65
    ) +
    geom_polygon(
      data = lens_df,
      aes(x, y),
      fill = overlap_fill,
      color = "black",
      linewidth = 0.4,
      alpha = 0.98
    ) +
    annotate("text", x = 0, y = 0.98, label = "GSE103119", size = 4.1, fontface = "bold") +
    annotate("text", x = 0, y = -0.98, label = "GSE196399", size = 4.1, fontface = "bold") +
    annotate("text", x = 0, y = 0, label = overlap_label, size = 4.0, fontface = "bold") +
    coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.65, 1.65), expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin = margin(20, 20, 20, 20)
    )
}

############################################################
# D: Sankey / alluvial plot for direction concordance
############################################################

plot_direction_sankey <- function(deg103_gene,
                                  deg196_gene,
                                  use_deg_only = TRUE) {
  
  common_genes <- intersect(deg103_gene$gene, deg196_gene$gene)
  
  df <- data.frame(
    gene = common_genes,
    logFC_103 = deg103_gene$logFC[match(common_genes, deg103_gene$gene)],
    logFC_196 = deg196_gene$logFC[match(common_genes, deg196_gene$gene)],
    padj_103 = deg103_gene$adj.P.Val[match(common_genes, deg103_gene$gene)],
    padj_196 = deg196_gene$adj.P.Val[match(common_genes, deg196_gene$gene)]
  )
  
  df <- df[is.finite(df$logFC_103) & is.finite(df$logFC_196), ]
  df <- df[df$logFC_103 != 0 & df$logFC_196 != 0, ]
  
  if (use_deg_only) {
    df <- df[
      (df$padj_103 < padj_cutoff & abs(df$logFC_103) >= logfc_cutoff) |
        (df$padj_196 < padj_cutoff & abs(df$logFC_196) >= logfc_cutoff),
    ]
  }
  df$GSE103119 <- ifelse(df$logFC_103 > 0, "Up", "Down")
  df$GSE196399 <- ifelse(df$logFC_196 > 0, "Up", "Down")
  
  df$GSE103119 <- factor(df$GSE103119, levels = c("Up", "Down"))
  df$GSE196399 <- factor(df$GSE196399, levels = c("Up", "Down"))
  
  flow_df <- df %>%
    count(GSE103119, GSE196399, name = "n") %>%
    mutate(
      transition = paste0(GSE103119, " → ", GSE196399),
      concordance = ifelse(GSE103119 == GSE196399, "Concordant", "Discordant"),
      pct = n / sum(n)
    )
  
  concordant_pct <- 100 * sum(flow_df$pct[flow_df$concordance == "Concordant"])
  discordant_pct <- 100 * sum(flow_df$pct[flow_df$concordance == "Discordant"])
  
  flow_df$GSE103119 <- factor(flow_df$GSE103119, levels = c("Up", "Down"))
  flow_df$GSE196399 <- factor(flow_df$GSE196399, levels = c("Up", "Down"))
  
  ggplot(
    flow_df,
    aes(axis1 = GSE103119, axis2 = GSE196399, y = pct)
  ) +
    geom_alluvium(
      aes(fill = transition),
      width = 0.18,
      alpha = 0.82,
      knot.pos = 0.45
    ) +
    geom_stratum(
      width = 0.22,
      fill = "grey96",
      color = "grey35",
      linewidth = 0.5
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 4.2,
      fontface = "bold"
    ) +
    scale_x_discrete(
      limits = c("GSE103119", "GSE196399"),
      expand = c(0.16, 0.08)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    scale_fill_manual(
      values = c(
        "Up → Up" = col_cap,
        "Up → Down" = "grey70",
        "Down → Up" = "grey55",
        "Down → Down" = col_control
      ),
      breaks = c(
        "Up → Up",
        "Up → Down",
        "Down → Up",
        "Down → Down"
      )
    ) +
    labs(
      title = "Direction concordance of DEGs",
      subtitle = paste0(
        "Concordant = ",
        round(concordant_pct, 1),
        "%; Discordant = ",
        round(discordant_pct, 1),
        "%"
      ),
      x = NULL,
      y = "Proportion of genes"
    ) +
    theme_pub() +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      plot.margin = margin(25, 20, 25, 20)
    )
}


## =========================
## 8) Run analysis
## =========================

g103 <- load_gse103119_from_original_file(gse103119_file)
g196 <- load_gse196399_count(gse196399_file)

deg103_probe <- run_limma_matrix(
  expr = g103$expr,
  group = g103$group,
  dataset_name = "GSE103119"
)

tmp196 <- run_voom_count(
  counts = g196$counts,
  group = g196$group,
  dataset_name = "GSE196399"
)

deg196_raw <- tmp196$res
expr196_voom <- tmp196$voom_expr

deg103_gene <- annotate_gse103119_probe(deg103_probe)
deg196_gene <- annotate_gse196399_ensembl(deg196_raw)

write.csv(deg103_gene, file.path(out_dir, "GSE103119_gene_level_all.csv"), row.names = FALSE)
write.csv(subset(deg103_gene, adj.P.Val < padj_cutoff & abs(logFC) >= logfc_cutoff),
          file.path(out_dir, "GSE103119_gene_level_DEG.csv"), row.names = FALSE)

write.csv(deg196_gene, file.path(out_dir, "GSE196399_gene_level_all.csv"), row.names = FALSE)
write.csv(subset(deg196_gene, adj.P.Val < padj_cutoff & abs(logFC) >= logfc_cutoff),
          file.path(out_dir, "GSE196399_gene_level_DEG.csv"), row.names = FALSE)

## =========================
## 9) Figure2 panels: revised A | B / C | D layout
## =========================

make_panel_title <- function(label, title, size = 15) {
  
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.5,
      label = paste(label, title),
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = size / ggplot2::.pt
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 2, 0)
    )
}


add_panel_title <- function(p, label, title) {
  
  make_panel_title(label, title) / p +
    plot_layout(heights = c(0.08, 1))
}

## A. Volcano plots: DEG screening in two cohorts
panel_title_theme <- theme(
  plot.title = element_text(
    hjust = 0,
    face = "bold",
    size = 15,
    margin = margin(b = 8)
  )
)

## A. Volcano plots: DEG screening in two cohorts
A1 <- plot_volcano(deg103_gene, "GSE103119")
A2 <- plot_volcano(deg196_gene, "GSE196399")

A_panel <- add_panel_title(
  A1 | A2,
  "A",
  "Differential expression screening"
)



## B. Cross-cohort logFC correlation
common_genes <- intersect(deg103_gene$gene, deg196_gene$gene)

df_cor <- data.frame(
  gene = common_genes,
  GSE103119_logFC = deg103_gene$logFC[match(common_genes, deg103_gene$gene)],
  GSE196399_logFC = deg196_gene$logFC[match(common_genes, deg196_gene$gene)]
)

B_panel <- add_panel_title(
  plot_logfc_correlation(df_cor) + labs(title = NULL),
  "B",
  "Cross-cohort logFC consistency"
)


## C. Venn plots: upregulated and downregulated shared genes
up103 <- deg103_gene$gene[
  deg103_gene$adj.P.Val < padj_cutoff &
    deg103_gene$logFC >= logfc_cutoff
]

up196 <- deg196_gene$gene[
  deg196_gene$adj.P.Val < padj_cutoff &
    deg196_gene$logFC >= logfc_cutoff
]

down103 <- deg103_gene$gene[
  deg103_gene$adj.P.Val < padj_cutoff &
    deg103_gene$logFC <= -logfc_cutoff
]

down196 <- deg196_gene$gene[
  deg196_gene$adj.P.Val < padj_cutoff &
    deg196_gene$logFC <= -logfc_cutoff
]

C_up <- plot_custom_venn(
  set1 = up103,
  set2 = up196,
  title = "Upregulated genes",
  top_fill = "#B2182B",
  bottom_fill = "#EF8A62",
  overlap_fill = "#FDE0DD"
)

C_down <- plot_custom_venn(
  set1 = down103,
  set2 = down196,
  title = "Downregulated genes",
  top_fill = "#2166AC",
  bottom_fill = "#67A9CF",
  overlap_fill = "#E0F3F8"
)

C_panel <- add_panel_title(
  C_up | C_down,
  "C",
  "Overlap of CAP-regulated genes"
)


## D. Sankey plot: direction concordance
D_panel <- plot_direction_sankey(
  deg103_gene = deg103_gene,
  deg196_gene = deg196_gene,
  use_deg_only = TRUE
) +
  labs(
    title = "D Direction concordance of DEGs",
    subtitle = NULL
  ) +
  panel_title_theme


## =========================
## 10) Final Figure2: A | B / C | D
## =========================

Figure2 <- (A_panel | B_panel) /
  (C_panel | D_panel) +
  plot_layout(widths = c(1.35, 1), heights = c(1, 1)) &
  theme(plot.margin = margin(20, 20, 20, 20))

Figure2

ggsave(file.path(out_dir, "Figure2_revised_grouped_ABCD.pdf"), Figure2, width = 18, height = 14)
ggsave(file.path(out_dir, "Figure2_revised_grouped_ABCD.png"), Figure2, width = 18, height = 14, dpi = 300)