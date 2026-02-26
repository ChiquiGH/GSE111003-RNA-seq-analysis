# ============================================================
# helperFunctions.R
# Funciones auxiliares para la práctica GSE111003
# (Mantener aquí SOLO funciones, sin setwd() ni objetos globales)
# ============================================================

# -----------------------------
# Paletas de colores (tuyas)
# -----------------------------
color.list <- function() {
  color.list.2 <- c(
    RColorBrewer::brewer.pal(12, "Paired"), "#d45b91", "#374738",
    RColorBrewer::brewer.pal(8, "Pastel2"),
    RColorBrewer::brewer.pal(8, "Pastel2"),
    "#333333", "#5D5D5D",
    "#888888", "#B3B3B3"
  )
  color.list.2[11] <- "#e3dc5b"
  color.list.2[15] <- "#60c4b4"
  color.list.2[16] <- "#7e8046"
  return(color.list.2)
}

color.numeric <- grDevices::colorRampPalette(
  colors = c("#cfcfcf", "#958da6", "#7452b3", "#070084")
)(100)

# -----------------------------
# PCA plot (tu función, con namespaces)
# Requiere: ggplot2, factoextra
# -----------------------------
plotPCA <- function(
    pcaObject, col.points, shape.points = NULL, palette,
    legend.col, point.size = 3, title = "", pcs = c(1, 2)
) {
  variance <- round(factoextra::get_eigenvalue(pcaObject)[pcs, 2], 1)
  
  p <- ggplot2::ggplot(
    data.frame(pcaObject[["x"]]),
    ggplot2::aes(
      x = .data[[paste0("PC", pcs[1])]],
      y = .data[[paste0("PC", pcs[2])]],
      color = col.points, shape = shape.points
    )
  ) +
    ggplot2::geom_point(size = point.size) +
    ggplot2::scale_color_manual(name = legend.col, values = palette) +
    ggplot2::xlab(paste0("PC", pcs[1], " (", variance[1], "%)")) +
    ggplot2::ylab(paste0("PC", pcs[2], " (", variance[2], "%)")) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21))) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  return(p)
}

# -----------------------------
# Lectura mmseq y construcción de matriz
# Requiere: data.table, stringr, purrr
# -----------------------------
get_sample_name <- function(f) {
  b <- basename(f)
  s <- stringr::str_match(b, "_bowtie_(.*?)_RNAseq")[, 2]
  if (is.na(s)) s <- stringr::str_match(b, "_bowtie_(.*?)\\.gene\\.mmseq")[, 2]
  s
}

read_one_mmseq <- function(f) {
  dt <- data.table::fread(f)
  sample <- get_sample_name(f)
  out <- dt[, .(feature_id, unique_hits)]
  data.table::setnames(out, c("feature_id", sample))
  out
}

build_count_matrix_mmseq <- function(dir_mmseq, pattern = "\\.gene\\.mmseq\\.txt(\\.gz)?$") {
  stopifnot(dir.exists(dir_mmseq))
  files <- list.files(dir_mmseq, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No se encontraron archivos mmseq en: ", dir_mmseq)
  
  dt_list <- purrr::map(files, read_one_mmseq)
  
  counts_dt <- purrr::reduce(
    dt_list,
    function(x, y) merge(x, y, by = "feature_id", all = TRUE)
  )
  counts_dt[is.na(counts_dt)] <- 0L
  
  rawCounts <- as.matrix(counts_dt[, -1])
  rownames(rawCounts) <- counts_dt$feature_id
  storage.mode(rawCounts) <- "integer"
  
  rawCounts
}

# -----------------------------
# Metadatos de muestra para GSE111003
# -----------------------------
make_samples_metadata <- function(sample_ids) {
  md <- tibble::tibble(Sample.ID = sample_ids) |>
    dplyr::mutate(
      Donor = stringr::str_extract(Sample.ID, "^HD\\d+"),
      Time = dplyr::case_when(
        stringr::str_detect(Sample.ID, "T0$") ~ "T0",
        stringr::str_detect(Sample.ID, "4h$") ~ "4h",
        stringr::str_detect(Sample.ID, "24h$") ~ "24h",
        stringr::str_detect(Sample.ID, "d6") ~ "d6",
        TRUE ~ NA_character_
      ),
      Treatment = dplyr::case_when(
        stringr::str_detect(Sample.ID, "_BG_") ~ "BG",
        stringr::str_detect(Sample.ID, "_RPMI") | stringr::str_detect(Sample.ID, "_T0") ~ "RPMI",
        stringr::str_detect(Sample.ID, "_LPS") ~ "LPS",
        TRUE ~ "OTHER"
      ),
      Condition = dplyr::case_when(
        stringr::str_detect(Sample.ID, "T0$") ~ "T0",
        stringr::str_detect(Sample.ID, "_RPMI_4h$") ~ "RPMI_4h",
        stringr::str_detect(Sample.ID, "_BG_4h$") ~ "BG_4h",
        stringr::str_detect(Sample.ID, "_RPMI_24h$") ~ "RPMI_24h",
        stringr::str_detect(Sample.ID, "_BG_24h$") ~ "BG_24h",
        stringr::str_detect(Sample.ID, "_RPMI_d6$") ~ "RPMI_d6",
        stringr::str_detect(Sample.ID, "_RPMI_d6_4h$") ~ "RPMI_d6_4h",
        stringr::str_detect(Sample.ID, "_RPMI_d6_LPS$") ~ "RPMI_d6_LPS",
        TRUE ~ "OTHER"
      )
    ) |>
    tibble::column_to_rownames("Sample.ID")
  
  md
}

subset_phase_BG_early <- function(rawCounts, samplesMetadata) {
  keep <- rownames(samplesMetadata)[samplesMetadata$Condition %in% c("T0", "RPMI_4h", "BG_4h", "RPMI_24h", "BG_24h")]
  list(
    counts = rawCounts[, keep, drop = FALSE],
    meta = samplesMetadata[keep, , drop = FALSE]
  )
}

# -----------------------------
# Densidades para comparar distribuciones
# (estilo clase)
# -----------------------------
plotDensities2 <- function(
    matrix,
    title = "",
    xlab = "",
    ylim = 0.27,
    cols = NULL,
    cutoff = NULL
) {
  nsamples <- ncol(matrix)
  grDevices::plot(stats::density(matrix[, 1]),
                  col = cols[1], lwd = 2, las = 1, ylim = c(0, ylim),
                  main = "", xlab = "")
  graphics::grid()
  graphics::title(main = title, xlab = xlab)
  if (!is.null(cutoff)) graphics::abline(v = cutoff, lty = 3)
  
  if (nsamples >= 2) {
    for (i in 2:nsamples) {
      den <- stats::density(matrix[, i])
      graphics::lines(den$x, den$y, col = cols[i], lwd = 2)
    }
  }
}