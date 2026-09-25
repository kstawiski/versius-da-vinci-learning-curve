# PostCDDP house style for manuscript6 figures (base_size 8, Open Sans, final print size).
suppressPackageStartupMessages({ library(ggplot2); library(ragg); library(systemfonts) })
if (!("Open Sans" %in% unique(systemfonts::system_fonts()$family))) stop("Open Sans required")
theme_postcddp <- function(base_size = 8, base_family = "Open Sans", grid_col = "#d7d7d7") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(plot.title = element_text(face = "bold", size = base_size + 2, color = "#1a1a1a"),
          axis.title = element_text(face = "bold", size = base_size + 1, color = "#222222"),
          axis.text = element_text(size = base_size - 1, color = "#222222"),
          legend.title = element_text(face = "bold", size = base_size), legend.text = element_text(size = base_size - 1),
          strip.text = element_text(face = "bold", size = base_size, hjust = 0),
          panel.grid.major = element_line(color = grid_col, linewidth = 0.4), panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.key.size = unit(0.8, "lines"), panel.spacing = unit(4, "pt"),
          plot.margin = margin(2, 2, 2, 2, "mm"), plot.title.position = "plot", text = element_text(lineheight = 0.9))
}
PAL_PLATFORM <- c("Versius" = "#B45A1E", "da Vinci" = "#1E3C78", "Laparoscopic" = "#9A9A9A")
save_fig <- function(plot, stem, width, height) {
  out <- file.path(WS, "figures")
  dir.create(out, showWarnings = FALSE)
  ggsave(file.path(out, paste0(stem, ".png")), plot, width = width, height = height, units = "mm", dpi = 320, device = ragg::agg_png, bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = width, height = height, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(out, paste0(stem, ".svg")), plot, width = width, height = height, units = "mm", device = svglite::svglite, bg = "white")
  ggsave(file.path(out, paste0(stem, ".tiff")), plot, width = width, height = height, units = "mm", dpi = 600, device = ragg::agg_tiff, bg = "white", compression = "lzw")
}
