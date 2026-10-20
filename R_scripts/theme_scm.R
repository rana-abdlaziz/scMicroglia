# Shared style for the scMicroglia manuscript figures (R / ggplot2).
# source("theme_scm.R") then build panels with theme_scm() and save with save_fig().
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel); library(scales); library(ragg)
})

# ---- page geometry (Alzheimer's & Dementia / Wiley) ----------------------------------
FIG_W_MM  <- 183     # double-column width
FIG_W1_MM <- 88      # single-column width
FIG_HMAX_MM <- 230   # maximum height of a full-page figure

# ---- semantic colours (identical to the earlier Python figure set) --------------------
COL <- list(
  men = "#2a78d6", women = "#eb6834", dimorphic = "#4a3aa7", shared = "#52514e",
  up = "#e34948", down = "#2a78d6", neutral = "#f0efec",
  ink = "#0b0b0b", ink2 = "#52514e", muted = "#9a9893", grid = "#e6e5e1",
  highlight = "#eda100", hpc = "#2166ac", pfc = "#e08214",
  act = "#e34948", inh = "#2a78d6", uns = "#c9c7c1", both = "#4a3aa7"
)
CAT_COLS <- c("sex-shared" = COL$shared, "female-specific" = COL$women,
              "male-specific" = COL$men, "sex-dimorphic" = COL$dimorphic)
# diverging log2FC scale (centred at 0, clipped at +/-4 as in Figure 1D)
DIV_COLS <- c("#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f2a3a2", "#e34948", "#a32a29")
scale_fill_lfc <- function(limits = c(-4, 4), name = expression(log[2]*FC), ...)
  scale_fill_gradientn(colours = DIV_COLS, limits = limits, oob = scales::squish,
                       breaks = c(-4, -2, 0, 2, 4), labels = c("\u2264\u22124", "\u22122", "0", "2", "\u22654"),
                       name = name, ...)
SEQ_COLS <- c("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b")

# ---- typography: three sizes (pt) + panel letters ----------------------------------------
FS <- list(base = 7, small = 6.2, tick = 6, title = 7.5, tag = 10)
FONT <- "Arial"
pt <- function(x) x / ggplot2::.pt          # convert pt to geom_text size units

theme_scm <- function(base = FS$base, grid = FALSE) {
  th <- theme_classic(base_size = base, base_family = FONT) +
    theme(
      text = element_text(colour = COL$ink),
      axis.text = element_text(size = FS$tick, colour = COL$ink2),
      axis.title = element_text(size = base),
      axis.line = element_line(linewidth = 0.3, colour = COL$ink2),
      axis.ticks = element_line(linewidth = 0.3, colour = COL$ink2),
      axis.ticks.length = unit(1.5, "pt"),
      plot.title = element_text(size = FS$title, face = "bold", hjust = 0, margin = margin(b = 3, l = 11)),
      plot.subtitle = element_text(size = FS$small, colour = COL$ink2, hjust = 0, margin = margin(b = 3, l = 11)),
      plot.title.position = "plot",
      legend.title = element_text(size = FS$small),
      legend.text = element_text(size = FS$small),
      legend.key.size = unit(7, "pt"),
      legend.background = element_blank(), legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base, face = "bold", hjust = 0),
      plot.margin = margin(3, 3, 3, 3),
      plot.tag = element_text(size = FS$tag, face = "bold", family = FONT)
    )
  if (grid) th <- th + theme(panel.grid.major = element_line(colour = COL$grid, linewidth = 0.25))
  th
}
theme_void_scm <- function() theme_void(base_size = FS$base, base_family = FONT) +
  theme(plot.title = element_text(size = FS$title, face = "bold", hjust = 0, margin = margin(b = 3, l = 11)),
        plot.subtitle = element_text(size = FS$small, colour = COL$ink2, hjust = 0),
        plot.title.position = "plot", plot.margin = margin(3, 3, 3, 3),
        legend.text = element_text(size = FS$small), legend.title = element_text(size = FS$small),
        plot.tag = element_text(size = FS$tag, face = "bold", family = FONT))

# patchwork panel letters (upper-case bold, A&D style): add_tags(p1 | p2)
add_tags <- function(p, ...) p + plot_annotation(tag_levels = "A", ...) &
  theme(plot.tag = element_text(size = FS$tag, face = "bold", family = FONT, hjust = 0, vjust = 1), plot.tag.location = "plot", plot.tag.position = c(0, 1))

# Unicode minus for negative numbers on axes / labels
lab_minus <- function(x, ...) sub("^-", "\u2212", scales::label_number(...)(x))
fmt_minus <- function(s) gsub("-", "\u2212", s, fixed = TRUE)
# P-value formatter: 1.5 x 10^-4 style (plotmath string, use parse = TRUE)
fmt_p_math <- function(p, digits = 2) {
  vapply(p, function(v) {
    if (is.na(v)) return("NA")
    if (v >= 0.001) return(sprintf("'%s'", formatC(signif(v, digits), format = "fg", digits = digits)))
    e <- floor(log10(v)); m <- signif(v / 10^e, digits)
    sprintf("%s %%*%% 10^%d", formatC(m, format = "f", digits = digits - 1), e)
  }, character(1))
}

# ---- export: vector PDF (cairo) + 600-dpi TIFF (LZW) + 600-dpi PNG -----------------------
save_fig <- function(plot, name, width_mm = FIG_W_MM, height_mm, outdir = "out", dpi = 600) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(height_mm <= FIG_HMAX_MM + 1)
  f <- file.path(outdir, name)
  ggsave(paste0(f, ".pdf"), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(paste0(f, ".tif"), plot, width = width_mm, height = height_mm, units = "mm", dpi = dpi,
         device = ragg::agg_tiff, compression = "lzw")
  ggsave(paste0(f, ".png"), plot, width = width_mm, height = height_mm, units = "mm", dpi = dpi,
         device = ragg::agg_png)
  invisible(paste0(f, c(".pdf", ".tif", ".png")))
}
