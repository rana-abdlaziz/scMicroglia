# Figure 3 - pathway enrichment of sex-shared, sex-specific and sex-dimorphic DEGs (R / ggplot2 redraw).
# Reproduces python_reference_scripts/fig3.py. Panels are ggplots with fixed panel sizes that are placed
# at absolute page positions (mm from the top-left corner), mirroring the inch-based layout of the
# reference so that label columns, bars and text columns line up across panels.
# Run from a fresh R session in the working directory that contains inputs/ and theme_scm.R.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(grid); library(scales)
})
if (!dir.exists("inputs/source_data")) stop("inputs/source_data not found - unzip the figure inputs first")
source("theme_scm.R")
tmp_dev <- tempfile(fileext = ".png"); ragg::agg_png(tmp_dev)   # font-metric device for grob construction (avoids the Arial-less default pdf device)
SD <- "inputs/source_data/"

# ------------------------------------------------------------------------------ helpers
ABBR <- list(c("RNA polymerase II-specific", "RNA Pol II-specific"), c("DNA-binding transcription factor", "DNA-binding TF"),
             c("transcription factor", "TF"), c("positive regulation of", "pos. reg. of"),
             c("negative regulation of", "neg. reg. of"), c("regulation of", "reg. of"),
             c("involved in", "in"), c("establishment of", "estab. of"), c("nucleoside-triphosphatase", "NTPase"),
             c(" signaling pathway", " signaling"),
             c("AGE-RAGE signaling in diabetic complications", "AGE-RAGE signaling (diabetic compl.)"),
             c("structural constituent", "struct. constituent"), c("transmembrane transporter", "transm. transporter"),
             c("- Homo sapiens (human)", ""))
short <- function(s, n = 46) vapply(s, function(x) {
  for (a in ABBR) x <- gsub(a[1], a[2], x, fixed = TRUE)
  x <- trimws(x)
  if (nchar(x) <= n) x else paste0(sub("\\s+$", "", substr(x, 1, n - 1)), "\u2026")
}, character(1), USE.NAMES = FALSE)
ratio <- function(gr) vapply(strsplit(gr, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
nlog <- function(p) -log10(p)
FDR05 <- -log10(0.05)

# page geometry -----------------------------------------------------------------
W <- 183; H <- 229            # mm
RH <- 2.75                    # row pitch (mm) shared by every bar / dot panel
LAB <- FS$small               # term-label / annotation size (6.2 pt)
xlab_lfdr <- expression("\u2212" * log[10] ~ FDR)

# place a ggplot so that its PANEL has size w x h mm with top-left corner at (x, y) mm (y from page top)
place <- function(p, x, y, w, h) {
  g <- ggplotGrob(p + theme(panel.widths = unit(w, "mm"), panel.heights = unit(h, "mm")))
  pl <- g$layout[grepl("^panel", g$layout$name), ]
  l <- min(pl$l); t <- min(pl$t)
  lw <- if (l > 1) sum(g$widths[seq_len(l - 1)]) else unit(0, "mm")
  th <- if (t > 1) sum(g$heights[seq_len(t - 1)]) else unit(0, "mm")
  vp <- viewport(x = unit(x, "mm") - lw, y = unit(H - y, "mm") + th,
                 width = sum(g$widths), height = sum(g$heights), just = c("left", "top"))
  grobTree(g, vp = vp)
}
txt <- function(label, x, y, size = LAB, face = "plain", col = COL$ink, hjust = 0, vjust = 1, parse = FALSE) {
  lab <- if (parse) parse(text = label) else label
  textGrob(lab, x = unit(x, "mm"), y = unit(H - y, "mm"), hjust = hjust, vjust = vjust,
           gp = gpar(fontfamily = FONT, fontsize = size, fontface = face, col = col))
}
letter <- function(l, x, y) txt(l, x, y, size = FS$tag, face = "bold")
ptitle <- function(s, x, y) txt(s, x, y, size = FS$title, face = "bold")

# common theme for bar / dot panels: term labels at 6.2 pt in ink, no y ticks
th_panel <- function() theme_scm() + theme(
  axis.text.y = element_text(size = LAB, colour = COL$ink, margin = margin(r = 2)),
  axis.ticks.y = element_blank(),
  axis.text.x = element_text(size = FS$tick, margin = margin(t = 1)),
  axis.title.x = element_text(size = FS$small, margin = margin(t = 1.5)),
  axis.title.y = element_blank(), plot.margin = margin(0, 0, 0, 0), legend.position = "none")
# convert an offset in mm (from the panel left edge) to data units, given panel width w and x limits
mm2x <- function(mm, w, xlim) xlim[1] + mm / w * diff(xlim)

# ------------------------------------------------------------------------------ data
A  <- fread(paste0(SD, "Fig3A_GO_shared_HPC_PFC_dotplot.csv"))
B  <- fread(paste0(SD, "Fig3B_KEGG_core109.csv"))
C  <- fread(paste0(SD, "Fig3C_GOBP_core109_top20_recomputed.csv"))
D  <- fread(paste0(SD, "Fig3D_KEGG_HPC_male_female_specific.csv"))
Ed <- fread(paste0(SD, "Fig3E_GO_dimorphic_recomputed.csv"))
A[, ratio := ratio(GeneRatio)]
stopifnot(all(A$Count == A$n_up + A$n_down), !is.unsorted(B$p.adjust), !is.unsorted(C$p.adjust))

# ------------------------------------------------------------------------------ checks against the legend
chk <- list(
  A_rows_HPC_BP = nrow(A[region == "HPC" & ont == "BP" & group == "top"]) == 8 &&
                  nrow(A[region == "HPC" & ont == "BP" & group == "down-dominated"]) == 3,
  A_rows_PFC_BP = nrow(A[region == "PFC" & ont == "BP"]) == 11,
  A_rows_MF     = nrow(A[region == "HPC" & ont == "MF"]) == 8 && nrow(A[region == "PFC" & ont == "MF"]) == 8,
  A_FDR_lt_005  = all(A$p.adjust < 0.05),
  B_10_pathways = nrow(B) == 10, B_77_annot = all(sub(".*/", "", B$GeneRatio) == "77"),
  B_25_unique   = length(unique(unlist(strsplit(B$genes, ",")))) == 25, B_FDR = all(B$p.adjust < 0.05),
  C_20_terms    = nrow(C) == 20, C_105 = all(sub(".*/", "", C$GeneRatio) == "105"),
  C_FDR_range   = sprintf("%.2e", min(C$p.adjust)) == "8.84e-08" && sprintf("%.2e", max(C$p.adjust)) == "2.38e-05",
  D_male_5      = all(sub(".*/", "", D[sex == "HPC male-specific"]$GeneRatio) == "5"),
  D_female_19   = all(sub(".*/", "", D[sex == "HPC female-specific"]$GeneRatio) == "19"),
  D_FDR         = all(D$p.adjust < 0.05),
  E_FDR         = all(Ed$p.adjust < 0.05),
  E_genes_HPC   = setequal(unique(unlist(strsplit(Ed[region == "HPC"]$geneID, "/"))), c("SLC27A6", "SLC26A3")),
  E_genes_PFC   = setequal(unique(unlist(strsplit(Ed[region == "PFC"]$geneID, "/"))), "ZBTB16"),
  E_le2_genes   = all(Ed$Count <= 2))
print(unlist(chk))

# ------------------------------------------------------------------------------ panel A builders
NORM <- c(1.3, 12)                                     # colour limits (-log10 FDR), clipped above 12
seq_pal <- gradient_n_pal(SEQ_COLS)
fdr_col <- function(v) seq_pal(rescale(pmin(pmax(v, NORM[1]), NORM[2]), from = NORM))
dot_size <- function(cnt) sqrt(2 + 1.1 * cnt) / (0.75 * .pt) * 1.1   # diameter ~ sqrt(marker area of the reference)

ydots <- function(d) {
  n <- nrow(d); y <- rev(seq_len(n) - 1)
  ndn <- sum(d$group == "down-dominated")
  if (ndn) y[(n - ndn + 1):n] <- y[(n - ndn + 1):n] - 1
  list(y = y, ndn = ndn, ylim = c(min(y) - 0.55, max(y) + 0.55))
}
A_dot <- function(d, ont, w) {
  yy <- ydots(d); d$y <- yy$y
  d$fill <- fdr_col(nlog(d$p.adjust)); d$sz <- dot_size(d$Count)
  xlim <- c(0, max(d$ratio) * 1.2)
  br <- pretty(xlim, n = 3); br <- br[br <= xlim[2]]
  p <- ggplot(d, aes(ratio, y)) +
    geom_vline(xintercept = br, colour = COL$grid, linewidth = 0.25) +
    geom_point(aes(fill = fill, size = sz), shape = 21, colour = COL$ink2, stroke = 0.25) +
    scale_fill_identity() + scale_size_identity() +
    scale_y_continuous(breaks = d$y, labels = short(d$Description, 40), expand = c(0, 0)) +
    scale_x_continuous(breaks = br, labels = label_number(), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim, ylim = yy$ylim, clip = "off") +
    labs(x = "GeneRatio") + th_panel() +
    annotate("text", x = xlim[2], y = yy$ylim[2] + 0.25, label = paste0("GO-", ont), hjust = 1, vjust = 0,
             size = pt(LAB), fontface = "bold", colour = COL$ink2, family = FONT)
  if (yy$ndn) {
    ysep <- max(yy$y[(nrow(d) - yy$ndn + 1):nrow(d)]) + 1
    p <- p + annotate("segment", x = xlim[1], xend = xlim[2], y = ysep - 0.45, yend = ysep - 0.45,
                      colour = COL$muted, linewidth = 0.25) +
      annotate("text", x = mm2x(-1.2, w, xlim), y = ysep + 0.08, label = "best terms dominated by decreased genes",
               hjust = 1, size = pt(FS$tick), fontface = "italic", colour = COL$ink2, family = FONT)
  }
  list(p = p, h = diff(yy$ylim) * RH, ylim = yy$ylim, y = yy$y)
}
A_bar <- function(d, yinfo, header) {
  d$y <- yinfo$y; d$fu <- d$n_up / d$Count
  p <- ggplot(d) +
    geom_rect(aes(xmin = 0, xmax = fu, ymin = y - 0.32, ymax = y + 0.32), fill = COL$up) +
    geom_rect(aes(xmin = fu, xmax = 1, ymin = y - 0.32, ymax = y + 0.32), fill = COL$down) +
    geom_text(aes(x = 1.06, y = y, label = paste0(n_up, "/", n_down)), hjust = 0, size = pt(FS$tick),
              colour = COL$ink2, family = FONT) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = yinfo$ylim, clip = "off") + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  if (header) p <- p + annotate("text", x = 0, y = yinfo$ylim[2] + 0.25, label = "up/down", hjust = 0, vjust = 0,
                                size = pt(FS$tick), colour = COL$ink2, family = FONT)
  p
}

grobs <- list(rectGrob(gp = gpar(fill = "white", col = NA)))
add <- function(g) grobs[[length(grobs) + 1]] <<- g

# ---------------- A ----------------------------------------------------------------
DW <- 23; UW <- 7.6
xcol <- c(HPC = 44, PFC = 135)
yA_BP <- 9.5
add(letter("A", 0, 0.4))
nset <- c(HPC = "253", PFC = "1,128")
yA_bottom <- 0
for (reg in c("HPC", "PFC")) {
  add(ptitle(sprintf("%s sex-shared genes (n = %s)", reg, nset[reg]), xcol[reg] - 38.4, 0.9))
  yy <- yA_BP
  for (o in c("BP", "MF")) {
    d <- as.data.frame(A[region == reg & ont == o])
    a <- A_dot(d, o, DW)
    add(place(a$p, xcol[reg], yy, DW, a$h))
    add(place(A_bar(d, a, header = o == "BP"), xcol[reg] + DW + 1.5, yy, UW, a$h))
    if (o == "BP") { if (reg == "HPC") yMF <- yy + a$h + 10.0 }
    else yA_bottom <- max(yA_bottom, yy + a$h)
    yy <- yMF
  }
}

# legend strip for A --------------------------------------------------------------
yL <- yA_bottom + 9.0             # vertical centre of the strip (mm from top)
cb_x0 <- 27; cb_x1 <- 53; cb_h <- 1.8
v <- seq(NORM[1], NORM[2], length.out = 120)
cb <- data.frame(x0 = cb_x0 + (seq_along(v) - 1) / length(v) * (cb_x1 - cb_x0))
cb$x1 <- cb$x0 + (cb_x1 - cb_x0) / length(v) + 0.02; cb$fill <- fdr_col(v)
tri <- data.frame(x = c(cb_x1, cb_x1 + 1.8, cb_x1), y = c(-cb_h / 2, 0, cb_h / 2))
outline <- data.frame(x = c(cb_x0, cb_x1, cb_x1 + 1.8, cb_x1, cb_x0), y = c(-cb_h / 2, -cb_h / 2, 0, cb_h / 2, cb_h / 2))
tk <- c(5, 10); tkx <- cb_x0 + rescale(tk, from = NORM) * (cb_x1 - cb_x0)
gk <- data.frame(cnt = c(5, 25, 50, 75), x = 76 + (0:3) * 10)
leg <- ggplot() +
  geom_rect(data = cb, aes(xmin = x0, xmax = x1, ymin = -cb_h / 2, ymax = cb_h / 2, fill = I(fill))) +
  geom_polygon(data = tri, aes(x, y), fill = SEQ_COLS[length(SEQ_COLS)]) +
  geom_path(data = outline, aes(x, y), linewidth = 0.25, colour = COL$ink2) +
  annotate("segment", x = tkx, xend = tkx, y = -cb_h / 2, yend = -cb_h / 2 - 0.6, linewidth = 0.25, colour = COL$ink2) +
  annotate("text", x = tkx, y = -cb_h / 2 - 0.9, label = tk, vjust = 1, size = pt(FS$tick), colour = COL$ink2, family = FONT) +
  annotate("text", x = cb_x0 - 1.2, y = 0, label = "'\u2212'*log[10]~FDR", parse = TRUE, hjust = 1,
           size = pt(LAB), colour = COL$ink, family = FONT) +
  annotate("text", x = 66, y = 0, label = "genes", hjust = 0, size = pt(LAB), colour = COL$ink, family = FONT) +
  geom_point(data = gk, aes(x, 0, size = I(dot_size(cnt))), shape = 21, fill = "#d8d6d0", colour = COL$ink2, stroke = 0.25) +
  geom_text(data = gk, aes(x + 2.4, 0, label = cnt), hjust = 0, size = pt(FS$tick), colour = COL$ink2, family = FONT) +
  annotate("rect", xmin = 117, xmax = 119.6, ymin = -0.9, ymax = 0.9, fill = COL$up) +
  annotate("text", x = 120.6, y = 0, label = "increased", hjust = 0, size = pt(LAB), colour = COL$ink, family = FONT) +
  annotate("rect", xmin = 135, xmax = 137.6, ymin = -0.9, ymax = 0.9, fill = COL$down) +
  annotate("text", x = 138.6, y = 0, label = "decreased member genes (n up/down)", hjust = 0, size = pt(LAB),
           colour = COL$ink, family = FONT) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, W), ylim = c(-2, 2), clip = "off") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
add(place(leg, 0, yL - 2, W, 4))

# ---------------- B ----------------------------------------------------------------
yB <- yL + 5.5
add(letter("B", 0, yB)); add(ptitle("KEGG pathways of the 109 core genes (77 KEGG-annotated)", 5.6, yB + 0.5))
d <- as.data.frame(B); n <- nrow(d); d$y <- rev(seq_len(n) - 1)
wB <- 32; xlimB <- c(0, 3); ylimB <- c(-0.55, n - 1 + 0.55)
xn <- mm2x(wB + 1.4, wB, xlimB); xg <- mm2x(wB + 7.6, wB, xlimB)
pB <- ggplot(d) +
  geom_col(aes(nlog(p.adjust), y), orientation = "y", width = 0.66, fill = COL$shared) +
  geom_vline(xintercept = FDR05, colour = COL$muted, linewidth = 0.3, linetype = "22") +
  geom_text(aes(xn, y, label = GeneRatio), hjust = 0, size = pt(LAB), colour = COL$ink, family = FONT) +
  geom_text(aes(xg, y, label = gsub(",", ", ", genes)), hjust = 0, size = pt(LAB), colour = COL$ink2,
            fontface = "italic", family = FONT) +
  annotate("text", x = c(xn, xg), y = ylimB[2] + 0.35, label = c("n/N", "member genes"), hjust = 0, vjust = 0,
           size = pt(LAB), fontface = "bold", colour = COL$ink2, family = FONT) +
  scale_y_continuous(breaks = d$y, labels = short(d$Description, 40), expand = c(0, 0)) +
  scale_x_continuous(breaks = 0:3, expand = c(0, 0)) +
  coord_cartesian(xlim = xlimB, ylim = ylimB, clip = "off") +
  labs(x = expression("\u2212" * log[10] ~ FDR ~ "(dashed line: FDR = 0.05)")) + th_panel()
yBp <- yB + 8.5; hB <- diff(ylimB) * RH
add(place(pB, 44, yBp, wB, hB))

# ---------------- C ----------------------------------------------------------------
yC <- yBp + hB + 8.5
add(letter("C", 0, yC)); add(ptitle("GO-BP of the 109 core genes (top 20, recomputed)", 5.6, yC + 0.5))
d <- as.data.frame(C); n <- nrow(d); d$y <- rev(seq_len(n) - 1)
wC <- 33; xlimC <- c(0, 9); ylimC <- c(-0.55, n - 1 + 0.55)
pC <- ggplot(d) +
  geom_col(aes(nlog(p.adjust), y), orientation = "y", width = 0.66, fill = COL$shared) +
  geom_text(aes(nlog(p.adjust) + 0.15, y, label = GeneRatio), hjust = 0, size = pt(LAB), colour = COL$ink2, family = FONT) +
  scale_y_continuous(breaks = d$y, labels = short(d$Description, 40), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 8, 2), expand = c(0, 0)) +
  coord_cartesian(xlim = xlimC, ylim = ylimC, clip = "off") + labs(x = xlab_lfdr) + th_panel()
yCp <- yC + 7.5; hC <- diff(ylimC) * RH
add(place(pC, 44, yCp, wC, hC))
add(txt("org.Hs.eg.db/GO.db 3.18; 105 of 107 mapped genes GO-BP-annotated;\nno GO-MF term at FDR < 0.05",
        5.6, yCp + hC + 8.2, size = LAB, col = COL$ink2))

# ---------------- D ----------------------------------------------------------------
xL <- 125.7; wDE <- 15.7; xHead <- xL - 30.6
add(letter("D", 90.2, yC)); add(ptitle("KEGG, HPC sex-specific genes", 95.8, yC + 0.5))
d <- as.data.frame(D); d$y <- c(3.6, 1.0, 0.0)
d$col <- ifelse(startsWith(d$sex, "HPC male"), COL$men, COL$women)
xlimD <- c(0, 2); ylimD <- c(-0.55, 4.2)
heads <- data.frame(y = c(4.45, 1.85), lab = c("male-specific (9 genes)", "female-specific (45 genes)"),
                    col = c(COL$men, COL$women))
xr <- mm2x(wDE + 1.6, wDE, xlimD); xgD <- mm2x(wDE + 7.4, wDE, xlimD)
pD <- ggplot(d) +
  geom_col(aes(nlog(p.adjust), y, fill = I(col)), orientation = "y", width = 0.66) +
  geom_vline(xintercept = FDR05, colour = COL$muted, linewidth = 0.3, linetype = "22") +
  geom_text(aes(xr, y, label = GeneRatio), hjust = 0, size = pt(LAB), colour = COL$ink, family = FONT) +
  geom_text(aes(xgD, y, label = genes), hjust = 0, size = pt(LAB), colour = COL$ink2, fontface = "italic", family = FONT) +
  geom_text(data = heads, aes(mm2x(xHead - xL, wDE, xlimD), y, label = lab), hjust = 0, size = pt(LAB),
            fontface = "bold", colour = COL$ink, family = FONT) +
  geom_tile(data = heads, aes(mm2x(xHead - xL - 1.6, wDE, xlimD), y, fill = I(col)),
            width = diff(xlimD) * 1.5 / wDE, height = 1.5 / RH) +
  scale_y_continuous(breaks = d$y, labels = short(d$Description, 34), expand = c(0, 0)) +
  scale_x_continuous(breaks = 0:2, expand = c(0, 0)) +
  coord_cartesian(xlim = xlimD, ylim = ylimD, clip = "off") + labs(x = xlab_lfdr) + th_panel()
yDp <- yCp + 1.0; hD <- diff(ylimD) * RH
add(place(pD, xL, yDp, wDE, hD))

# ---------------- E ----------------------------------------------------------------
yE <- yDp + hD + 9.0
add(letter("E", 90.2, yE)); add(ptitle("GO, sex-dimorphic genes (recomputed)", 95.8, yE + 0.5))
add(txt(paste0('"HPC: "*italic("SLC27A6")*", "*italic("SLC26A3")*", "*italic("INO80D")*" ("*italic("BTBD11")*',
               '" unmapped); PFC: "*italic("ZBTB16")*"."'), 95.8, yE + 4.6, parse = TRUE))
add(txt("1\u20132 genes per term: restates single-gene annotations", 95.8, yE + 7.4))
d <- as.data.frame(Ed)
ys <- numeric(0); yv <- 0; hd <- list(); prev <- NULL
for (i in rev(seq_len(nrow(d)))) {
  key <- paste(d$region[i], d$ont[i])
  if (!is.null(prev) && key != prev) { hd[[length(hd) + 1]] <- c(yv - 0.1, prev); yv <- yv + 1.25 }
  ys <- c(ys, yv); yv <- yv + 1; prev <- key
}
hd[[length(hd) + 1]] <- c(yv - 0.1, prev)
d$y <- rev(ys)
heads <- data.frame(y = as.numeric(sapply(hd, `[`, 1)) + 0.05,
                    lab = sub("^(HPC|PFC) (BP|MF)$", "\\1 \u00b7 GO-\\2", sapply(hd, `[`, 2)))
xlimE <- c(0, 3.2); ylimE <- c(-0.55, yv + 0.3)
xr <- mm2x(wDE + 1.6, wDE, xlimE); xgE <- mm2x(wDE + 7.4, wDE, xlimE)
pE <- ggplot(d) +
  geom_col(aes(nlog(p.adjust), y), orientation = "y", width = 0.66, fill = COL$dimorphic) +
  geom_vline(xintercept = FDR05, colour = COL$muted, linewidth = 0.3, linetype = "22") +
  geom_text(aes(xr, y, label = GeneRatio), hjust = 0, size = pt(LAB), colour = COL$ink, family = FONT) +
  geom_text(aes(xgE, y, label = gsub("/", ", ", geneID)), hjust = 0, size = pt(LAB), colour = COL$ink2,
            fontface = "italic", family = FONT) +
  geom_text(data = heads, aes(mm2x(xHead - xL, wDE, xlimE), y, label = lab), hjust = 0, size = pt(LAB),
            fontface = "bold", colour = COL$ink, family = FONT) +
  scale_y_continuous(breaks = d$y, labels = short(d$Description, 34), expand = c(0, 0)) +
  scale_x_continuous(breaks = 0:3, expand = c(0, 0)) +
  coord_cartesian(xlim = xlimE, ylim = ylimE, clip = "off") + labs(x = xlab_lfdr) + th_panel()
yEp <- yE + 11.0; hE <- diff(ylimE) * RH
add(place(pE, xL, yEp, wDE, hE))
cat(sprintf("bottom of E panel: %.1f mm; C note top: %.1f mm; page H = %d mm\n", yEp + hE, yCp + hC + 8.2, H))

fig <- gTree(children = do.call(gList, grobs))
save_fig(fig, "Figure3", width_mm = W, height_mm = H, outdir = "out")
invisible(gc())
invisible(dev.off()); unlink(tmp_dev)
