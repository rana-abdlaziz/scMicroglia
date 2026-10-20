# Figure 2 - Sex modulates the magnitude, not the identity, of the shared program (R / ggplot2 redraw)
# Run from the working directory containing theme_scm.R and inputs/ (unzipped scm_figure_inputs.zip).
suppressPackageStartupMessages({
  library(grid); library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
  library(ggstar); library(ggtext); library(ggh4x)
})
source("theme_scm.R")
SRC <- "inputs/source_data"
LOWGREY <- "#bdbbb6"

# ---- local helpers ---------------------------------------------------------------------
mm <- function(x) unit(x, "mm")
fix_panel <- function(w, h) force_panelsizes(cols = mm(w), rows = mm(h))
num <- function(v, d = 3, plus = FALSE) fmt_minus(sprintf(paste0("%", if (plus) "+" else "", ".", d, "f"), v))
# P value as markdown: 0.83 / 0.049 / 1.51 x 10^-4 (3 significant digits below 0.01, matching the legends)
fmt_p_md <- function(p) {
  if (p >= 0.1) return(sprintf("%.2f", p))
  if (p >= 0.01) return(sprintf("%.3f", p))
  e <- floor(log10(p)); m <- p / 10^e
  sprintf("%.2f \u00d7 10<sup>%s</sup>", m, fmt_minus(as.character(e)))
}
xlab_e <- function(who, extra = "") bquote(log[2] * FC * ", " * .(who) * .(extra))
ylab_abs <- function(who) bquote("|" * log[2] * FC * "|, " * .(who))
BIAS_LAB <- expression("|" * log[2] * FC[F] * "| \u2212 |" * log[2] * FC[M] * "|")
th_tag <- theme(plot.tag = element_text(size = FS$tag, face = "bold", family = FONT, hjust = 0, vjust = 1),
                plot.tag.position = c(0, 1), plot.tag.location = "plot")

# =========================================================================================
# A  women vs men log2FC of sex-shared genes (HPC, PFC)
# =========================================================================================
a <- fread(file.path(SRC, "Figure2A.csv"))
stA <- a[, .(n = .N, rho = cor(lfc_m, lfc_f, method = "spearman"), r = cor(lfc_m, lfc_f)), by = region]
panelA <- function(reg, lim, brk, tag = NULL) {
  s <- a[region == reg]; st <- stA[region == reg]
  lab <- sprintf("%s, n = %s\nSpearman \u03c1 = %.3f\nPearson r = %.3f", reg,
                 format(st$n, big.mark = ","), st$rho, st$r)
  ggplot(s, aes(lfc_m, lfc_f)) +
    geom_hline(yintercept = 0, colour = COL$grid, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = COL$grid, linewidth = 0.3) +
    annotate("segment", x = lim[1], y = lim[1], xend = lim[2], yend = lim[2], colour = COL$ink2, linewidth = 0.3) +
    geom_point(size = 0.55, colour = COL$shared, alpha = 0.45, stroke = 0) +
    annotate("text", x = lim[1] + 0.04 * diff(lim), y = lim[2] - 0.03 * diff(lim), label = lab,
             hjust = 0, vjust = 1, size = pt(FS$tick), lineheight = 0.95, family = FONT, colour = COL$ink) +
    scale_x_continuous(breaks = brk, labels = lab_minus) +
    scale_y_continuous(breaks = brk, labels = lab_minus) +
    coord_cartesian(xlim = lim, ylim = lim, expand = FALSE) +
    labs(x = xlab_e("men"), y = xlab_e("women"), tag = tag) +
    theme_scm() + th_tag + theme(plot.margin = margin(2, 4, 2, 2)) + fix_panel(34, 34)
}
pA1 <- panelA("HPC", c(-3.5, 11.2), seq(0, 10, 5))
pA2 <- panelA("PFC", c(-4, 7.5), seq(-4, 6, 2))

# =========================================================================================
# B  |log2FC| women vs men, 109 core genes, log-log, 38 high-effect genes labelled
# =========================================================================================
b <- fread(file.path(SRC, "Figure2B.csv"))
SHP <- c("other" = 15L, "state marker" = 11L, "feedback regulator" = 28L, "state marker + feedback regulator" = 1L)
SZ  <- c("other" = 1.1, "state marker" = 1.7, "feedback regulator" = 1.7, "state marker + feedback regulator" = 2.3)
b[, `:=`(shp = SHP[marker_class], sz = SZ[marker_class],
         fillc = fifelse(stratum == "low", LOWGREY, fifelse(marker_class == "other", COL$ink2, COL$ink)))]
BLIM <- c(0.07, 14); BXLIM <- c(0.07, 40); B_W <- 61; B_H <- 76
PITCH_MM <- 2.32                                        # vertical label pitch (6-pt italic labels)
pitch_dec <- PITCH_MM / B_H * log10(BLIM[2] / BLIM[1])  # same pitch in decades of the y axis
place_column <- function(want, top, bot, pitch) {      # deterministic call-out column (as in fig2.py)
  pos <- numeric(length(want))
  for (j in seq_along(want)) pos[j] <- if (j == 1) min(want[j], top) else min(want[j], pos[j - 1] - pitch)
  for (j in rev(seq_along(want))) pos[j] <- max(pos[j], if (j == length(want)) bot else pos[j + 1] + pitch)
  stopifnot(pos[1] <= top + 1e-9, all(-diff(pos) >= pitch - 1e-9))
  pos
}
statB <- list()
panelB <- function(reg, first, tag = NULL) {
  d <- b[region == reg]; hi <- d[stratum == "high"]
  dd <- hi$abs_lfc_women - hi$abs_lfc_men
  w <- wilcox.test(hi$abs_lfc_women, hi$abs_lfc_men, paired = TRUE, exact = TRUE)
  statB[[reg]] <<- data.table(region = reg, n = nrow(hi), n_F_gt_M = sum(dd > 0), n_M_gt_F = sum(dd < 0),
                              median = median(dd), P = w$p.value)
  labs_dt <- rbindlist(lapply(list(list(side = "L", sel = dd >= 0, xa = 1.08, ylo = 0.72, yhi = 13.0),
                                   list(side = "R", sel = dd < 0, xa = 8.2, ylo = 0.085, yhi = 7.5)), function(s) {
    h <- hi[s$sel][order(-abs_lfc_women)]
    pos <- place_column(log10(h$abs_lfc_women), log10(s$yhi), log10(s$ylo), pitch_dec)
    h[, .(gene, x = abs_lfc_men, y = abs_lfc_women, side = s$side, lx = s$xa, ly = 10^pos,
          x0 = if (s$side == "L") s$xa * 1.06 else s$xa / 1.06)]
  }))
  sub <- sprintf("High-effect genes (n = %d): women &gt; men %d, men &gt; women %d<br>median \u0394 = %s, Wilcoxon *P* = %s",
                 nrow(hi), sum(dd > 0), sum(dd < 0), num(median(dd), 3, plus = TRUE), fmt_p_md(w$p.value))
  tk <- c(0.1, 0.5, 1.5, 5, 10)
  ggplot(d, aes(abs_lfc_men, abs_lfc_women)) +
    annotate("segment", x = BXLIM[1], y = BXLIM[1], xend = BLIM[2], yend = BLIM[2], colour = COL$ink2, linewidth = 0.3) +
    geom_segment(data = labs_dt, aes(x = x0, y = ly, xend = x, yend = y), inherit.aes = FALSE,
                 colour = COL$muted, linewidth = 0.2) +
    geom_star(data = d[stratum == "low"], aes(starshape = shp, size = sz, fill = fillc),
              colour = "white", starstroke = 0.25) +
    geom_star(data = d[stratum == "high"], aes(starshape = shp, size = sz, fill = fillc),
              colour = "white", starstroke = 0.25) +
    geom_text(data = labs_dt, aes(x = lx, y = ly, label = gene, hjust = ifelse(side == "L", 1, 0)),
              inherit.aes = FALSE, size = pt(FS$tick), fontface = "italic", family = FONT, colour = COL$ink) +
    scale_starshape_identity() + scale_size_identity() + scale_fill_identity() +
    scale_x_log10(breaks = tk, labels = format(tk, drop0trailing = TRUE, trim = TRUE)) +
    scale_y_log10(breaks = tk, labels = if (first) format(tk, drop0trailing = TRUE, trim = TRUE) else NULL) +
    coord_cartesian(xlim = BXLIM, ylim = BLIM, expand = FALSE, clip = "off") +
    labs(x = ylab_abs("men"), y = if (first) ylab_abs("women") else NULL, title = reg, subtitle = sub, tag = tag) +
    theme_scm() + th_tag +
    theme(plot.subtitle = element_markdown(size = FS$small, colour = COL$ink2, lineheight = 1.15,
                                           margin = margin(b = 4, l = 0)),
          plot.title = element_text(margin = margin(b = 2, l = 0)),
          plot.title.position = "panel",
          plot.margin = margin(2, 4, 2, if (first) 2 else 6)) +
    fix_panel(B_W, B_H)
}
pB1 <- panelB("HPC", TRUE)
pB2 <- panelB("PFC", FALSE)

# shared marker legend for B (drawn as a small strip)
legB <- data.table(lab = c("State marker", "Feedback regulator", "Both", "Other, high effect", "Other, low effect"),
                   shp = c(11L, 28L, 1L, 15L, 15L), sz = c(1.9, 1.9, 2.5, 1.3, 1.3),
                   fillc = c(COL$ink, COL$ink, COL$ink, COL$ink2, LOWGREY))
legB[, w := systemfonts::string_width(lab, family = FONT, size = FS$small) * 25.4 / 72]; legB[, x := cumsum(c(1.5, head(w, -1) + 7))]
LEGB_W <- ceiling(max(legB$x + 1.9 + legB$w)) + 1
pBleg <- ggplot(legB) +
  geom_star(aes(x = x, y = 0, starshape = shp, size = sz, fill = fillc), colour = "white", starstroke = 0.2) +
  geom_text(aes(x = x + 1.9, y = 0, label = lab), hjust = 0, size = pt(FS$small), family = FONT, colour = COL$ink) +
  scale_starshape_identity() + scale_size_identity() + scale_fill_identity() +
  coord_cartesian(xlim = c(0, LEGB_W), ylim = c(-2, 2), expand = FALSE, clip = "off") +
  theme_void_scm() + theme(plot.margin = margin(0, 0, 0, 0))

# =========================================================================================
# C  dumbbell: sex bias of the 38 high-effect core genes, HPC vs PFC
# =========================================================================================
hc <- fread(file.path(SRC, "Figure2C.csv"))
hc[, `:=`(bias_HPC_chk = abs(HPC_women) - abs(HPC_men), bias_PFC_chk = abs(PFC_women) - abs(PFC_men))]
stopifnot(max(abs(hc$bias_HPC - hc$bias_HPC_chk)) < 1e-9, max(abs(hc$bias_PFC - hc$bias_PFC_chk)) < 1e-9)
hc[, sw := (bias_HPC > 0) & (bias_PFC < 0)]
stopifnot(identical(hc$sw, as.logical(hc$switch_F_to_M)))
setorder(hc, -bias_HPC); hc[, x := .I]
wc <- wilcox.test(hc$bias_HPC, hc$bias_PFC, paired = TRUE, exact = TRUE)
nsw <- sum(hc$sw)
hc[, lab_md := fifelse(sw, sprintf("<b><i>%s</i></b>", gene),
                       sprintf("<span style='color:%s'><i>%s</i></span>", COL$ink2, gene))]
CY <- c(-2.6, 3.8); NX <- nrow(hc)
legC_y <- 3.5; lx0 <- 7
pC <- ggplot(hc) +
  geom_hline(yintercept = 0, colour = COL$ink2, linewidth = 0.3) +
  geom_segment(data = hc[sw == FALSE], aes(x = x, xend = x, y = bias_HPC, yend = bias_PFC),
               colour = "#c9c7c2", linewidth = 0.45, lineend = "butt") +
  geom_segment(data = hc[sw == TRUE], aes(x = x, xend = x, y = bias_HPC, yend = bias_PFC),
               colour = COL$men, linewidth = 0.85, lineend = "butt") +
  geom_point(aes(x, bias_PFC), shape = 16, size = 1.55, colour = COL$ink) +
  geom_point(aes(x, bias_HPC), shape = 21, size = 1.55, stroke = 0.45, fill = "white", colour = COL$ink) +
  # in-panel key
  annotate("point", x = lx0, y = legC_y, shape = 21, size = 1.55, stroke = 0.45, fill = "white", colour = COL$ink) +
  annotate("text", x = lx0 + 0.45, y = legC_y, label = "HPC", hjust = 0, size = pt(FS$small), family = FONT) +
  annotate("point", x = lx0 + 3.2, y = legC_y, shape = 16, size = 1.55, colour = COL$ink) +
  annotate("text", x = lx0 + 3.65, y = legC_y, label = "PFC", hjust = 0, size = pt(FS$small), family = FONT) +
  annotate("segment", x = lx0 + 6.3, xend = lx0 + 7.5, y = legC_y, yend = legC_y, colour = COL$men, linewidth = 0.85) +
  annotate("text", x = lx0 + 7.9, y = legC_y, hjust = 0, size = pt(FS$small), family = FONT,
           label = sprintf("Female-biased in HPC \u2192 male-biased in PFC (n = %d; bold)", nsw)) +
  annotate("text", x = NX + 0.4, y = CY[2] - 0.08, label = "\u2191 larger in women", hjust = 1, vjust = 1,
           size = pt(FS$small), family = FONT, colour = COL$ink2) +
  annotate("text", x = NX + 0.4, y = CY[1] + 0.08, label = "\u2193 larger in men", hjust = 1, vjust = 0,
           size = pt(FS$small), family = FONT, colour = COL$ink2) +
  annotate("richtext", x = NX + 0.4, y = 2.95, hjust = 1, vjust = 1, size = pt(FS$small), family = FONT,
           fill = NA, label.colour = NA, label.padding = unit(0, "pt"), lineheight = 1.2, colour = COL$ink,
           label = sprintf("Paired Wilcoxon *P* = %s<br>median HPC %s, PFC %s", fmt_p_md(wc$p.value),
                           num(median(hc$bias_HPC)), num(median(hc$bias_PFC)))) +
  scale_x_continuous(breaks = hc$x, labels = hc$lab_md) +
  scale_y_continuous(breaks = -2:3, labels = lab_minus) +
  coord_cartesian(xlim = c(0.3, NX + 0.7), ylim = CY, expand = FALSE) +
  labs(x = NULL, y = BIAS_LAB) +
  theme_scm() + th_tag +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = FS$tick, colour = COL$ink),
        axis.ticks.x = element_blank(), plot.margin = margin(4, 4, 2, 2)) +
  fix_panel(166, 44)

# =========================================================================================
# D  58 non-shared HPC genes: colour = HPC category, shape = PFC status
# =========================================================================================
dD <- fread(file.path(SRC, "Figure2D.csv"))
CATS3 <- c("female-specific", "male-specific", "sex-dimorphic")
STATUS <- c("PFC sex-shared (filtered list)", "PFC sex-shared, unfiltered list only", "PFC other category",
            "absent from PFC lists")
stopifnot(all(dD$pfc_status %in% STATUS), all(dD$category %in% CATS3))
nCat <- dD[, .N, by = category][match(CATS3, category), N]
nSt <- dD[, .N, by = pfc_status][match(STATUS, pfc_status), N]
catLab <- setNames(sprintf("%s (%d)", CATS3, nCat), CATS3)
stLab <- setNames(sprintf("%s (%d)", STATUS, nSt), STATUS)
ST_SHP <- setNames(c(15L, 13L, 11L, 29L), STATUS); ST_ANG <- setNames(c(0, 0, 0, 45), STATUS)
ST_SZ <- setNames(c(1.9, 1.7, 1.9, 2.3), STATUS)
dD[, `:=`(category = factor(category, CATS3), pfc_status = factor(pfc_status, STATUS),
          ang = ST_ANG[as.character(pfc_status)], sz = ST_SZ[as.character(pfc_status)])]
LBL_D <- c("PRDM11", "CXCL12", "DHRS9", "MTRNR2L8", "OTUD1", "GABBR1", "BCL10", "TICAM1", "SLC27A6", "BTBD11",
           "INO80D", "SLC26A3")
stopifnot(all(LBL_D %in% dD$gene))
DX <- c(-2.3, 2.5); DY <- c(-0.95, 2.2)
# three labels next to the dense female-specific cluster (TICAM1, BCL10; sex-dimorphic INO80D) are placed explicitly
MAN_D <- data.table(gene = c('TICAM1', 'BCL10', 'INO80D'), lx = c(-0.17, 0.72, -0.62), ly = c(0.83, 0.52, 0.34), hj = c(1, 0, 0.5))
MAN_D[, `:=`(sx = lx + c(0.03, -0.03, 0), sy = ly + c(0, 0, 0.05))]
MAN_D <- merge(MAN_D, dD[, .(gene, lfc_m, lfc_f)], by = 'gene')
pD <- ggplot(dD, aes(lfc_m, lfc_f)) +
  geom_hline(yintercept = 0, colour = COL$grid, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = COL$grid, linewidth = 0.3) +
  annotate("segment", x = DX[1], y = DX[1], xend = DX[2], yend = DX[2], colour = COL$ink2, linewidth = 0.3) +
  geom_star(aes(fill = category, starshape = pfc_status, angle = ang, size = sz), colour = "white", starstroke = 0.3) +
  geom_segment(data = MAN_D, aes(x = sx, y = sy, xend = lfc_m, yend = lfc_f),
               colour = COL$muted, linewidth = 0.2) +
  geom_text(data = MAN_D, aes(lx, ly, label = gene, hjust = hj), size = pt(FS$tick), fontface = 'italic',
            family = FONT, colour = COL$ink) +
  geom_text_repel(data = dD[gene %in% setdiff(LBL_D, MAN_D$gene)], aes(label = gene), size = pt(FS$tick), fontface = "italic",
                  family = FONT, colour = COL$ink, seed = 7, box.padding = 0.3, point.padding = 0.15,
                  min.segment.length = 0, segment.colour = COL$muted, segment.size = 0.2, max.overlaps = Inf,
                  force = 2, max.time = 2, max.iter = 20000) +
  scale_fill_manual(values = CAT_COLS[CATS3], guide = "none") +
  scale_starshape_manual(values = ST_SHP, guide = "none") + scale_size_identity() +
  scale_x_continuous(breaks = -2:2, labels = lab_minus) +
  scale_y_continuous(breaks = seq(-0.5, 2, 0.5), labels = function(x) fmt_minus(sprintf("%.1f", x))) +
  coord_cartesian(xlim = DX, ylim = DY, expand = FALSE) +
  labs(x = xlab_e("men", " (HPC)"), y = xlab_e("women", " (HPC)")) +
  theme_scm() + th_tag +
  theme(plot.margin = margin(4, 2, 2, 2)) +
  fix_panel(55, 54)

# =========================================================================================
# manual legend for D (data units = mm inside a viewport of LEGD_W x LEGD_H mm)
txt_mm <- function(s, bold = FALSE) systemfonts::string_width(s, family = FONT, size = FS$small, weight = if (bold) "bold" else "normal") * 25.4 / 72
legD <- rbind(
  data.table(row = 0, type = 'title', lab = 'HPC category (color)', shp = NA_integer_, ang = 0, fillc = NA_character_, sz = NA_real_),
  data.table(row = 1:3, type = 'key', lab = unname(catLab), shp = 15L, ang = 0, fillc = unname(CAT_COLS[CATS3]), sz = 2.2),
  data.table(row = 4.5, type = 'title', lab = 'Status in PFC (shape)', shp = NA_integer_, ang = 0, fillc = NA_character_, sz = NA_real_),
  data.table(row = 5.5:8.5, type = 'key', lab = unname(stLab), shp = unname(ST_SHP), ang = unname(ST_ANG), fillc = COL$ink2,
             sz = c(2.2, 2.0, 2.2, 2.6)))
LEGD_ROW <- 3.4; legD[, y := -1.6 - row * LEGD_ROW]
LEGD_W <- ceiling(max(4 + txt_mm(legD$lab, bold = TRUE))) + 1; LEGD_H <- ceiling(-min(legD$y) + 2)
pDleg <- ggplot(legD) +
  geom_text(data = legD[type == 'title'], aes(0, y, label = lab), hjust = 0, fontface = 'bold', size = pt(FS$small),
            family = FONT, colour = COL$ink) +
  geom_star(data = legD[type == 'key'], aes(1.4, y, starshape = shp, angle = ang, fill = fillc, size = sz),
            colour = 'white', starstroke = 0.3) +
  geom_text(data = legD[type == 'key'], aes(3.3, y, label = lab), hjust = 0, size = pt(FS$small), family = FONT,
            colour = COL$ink) +
  scale_starshape_identity() + scale_size_identity() + scale_fill_identity() +
  coord_cartesian(xlim = c(0, LEGD_W), ylim = c(-LEGD_H, 0), expand = FALSE, clip = 'off') +
  theme_void_scm() + theme(plot.margin = margin(0, 0, 0, 0))

# E  PFC non-shared genes retained (8) vs removed (24) by the effect-size filter, symlog axes
# =========================================================================================
e <- fread(file.path(SRC, "Figure2E.csv"))
nKeep <- sum(e$list == "retained after filter"); nRem <- sum(e$list == "removed by filter")
stopifnot(nKeep == 8, nRem == 24)
LS <- 1.2 / (1 - 1 / 10)                                       # matplotlib symlog: linthresh 1, linscale 1.2
symlog <- function(v) ifelse(abs(v) <= 1, v * LS, sign(v) * (LS + log10(abs(v))))
e[, `:=`(tx = symlog(lfc_m), ty = symlog(lfc_f), category = factor(category, CATS3))]
EL <- symlog(c(-4, 18)); E_W <- 55; E_H <- 54
ux <- diff(EL) / E_W * 0.3528; uy <- diff(EL) / E_H * 0.3528          # data units per point
OFF_E <- data.table(gene = c("ZBTB20", "SMAP2", "VWC2", "EDIL3", "ZBTB16", "HLA-DQB1", "CCL4"),
                    dx = c(-22, -22, 7, 7, 7, 7, -6), dy = c(12, -10, 4, 3, -7, -6, 0),
                    hj = c(1, 1, 0, 0, 0, 0, 1))
lE <- merge(OFF_E, e[, .(gene, tx, ty)], by = "gene")
lE[, `:=`(lx = tx + dx * ux, ly = ty + dy * uy, lead = abs(dx) > 10)]
tkE <- c(-3, -1, 0, 1, 3, 10)
keyE <- data.table(x = EL[1] + 0.25, y = EL[2] - c(0.62, 0.95), lab = c(sprintf("Retained (%d)", nKeep),
                                                                    sprintf("Removed by filter (%d)", nRem)))
pE <- ggplot(e, aes(tx, ty)) +
  geom_hline(yintercept = 0, colour = COL$grid, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = COL$grid, linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = COL$ink2, linewidth = 0.3) +
  geom_segment(data = lE[lead == TRUE], aes(x = lx + 0.03, y = ly, xend = tx - 0.04, yend = ty),
               colour = COL$muted, linewidth = 0.2) +
  geom_star(data = e[list == "removed by filter"], aes(colour = category), starshape = 13L, fill = "white",
            size = 1.6, starstroke = 0.55) +
  geom_star(data = e[list == "retained after filter"], aes(fill = category), starshape = 15L, colour = "white",
            size = 2.0, starstroke = 0.3) +
  geom_text(data = lE, aes(lx, ly, label = gene, hjust = hj), size = pt(FS$tick), fontface = "italic",
            family = FONT, colour = COL$ink) +
  annotate("text", x = EL[1] + 0.12, y = EL[2] - 0.28, label = "Colors as in D", hjust = 0, size = pt(FS$small),
           family = FONT, fontface = "bold") +
  annotate("star", x = keyE$x[1], y = keyE$y[1], starshape = 15L, fill = COL$ink2, colour = "white", size = 2.0,
           starstroke = 0.3) +
  annotate("star", x = keyE$x[2], y = keyE$y[2], starshape = 13L, fill = "white", colour = COL$ink2, size = 1.6,
           starstroke = 0.55) +
  annotate("text", x = keyE$x + 0.2, y = keyE$y, label = keyE$lab, hjust = 0, size = pt(FS$small), family = FONT) +
  annotate("text", x = EL[2] - 0.1, y = EL[1] + 0.1, label = "symmetric-log axes\n(linear within \u00b11)",
           hjust = 1, vjust = 0, size = pt(FS$small), family = FONT, colour = COL$muted, lineheight = 0.95) +
  scale_fill_manual(values = CAT_COLS[CATS3], guide = "none") +
  scale_colour_manual(values = CAT_COLS[CATS3], guide = "none") +
  scale_x_continuous(breaks = symlog(tkE), labels = fmt_minus(as.character(tkE))) +
  scale_y_continuous(breaks = symlog(tkE), labels = fmt_minus(as.character(tkE))) +
  coord_cartesian(xlim = EL, ylim = EL, expand = FALSE) +
  labs(x = xlab_e("men", " (PFC)"), y = xlab_e("women", " (PFC)")) +
  theme_scm() + th_tag + theme(plot.margin = margin(4, 4, 2, 8)) +
  fix_panel(E_W, E_H)

# =========================================================================================
# layout: grobs with fixed panel sizes placed at explicit mm positions (panel edges aligned)
# =========================================================================================
FIG_H <- 229
measure <- function(g) {
  pp <- g$layout[grepl('^panel', g$layout$name), ]
  wmm <- convertWidth(g$widths, 'mm', TRUE); hmm <- convertHeight(g$heights, 'mm', TRUE)
  list(g = g, w = sum(wmm), h = sum(hmm), left = sum(wmm[seq_len(min(pp$l) - 1)]),
       right = sum(wmm[-seq_len(max(pp$r))]), top = sum(hmm[seq_len(min(pp$t) - 1)]),
       bottom = sum(hmm[-seq_len(max(pp$b))]))
}
cairo_pdf(tmp_pdf <- tempfile(fileext = '.pdf'), width = 183 / 25.4, height = FIG_H / 25.4)
G <- lapply(list(A1 = pA1, A2 = pA2, B1 = pB1, B2 = pB2, C = pC, D = pD, E = pE), function(p) measure(ggplotGrob(p)))
invisible(dev.off()); unlink(tmp_pdf)
pos <- list()
put <- function(name, x, y) pos[[name]] <<- c(x = x, y = y)   # x, y = bottom-left corner of the grob (mm)
# row 1: B panels hang from the top; A1/A2 panels align with the B panel top/bottom
yB <- FIG_H - 0.5 - G$B1$h
put('A1', 0, yB + G$B1$bottom + B_H - (G$A1$h - G$A1$top))
put('A2', 0, yB + G$B1$bottom - G$A2$bottom)
xB1 <- max(G$A1$w, G$A2$w) + 1
put('B1', xB1, yB); put('B2', xB1 + G$B1$w + 1, yB)
leg_y <- yB - 5
cur <- leg_y - 1
put('C', 0, cur - G$C$h); cur <- cur - G$C$h - 1
yD <- cur - G$D$h; put('D', 0, yD); put('E', 183 - G$E$w, cur - G$E$h)
bottom_used <- cur - max(G$D$h, G$E$h)
xDleg <- G$D$w + 2
widest <- c(row1 = xB1 + G$B1$w + 1 + G$B2$w, C = G$C$w, DE = xDleg + LEGD_W + G$E$w, Bleg = xB1 + G$B1$left + LEGB_W)
print(widest); print(bottom_used); stopifnot(widest <= 183 + 0.01, bottom_used >= -0.01)
vp_mm <- function(x, y, w, h) viewport(x = mm(x), y = mm(y), width = mm(w), height = mm(h), just = c('left', 'bottom'))
kids <- lapply(names(pos), function(n) grobTree(G[[n]]$g, vp = vp_mm(pos[[n]]['x'], pos[[n]]['y'], G[[n]]$w, G[[n]]$h)))
legB_x <- xB1 + G$B1$left
kids[[length(kids) + 1]] <- grobTree(ggplotGrob(pBleg), vp = vp_mm(legB_x, leg_y, LEGB_W, 4))
yDleg <- yD + G$D$bottom + 54 / 2 - LEGD_H / 2
kids[[length(kids) + 1]] <- grobTree(ggplotGrob(pDleg), vp = vp_mm(xDleg, yDleg, LEGD_W, LEGD_H))
# panel letters (bold, upper case, top-left of each panel block)
tagg <- function(l, x, y) textGrob(l, x = mm(x), y = mm(y), just = c('left', 'top'),
                                   gp = gpar(fontsize = FS$tag, fontface = 'bold', fontfamily = FONT, col = COL$ink))
top_of <- function(n) pos[[n]]['y'] + G[[n]]$h
kids <- c(kids, list(tagg('A', 0.5, FIG_H - 0.5), tagg('B', xB1 - 1.5, FIG_H - 0.5), tagg('C', 0.5, top_of('C')),
                     tagg('D', 0.5, top_of('D')), tagg('E', pos$E['x'] - 1.5, top_of('E'))))
fig <- grobTree(rectGrob(gp = gpar(col = NA, fill = 'white')), do.call(gList, kids))
cat(sprintf('free bottom space %.1f mm; widths row1 %.1f, C %.1f, DE %.1f mm\n', bottom_used, widest[1], widest[2], widest[3]))

dir.create('out', showWarnings = FALSE)
save_fig(fig, 'Figure2', width_mm = 183, height_mm = FIG_H)
invisible(gc())


# ---- numeric checks against the legend ---------------------------------------------------
cat("A:\n"); print(stA)
cat("B:\n"); print(rbindlist(statB))
cat(sprintf("C: n = %d, switches = %d (%s); paired Wilcoxon P = %.4g; medians HPC %.4f, PFC %.4f\n",
            nrow(hc), nsw, paste(hc[sw == TRUE, gene], collapse = ", "), wc$p.value,
            median(hc$bias_HPC), median(hc$bias_PFC)))
cat("D: categories", paste(catLab, collapse = "; "), "| status", paste(stLab, collapse = "; "), "\n")
cat(sprintf("E: retained %d, removed %d; CCL4 women %.2f men %.2f\n", nKeep, nRem,
            e[gene == "CCL4", lfc_f], e[gene == "CCL4", lfc_m]))
