# Supplementary Figure S6 - overlap of sex-classified microglial DEGs with sex-classified bulk-proteome DEPs
# Rebuilt from Supplementary Table S10 blocks A, C, D (no Python reference script exists for this figure).
# Run from the working directory containing theme_scm.R and inputs/ (unzipped scm_figure_inputs.zip).
suppressPackageStartupMessages({ library(grid); library(data.table); library(ggplot2); library(patchwork); library(ggh4x) })
source("theme_scm.R")
SUP <- "inputs/supp_csv"
mm <- function(x) unit(x, "mm")

summ <- fread(file.path(SUP, "S10_A_summary.csv"))
hpc <- fread(file.path(SUP, "S10_C_HPC.csv"))
pfc <- fread(file.path(SUP, "S10_D_PFC_filtered.csv"))
ov <- rbind(hpc, pfc)
stopifnot(all(ov$deg_list == "filtered"))

# =========================================================================================
# A  unique genes per region x RNA class that are DEPs of the same / any sex class
# =========================================================================================
CLS <- c("sex-shared", "male-specific", "female-specific", "sex-dimorphic")
CLS_LAB <- c("sex-shared" = "shared", "male-specific" = "male-\nspecific", "female-specific" = "female-\nspecific",
             "sex-dimorphic" = "dimorphic")
sz <- summ[deg_list == "filtered" & category %in% CLS, .(region, category, n_DEGs, n_overlap)]
cnt <- ov[, .(any = uniqueN(gene), same = uniqueN(gene[same_class == "yes"])), by = .(region, category = rna_category)]
A <- merge(CJ(region = c("HPC", "PFC"), category = CLS), sz, by = c("region", "category"), all.x = TRUE)
A <- merge(A, cnt, by = c("region", "category"), all.x = TRUE)
A[is.na(any), `:=`(any = 0L, same = 0L)]
stopifnot(all(A$same == A$n_overlap))                 # Table S10 Block A n_overlap = same-class unique genes
A[, `:=`(diff = any - same, k = match(category, CLS))]
A[, x := k + fifelse(region == "PFC", 5, 0)]
A[, top_lab := fifelse(any == 0, "0", sprintf("%d / %d", same, any))]
stopifnot(all(A[category != 'sex-shared', same] == 0))   # basis of the panel-A title
A[, xlab := sprintf("%s\n(%s)", CLS_LAB[category], prettyNum(n_DEGs, big.mark = ","))]
setorder(A, x)
tot <- A[, .(n = sum(n_DEGs)), by = region]
bars <- rbind(A[, .(x, ymin = 0, ymax = same, part = "same")], A[, .(x, ymin = same, ymax = any, part = "diff")])[ymax > ymin]
PART_COL <- c(same = COL$up, diff = COL$muted)
YA <- c(0, 66)
pA <- ggplot() +
  geom_rect(data = bars, aes(xmin = x - 0.35, xmax = x + 0.35, ymin = ymin, ymax = ymax, fill = part)) +
  geom_text(data = A, aes(x, any + 1.2, label = top_lab), vjust = 0, size = pt(FS$small), family = FONT, colour = COL$ink) +
  annotate("text", x = c(2.5, 7.5), y = YA[2] - 1, vjust = 1, size = pt(FS$base), family = FONT, colour = COL$ink,
           label = sprintf("%s (%s classified DEGs)", c("Hippocampus", "Prefrontal cortex"),
                           prettyNum(tot[match(c("HPC", "PFC"), region), n], big.mark = ","))) +
  scale_fill_manual(values = PART_COL, breaks = c("same", "diff"),
                    labels = c("DEP of the same sex class", "DEP of a different sex class"), name = NULL) +
  scale_x_continuous(breaks = A$x, labels = A$xlab) +
  scale_y_continuous(breaks = seq(0, 60, 10)) +
  coord_cartesian(xlim = c(0.4, 9.6), ylim = YA, expand = FALSE) +
  labs(x = NULL, y = "DEGs that are also\nsex-classified DEPs (n)",
       title = "Few microglial DEGs are sex-classified bulk-proteome DEPs; none of the same class among sex-specific or dimorphic DEGs") +
  theme_scm() +
  theme(axis.text.x = element_text(size = FS$tick, colour = COL$ink, lineheight = 0.95), axis.ticks.x = element_blank(),
        axis.title.y = element_text(lineheight = 1.0),
        legend.position = "inside", legend.position.inside = c(0.03, 0.80), legend.justification = c(0, 1),
        legend.key.size = unit(7, "pt"), legend.key.spacing.y = unit(2, "pt"), legend.text = element_text(size = FS$small)) +
  force_panelsizes(rows = mm(52))

# =========================================================================================
# B, C  RNA vs protein log2FC for the 21 sex-shared genes that are sex-shared DEPs
# =========================================================================================
s <- ov[rna_category == "sex-shared" & same_class == "yes"]
stopifnot(nrow(s) == 21, uniqueN(s$gene) == 21, sum(s$region == "HPC") == 5, sum(s$region == "PFC") == 16)
chk <- s[, .(m = (sign(log2FC_men) == sign(prot_log2FC_men)) == (concordant_men == "yes"),
             w = (sign(log2FC_women) == sign(prot_log2FC_women)) == (concordant_women == "yes"))]
stopifnot(all(chk$m), all(chk$w))
REG_COL <- c(HPC = COL$hpc, PFC = COL$pfc); REG_SHP <- c(HPC = 21L, PFC = 22L)
BX <- c(-0.7, 2.5); BY <- c(-0.33, 0.33); P_W <- 70; P_H <- 58
ux <- diff(BX) / P_W; uy <- diff(BY) / P_H                  # data units per mm
LBL <- data.table(gene = c("HSPB1", "SLC4A7", "RASAL2", "APP", "RASSF8", "GK"),
                  dx = c(0, 1, -1, 6, -3, -2), dy = c(3.2, -3.2, -4.5, 5, 4.2, -5.5))  # label offsets (mm)
RIBO <- c("RPL10", "RPL11", "RPL18", "RPS3")
stats_BC <- list()
panelBC <- function(sex) {
  d <- copy(s)
  d[, `:=`(rx = get(paste0("log2FC_", sex)), py = get(paste0("prot_log2FC_", sex)),
           conc = get(paste0("concordant_", sex)) == "yes")]
  d[, fillc := fifelse(conc, REG_COL[region], "white")]
  nsame <- sum(d$conc)
  stats_BC[[sex]] <<- data.table(sex = sex, same_sign = nsame, n = nrow(d))
  l <- merge(LBL, d[, .(gene, rx, py)], by = "gene")
  l[, `:=`(lx = rx + dx * ux, ly = py + dy * uy, vj = fifelse(dy > 0, 0, 1),
           hj = fifelse(dx >= 5, 0, fifelse(dx <= -5, 1, 0.5)), lead = sqrt(dx^2 + dy^2) > 4.4)]
  rb <- d[gene %in% RIBO]; stopifnot(nrow(rb) == 4)
  cx <- 1.05; cy <- if (sex == "men") -0.2 else -0.19
  ttl <- sprintf("%s: same sign for %d of %d shared genes", if (sex == "men") "Men" else "Women", nsame, nrow(d))
  ggplot(d, aes(rx, py)) +
    geom_hline(yintercept = 0, colour = COL$muted, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = COL$muted, linewidth = 0.3) +
    geom_segment(data = l[lead == TRUE], aes(x = lx - (hj == 0) * 0.6 * ux, y = ly - sign(dy) * 0.8 * uy, xend = rx, yend = py + sign(dy) * 1.0 * uy),
                 colour = COL$ink2, linewidth = 0.2) +
    geom_segment(data = rb, aes(x = cx - 0.06, y = cy + 1.2 * uy, xend = rx + 0.02, yend = py - 0.9 * uy),
                 colour = COL$ink2, linewidth = 0.2) +
    geom_point(aes(shape = region, colour = region, fill = I(fillc)), size = 1.9, stroke = 0.6) +
    geom_text(data = l, aes(lx, ly, label = gene, vjust = vj, hjust = hj), size = pt(FS$tick), fontface = "italic",
              family = FONT, colour = COL$ink) +
    annotate("text", x = cx, y = cy, label = paste0(paste(RIBO[1:2], collapse = ", "), ",\n", paste(RIBO[3:4], collapse = ", ")),
             vjust = 1, hjust = 0.5, size = pt(FS$tick), fontface = "italic", family = FONT, colour = COL$ink, lineheight = 0.95) +
    scale_shape_manual(values = REG_SHP, guide = "none") + scale_colour_manual(values = REG_COL, guide = "none") +
    scale_x_continuous(breaks = seq(-0.5, 2.5, 0.5), labels = function(v) fmt_minus(sprintf("%.1f", v))) +
    scale_y_continuous(breaks = seq(-0.3, 0.3, 0.1), labels = function(v) fmt_minus(sprintf("%.1f", v))) +
    coord_cartesian(xlim = BX, ylim = BY, expand = FALSE) +
    labs(x = bquote("RNA " * log[2] * FC * ", " * .(sex) * " (microglia)"),
         y = bquote("protein " * log[2] * FC * ", " * .(sex) * " (bulk cortex)"), title = ttl) +
    theme_scm() + force_panelsizes(rows = mm(P_H), cols = mm(P_W))
}
pB <- panelBC("men"); pC <- panelBC("women")

# shared marker legend: grid grob with mm offsets from the centre of its strip (independent of strip width)
keys <- data.table(lab = c("HPC", "PFC", "same direction (RNA vs protein)", "opposite direction"),
                   shp = c(21L, 22L, 21L, 21L), col = c(COL$hpc, COL$pfc, COL$ink2, COL$ink2),
                   fil = c(COL$hpc, COL$pfc, COL$ink2, "white"))
keys[, w := systemfonts::string_width(lab, family = FONT, size = FS$small) * 25.4 / 72]
keys[, x := cumsum(c(0, head(w, -1) + 9))]
LW <- max(keys$x + 2.4 + keys$w); keys[, x := x - LW / 2]
leg_grob <- do.call(gList, c(
  lapply(seq_len(nrow(keys)), function(i) grid::pointsGrob(x = unit(0.5, "npc") + mm(keys$x[i]), y = unit(0.5, "npc"),
    pch = keys$shp[i], gp = grid::gpar(col = keys$col[i], fill = keys$fil[i], fontsize = 1.9 * .pt + 0.6 * .stroke / 2,
                                     lwd = 0.6 * .stroke / 2))),
  lapply(seq_len(nrow(keys)), function(i) grid::textGrob(keys$lab[i], x = unit(0.5, "npc") + mm(keys$x[i] + 2.4),
    y = unit(0.5, "npc"), hjust = 0, gp = grid::gpar(fontsize = FS$small, fontfamily = FONT, col = COL$ink)))))
pLeg <- ggplot() + annotation_custom(grid::grobTree(leg_grob)) + theme_void_scm() +
  theme(plot.margin = margin(4, 0, 0, 0)) + force_panelsizes(rows = mm(4))

fig <- pA + pB + pC + pLeg + plot_layout(design = "AA\nBC\nDD") +
  plot_annotation(tag_levels = list(c("A", "B", "C", ""))) &
  theme(plot.tag = element_text(size = FS$tag, face = "bold", family = FONT, hjust = 0, vjust = 1),
        plot.tag.location = "plot", plot.tag.position = c(0, 1))
dir.create("out", showWarnings = FALSE)
fig <- fig & theme(plot.background = element_rect(fill = 'white', colour = NA))
save_fig(fig, "Supplementary_Figure_S6", width_mm = 183, height_mm = 150)
invisible(gc())

# ---- numeric checks against the SI legend ------------------------------------------------
print(A[, .(region, category, n_DEGs, same, any, n_overlap)])
print(tot)
print(rbindlist(stats_BC))
cat("rows in D_PFC:", nrow(pfc), " unique genes:", uniqueN(pfc$gene), "\n")
