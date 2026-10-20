# Figure 1 - study design, DEG categories, cross-region overlap and the 109-gene core.
# Self-contained: run from the working directory containing inputs/ and theme_scm.R
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(ggforce)
})
source("theme_scm.R")

TAGTH <- theme(plot.tag = element_text(size = FS$tag, face = "bold", family = FONT, hjust = 0, vjust = 1),
               plot.tag.position = c(0, 1), plot.tag.location = "plot")
CHECKS <- list()
chk <- function(name, got, expected) {
  ok <- isTRUE(all.equal(as.numeric(got), as.numeric(expected)))
  CHECKS[[name]] <<- sprintf("%s: got %s, legend %s -> %s", name, paste(got, collapse = "/"),
                             paste(expected, collapse = "/"), if (ok) "PASS" else "MISMATCH")
}
CAT_ORDER <- c("sex-shared", "female-specific", "male-specific", "sex-dimorphic")

# =============================================================================
# A  workflow schematic (coordinates in mm)
# =============================================================================
srcA <- read.csv("inputs/source_data/Figure1A.csv", encoding = "UTF-8")
AW <- 183; AH <- 53
BW <- 42.6; BH <- 18.6; xs <- c(0.6, 46.8, 93.0, 139.2); yrow <- c(28.9, 2.2)
FILL <- "#f7f6f3"
pos <- data.frame(step = 1:8, x = rep(xs, 2), y = rep(yrow, each = 4))
boxes <- pos %>% rowwise() %>%
  reframe(step = step, px = c(x, x + BW, x + BW, x), py = c(y, y, y + BH, y + BH))
# body lines: plain text unless a math expression is needed (log2 subscript, P_adj)
body <- list(
  c("Sun et al. 2023", "snRNA-seq, microglial nuclei", "AD and non-AD donors"),
  c("HPC: 86 samples (66 identified donors)", "PFC: 469 samples (423 donors)", "Women and men separately"),
  c("12 microglial states", "Module scores per nucleus", "Hungarian cluster\u2192state assignment"),
  c("Inflammatory MG2+MG8+MG10 (AD)", "vs homeostatic MG0 (non-AD)", "Poisson GLM; age as latent variable"),
  NULL,
  c("MATH:'Shared: |'*log[2]*'FC| >1.5 or <0.5'", "  in both sexes; specific: \u22652\u00d7",
    "  other sex; PFC 1,967 \u2192 1,136"),
  c("HPC sex-shared \u2229 PFC sex-shared", "109 genes", "(101 increased, 8 decreased)"),
  c("ORA (GO, KEGG), gseGO", "MetaCore TF network", "decoupleR TF activity"))
# consistency with the panel-A source table
for (s in c(1:4, 7, 8)) stopifnot(paste(body[[s]], collapse = "; ") == srcA$details[s])
stopifnot(grepl("1,967", srcA$details[6]), grepl("1,136", srcA$details[6]))
txt <- do.call(rbind, lapply(seq_along(body), function(s) {
  if (is.null(body[[s]])) return(NULL)
  data.frame(step = s, j = seq_along(body[[s]]), lab = body[[s]])
})) %>% left_join(pos, by = "step") %>%
  mutate(tx = x + 2.0, ty = y + BH - 7.0 - (j - 1) * 3.25,
         math = startsWith(lab, "MATH:"), lab = sub("^MATH:", "", lab))
ttl <- pos %>% mutate(title = srcA$title, cx = x + 2.9, cy = y + BH - 2.9)
cats <- data.frame(lab = CAT_ORDER, k = 0:3) %>%
  mutate(x = xs[1] + 2.4 + (k %% 2) * 20.5, y = yrow[2] + BH - 7.4 - (k %/% 2) * 4.0)
arr <- arrow(type = "closed", length = unit(1.3, "mm"), angle = 25)
segs <- expand.grid(i = 1:3, r = 1:2) %>%
  mutate(x = xs[i] + BW + 0.3, xend = xs[i + 1] - 0.3, y = yrow[r] + BH / 2)
mid <- (yrow[2] + BH + yrow[1]) / 2
pA <- ggplot() +
  geom_shape(data = boxes, aes(px, py, group = step), fill = FILL, colour = COL$ink2,
             linewidth = 0.3, radius = unit(1.4, "mm")) +
  geom_point(data = ttl, aes(cx, cy), shape = 16, size = 3.6, colour = COL$ink) +
  geom_text(data = ttl, aes(cx, cy, label = step), colour = "white", fontface = "bold",
            size = pt(6), family = FONT) +
  geom_text(data = ttl, aes(cx + 2.8, cy, label = title), hjust = 0, fontface = "bold",
            size = pt(7), family = FONT, colour = COL$ink) +
  geom_text(data = filter(txt, !math), aes(tx, ty, label = lab), hjust = 0, size = pt(6),
            family = FONT, colour = COL$ink2) +
  geom_text(data = filter(txt, math), aes(tx, ty, label = lab), hjust = 0, size = pt(6),
            family = FONT, colour = COL$ink2, parse = TRUE) +
  geom_tile(data = cats, aes(x + 1, y, fill = lab), width = 2.1, height = 2.1) +
  geom_text(data = cats, aes(x + 3.2, y, label = lab), hjust = 0, size = pt(6), family = FONT,
            colour = COL$ink2) +
  annotate("text", x = xs[1] + 2.0, y = yrow[2] + 2.3, hjust = 0, size = pt(6), family = FONT,
           colour = COL$muted, parse = TRUE,
           label = "italic(P)[adj]*', |'*log[2]*'FC| and sign per sex'") +
  geom_segment(data = segs, aes(x, y, xend = xend, yend = y), colour = COL$ink2,
               linewidth = 0.4, arrow = arr) +
  annotate("path", x = c(xs[4] + BW / 2, xs[4] + BW / 2, xs[1] + BW / 2),
           y = c(yrow[1] - 0.2, mid, mid), colour = COL$ink2, linewidth = 0.4) +
  annotate("segment", x = xs[1] + BW / 2, xend = xs[1] + BW / 2, y = mid, yend = yrow[2] + BH + 0.3,
           colour = COL$ink2, linewidth = 0.4, arrow = arr) +
  scale_fill_manual(values = CAT_COLS, guide = "none") +
  scale_x_continuous(limits = c(0, AW), expand = c(0, 0)) + scale_y_continuous(limits = c(0, AH), expand = c(0, 0)) + coord_cartesian(clip = "off") +
  labs(tag = "A") + theme_void_scm() + TAGTH + theme(plot.margin = margin(0, 0, 0, 0))

# =============================================================================
# B  numbers of classified DEGs per category (log scale)
# =============================================================================
cnt <- read.csv("inputs/source_data/Figure1B.csv")
cnt$panel <- paste0(cnt$region, ", ", cnt$list)
tot <- cnt %>% group_by(panel) %>% summarise(n = sum(total), .groups = "drop")
chk("B HPC filtered n", tot$n[tot$panel == "HPC, filtered"], 311)
chk("B PFC before filter n", tot$n[tot$panel == "PFC, before filter"], 1967)
chk("B PFC after filter n", tot$n[tot$panel == "PFC, after filter"], 1136)
stopifnot(all(cnt$increased + cnt$decreased == cnt$total))
plev <- c("HPC, filtered", "PFC, before filter", "PFC, after filter")
strip_lab <- setNames(sprintf("%s\nn = %s", plev, format(tot$n[match(plev, tot$panel)], big.mark = ",", trim = TRUE)), plev)
OFF <- 0.2
bl <- cnt %>% pivot_longer(c(increased, decreased), names_to = "dir", values_to = "n") %>%
  mutate(panel = factor(panel, plev), yy = match(category, CAT_ORDER) + ifelse(dir == "increased", -OFF, OFF),
         dir = factor(dir, c("increased", "decreased"), c("Increased", "Decreased")),
         lab = format(n, big.mark = ",", trim = TRUE))
pB <- ggplot(bl, aes(y = yy)) +
  geom_hline(yintercept = 1:4, colour = COL$grid, linewidth = 0.3) +
  geom_point(data = filter(bl, n > 0), aes(x = n, shape = dir, fill = dir), size = 1.9, colour = "white",
             stroke = 0.25) +
  geom_text(data = filter(bl, n >= 100), aes(x = n / 1.7, label = lab), hjust = 1, size = pt(6), family = FONT,
            colour = COL$ink) +
  geom_text(data = filter(bl, n > 0, n < 100), aes(x = n * 1.7, label = lab), hjust = 0, size = pt(6),
            family = FONT, colour = COL$ink) +
  geom_text(data = filter(bl, n == 0), aes(x = 0.62, label = "0"), size = pt(6), family = FONT,
            colour = COL$muted) +
  facet_wrap(~panel, nrow = 1, labeller = as_labeller(strip_lab)) +
  scale_x_log10(limits = c(0.45, 6000), breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1,000"),
                expand = c(0, 0)) +
  scale_y_reverse(breaks = 1:4, labels = CAT_ORDER, limits = c(4.5, 0.45), expand = c(0, 0)) +
  scale_shape_manual(values = c(Increased = 24, Decreased = 25), name = NULL) +
  scale_fill_manual(values = c(Increased = COL$up, Decreased = COL$down), name = NULL) +
  labs(x = "Number of DEGs (log scale)", y = NULL, tag = "B",
       caption = paste0("Direction refers to: sex-shared, both sexes; female-specific, women;\n",
                        "male-specific and sex-dimorphic, men. Filter applied to PFC only (unfiltered\n",
                        "HPC list not retained).")) +
  guides(shape = guide_legend(override.aes = list(size = 2.4, colour = NA))) +
  theme_scm() + TAGTH +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = FS$base, colour = COL$ink2),
        strip.text = element_text(size = FS$base, face = "plain", hjust = 0.5, lineheight = 1.05,
                                  margin = margin(b = 2)),
        panel.spacing.x = unit(2.4, "mm"),
        legend.position = "top", legend.justification = c(0, 0.5), legend.location = "plot",
        legend.text = element_text(size = FS$small + 0.3), legend.key.spacing.x = unit(5, "mm"),
        legend.margin = margin(0, 0, 5, 28), legend.box.spacing = unit(0, "mm"),
        plot.caption = element_text(size = FS$tick, colour = COL$ink2, hjust = 0, lineheight = 1.15,
                                    margin = margin(t = 7)),
        plot.caption.position = "plot",
        plot.margin = margin(2, 4, 2, 1))

# =============================================================================
# C  UpSet plot of exclusive intersections among the 8 category sets
# =============================================================================
up <- read.csv("inputs/source_data/Figure1C.csv")[, c("intersection", "size")]
set_names <- c(paste("HPC", CAT_ORDER), paste("PFC", CAT_ORDER))
memb <- lapply(strsplit(up$intersection, " & ", fixed = TRUE), trimws)
stopifnot(all(unlist(memb) %in% set_names))
up <- up %>% mutate(ix = row_number(), nset = lengths(memb)) %>% arrange(desc(size), nset) %>% mutate(ix = row_number())
memb <- lapply(strsplit(up$intersection, " & ", fixed = TRUE), trimws)
set_size <- sapply(set_names, function(s) sum(up$size[sapply(memb, function(m) s %in% m)]))
chk("C HPC set sizes (shared/female/male/dimorphic)", set_size[1:4], c(253, 45, 9, 4))
chk("C PFC filtered set sizes", set_size[5:8], c(1128, 2, 5, 1))
core_ix <- which(up$intersection == "HPC sex-shared & PFC sex-shared")
chk("C cross-region core size", up$size[core_ix], 109)
ncol_u <- nrow(up)
up$core <- seq_len(ncol_u) == core_ix
XL <- c(-1.25, ncol_u + 0.6)
pCbar <- ggplot(up, aes(ix, size)) +
  geom_col(aes(fill = core), width = 0.62) +
  geom_text(aes(label = format(size, big.mark = ",", trim = TRUE)), vjust = -0.45, size = pt(6),
            family = FONT, colour = COL$ink) +
  annotate("segment", x = core_ix + 0.15, y = up$size[core_ix] + 95, xend = core_ix + 1.75, yend = 520,
           colour = COL$ink2, linewidth = 0.3) +
  annotate("text", x = core_ix + 1.85, y = 545, label = "cross-region core", hjust = 0, vjust = 0,
           size = pt(6.2), family = FONT, colour = COL$ink) +
  scale_fill_manual(values = c(`TRUE` = COL$ink, `FALSE` = COL$ink2), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1150), breaks = c(0, 500, 1000), labels = c("0", "500", "1,000"),
                     expand = c(0, 0)) +
  annotate("text", x = XL[1] - 2.6, y = 575, label = "Genes in exclusive\nintersection", angle = 90, size = pt(FS$base), family = FONT, colour = COL$ink, lineheight = 1.05) +
  coord_cartesian(xlim = XL, clip = "off") +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_scm() + TAGTH +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(lineheight = 1.05), plot.margin = margin(8, 4, 0, 0))
mat <- expand.grid(ix = seq_len(ncol_u), r = seq_along(set_names)) %>%
  mutate(on = mapply(function(i, r) set_names[r] %in% memb[[i]], ix, r))
lines_m <- mat %>% filter(on) %>% group_by(ix) %>% summarise(y0 = min(r), y1 = max(r), .groups = "drop") %>%
  filter(y1 > y0)
row_lab <- sprintf("%s (%s)", set_names, format(set_size, big.mark = ",", trim = TRUE))
sw <- data.frame(r = seq_along(set_names), cat = rep(CAT_ORDER, 2))
pCmat <- ggplot() +
  annotate("rect", xmin = -11.5, xmax = XL[2], ymin = 4.5, ymax = 8.5, fill = "#f4f3f0") +
  geom_point(data = mat, aes(ix, r), colour = COL$grid, size = 1.15) +
  geom_segment(data = lines_m, aes(x = ix, xend = ix, y = y0, yend = y1), colour = COL$ink, linewidth = 0.6) +
  geom_point(data = filter(mat, on), aes(ix, r), colour = COL$ink, size = 1.35) +
  geom_tile(data = sw, aes(x = -0.55, y = r, fill = cat), width = 0.42, height = 0.36) +
  scale_fill_manual(values = CAT_COLS, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(breaks = seq_along(set_names), labels = row_lab, limits = c(8.5, 0.5), expand = c(0, 0)) +
  coord_cartesian(xlim = XL, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_scm() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = FS$small + 0.3, colour = COL$ink, margin = margin(r = 0)),
        plot.margin = margin(1, 4, 2, 0))

# =============================================================================
# D  heatmap of the 109 core genes
# =============================================================================
core <- read.csv("inputs/source_data/Figure1D.csv", stringsAsFactors = FALSE)
GROUPS <- c("Chemokines/cytokines & receptors" = "Chemokines, cytokines & receptors",
            "Alarmins/complement/immune receptors & signalling" = "Alarmins & immune signalling",
            "AP-1/IEG transcription factors" = "AP-1 immediate-early TFs",
            "EGR family" = "EGR family", "NR4A family" = "NR4A family",
            "Other transcriptional/chromatin regulators" = "Other transcriptional regulators",
            "Apoptosis (BCL2 family) & autophagy" = "Apoptosis (BCL2) & autophagy",
            "Lipid / proteoglycan / glycan metabolism" = "Lipid, proteoglycan & glycan",
            "Cytoskeleton / small-GTPase regulation / motility" = "Cytoskeleton, GTPases & motility",
            "Vesicle trafficking / endolysosome / nuclear transport" = "Vesicle trafficking & endolysosome",
            "Solute/amino-acid/ion transport" = "Solute & ion transport",
            "Proteostasis / translation / DNA repair / other" = "Proteostasis, translation & other")
stopifnot(all(core$functional_group %in% names(GROUPS)))
CONTR <- c("HPC_men", "HPC_women", "PFC_men", "PFC_women")
core <- core %>% mutate(grank = match(functional_group, names(GROUPS)),
                        sumFC = rowSums(across(all_of(CONTR)))) %>%
  arrange(grank, desc(sumFC)) %>%
  mutate(high = stratum == "high", marker = !is.na(state_marker) & state_marker != "",
         fb = toupper(as.character(feedback_regulator)) == "TRUE")
sgn <- sign(as.matrix(core[, CONTR]))
chk("D core genes", nrow(core), 109)
chk("D high effect (square)", sum(core$high), 38)
chk("D state marker MG8/MG10 (circle)", sum(core$marker), 7)
chk("D feedback regulators (diamond)", sum(core$fb), 7)
chk("D same direction in all 4 contrasts", sum(apply(sgn, 1, function(s) length(unique(s)) == 1)), 109)
chk("D increased / decreased", c(sum(sgn[, 1] > 0), sum(sgn[, 1] < 0)), c(101, 8))
fbset <- c("NR4A1", "NR4A3", "BCL6", "BCL3", "SOCS3", "DUSP1", "ZFP36L1")
CHECKS[["D feedback gene identity"]] <- sprintf("D feedback genes = {%s} -> %s",
  paste(core$gene[core$fb], collapse = ", "), if (setequal(core$gene[core$fb], fbset)) "PASS" else "MISMATCH")
CHECKS[["D JUN not marked"]] <- sprintf("D JUN state-marker flag = %s -> %s",
  core$marker[core$gene == "JUN"], if (!core$marker[core$gene == "JUN"]) "PASS" else "MISMATCH")
hi_ok <- all(core$high == apply(abs(as.matrix(core[, CONTR])) > 1.5, 1, all))
CHECKS[["D high = |log2FC|>1.5 in all four"]] <- sprintf("D high-effect flag equals |log2FC|>1.5 in all 4 contrasts -> %s",
                                                     if (hi_ok) "PASS" else "MISMATCH")
CLIP <- 4
# row layout: each group = header line + gene lines
blocks <- list(names(GROUPS)[1:8], names(GROUPS)[9:12])
lay <- lapply(seq_along(blocks), function(b) {
  y <- 0; hd <- list(); gn <- list()
  for (g in blocks[[b]]) {
    hd[[g]] <- data.frame(y = y, lab = GROUPS[[g]]); y <- y + 1
    sub <- core[core$functional_group == g, ]
    gn[[g]] <- data.frame(sub, y = y + seq_len(nrow(sub)) - 1); y <- y + nrow(sub)
  }
  list(head = do.call(rbind, hd), genes = do.call(rbind, gn), n = y)
})
NL <- max(sapply(lay, `[[`, "n"))
# x geometry (mm-like units)
LBL <- 11.6; CW <- 4.3; CG <- 1.0
xcell <- LBL + 0.5 + c(0, CW, 2 * CW + CG, 3 * CW + CG) + CW / 2
xg <- xcell[4] + CW / 2 + c(2.4, 4.9, 7.4)
XMAX <- xg[3] + 1.6
YTOP <- -7.6
hdr_y <- -0.7
make_block <- function(b, legend = FALSE) {
  L <- lay[[b]]
  cells <- L$genes %>% select(gene, y, all_of(CONTR)) %>%
    pivot_longer(all_of(CONTR), names_to = "contr", values_to = "v") %>%
    mutate(x = xcell[match(contr, CONTR)])
  big <- filter(cells, abs(v) > CLIP) %>% mutate(lab = fmt_minus(sprintf("%.1f", v)))
  sym <- bind_rows(
    L$genes %>% filter(high) %>% transmute(x = xg[1], y, s = "high"),
    L$genes %>% filter(marker) %>% transmute(x = xg[2], y, s = "marker"),
    L$genes %>% filter(fb) %>% transmute(x = xg[3], y, s = "fb"))
  colh <- data.frame(x = xcell, lab = c("M", "F", "M", "F"))
  regh <- data.frame(x0 = xcell[c(1, 3)] - CW / 2 + 0.15, x1 = xcell[c(2, 4)] + CW / 2 - 0.15,
                     lab = c("HPC", "PFC"))
  p <- ggplot() +
    geom_tile(data = cells, aes(x, y, fill = v), width = CW - 0.3, height = 0.88) +
    geom_text(data = big, aes(x, y, label = lab), colour = "white", size = pt(6), family = FONT) +
    geom_text(data = L$head, aes(0, y + 0.08, label = lab), hjust = 0, fontface = "bold", size = pt(6),
              family = FONT, colour = COL$ink) +
    geom_text(data = L$genes, aes(LBL, y, label = gene), hjust = 1, fontface = "italic", size = pt(6),
              family = FONT, colour = COL$ink) +
    geom_point(data = sym, aes(x, y, shape = s, size = s), colour = COL$ink) +
    geom_text(data = colh, aes(x, hdr_y - 0.2, label = lab), size = pt(6.5), family = FONT, colour = COL$ink) +
    geom_segment(data = regh, aes(x = x0, xend = x1, y = hdr_y - 1.55, yend = hdr_y - 1.55),
                 colour = COL$ink2, linewidth = 0.3) +
    geom_text(data = regh, aes((x0 + x1) / 2, hdr_y - 2.45, label = lab), fontface = "bold",
              size = pt(7), family = FONT, colour = COL$ink) +
    annotate("text", x = xg, y = hdr_y + 0.3, label = c("High effect", "State marker", "Feedback reg."),
             angle = 90, hjust = 0, size = pt(6), family = FONT, colour = COL$ink2) +
    scale_shape_manual(values = c(high = 15, marker = 16, fb = 18), guide = "none") +
    scale_size_manual(values = c(high = 1.25, marker = 1.35, fb = 1.75), guide = "none") +
    scale_fill_lfc(guide = "none") +
    scale_x_continuous(limits = c(0, XMAX), expand = c(0, 0)) +
    scale_y_reverse(limits = c(NL - 0.5, YTOP), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void_scm() + theme(plot.margin = margin(0, 1.5, 0, 1.5))
  if (legend) {
    y0 <- L$n + 0.55
    grad <- data.frame(v = seq(-CLIP, CLIP, length.out = 120))
    gx0 <- 2.0; gx1 <- XMAX - 4.0
    grad$x <- gx0 + (grad$v + CLIP) / (2 * CLIP) * (gx1 - gx0)
    ticks <- data.frame(v = c(-4, -2, 0, 2, 4), lab = c("\u2264\u22124", "\u22122", "0", "2", "\u22654")) %>%
      mutate(x = gx0 + (v + CLIP) / (2 * CLIP) * (gx1 - gx0))
    keys <- data.frame(s = c("high", "marker", "fb"), y = y0 + 4.85 + 0:2 * 1.05,
                       lab = c("'High effect (|'*log[2]*'FC| > 1.5)'", "'State marker (MG8/MG10 list)'",
                               "'Negative-feedback regulator'"))
    p <- p +
      annotate("text", x = gx0 - 0.2, y = y0, hjust = 0, size = pt(6), family = FONT, colour = COL$ink,
               label = "log[2]*'FC, inflammatory-AD vs'", parse = TRUE) +
      annotate("text", x = gx0 - 0.2, y = y0 + 1.05, hjust = 0, size = pt(6), family = FONT, colour = COL$ink,
               label = "homeostatic non-AD") +
      geom_tile(data = grad, aes(x, y0 + 2.25, fill = v), width = (gx1 - gx0) / 119 * 1.02, height = 0.85) +
      geom_segment(data = ticks, aes(x = x, xend = x, y = y0 + 2.7, yend = y0 + 3.0), colour = COL$ink2,
                   linewidth = 0.3) +
      geom_text(data = ticks, aes(x, y0 + 3.55, label = lab), size = pt(6), family = FONT, colour = COL$ink2) +
      geom_point(data = keys, aes(gx0 + 0.6, y, shape = s, size = s), colour = COL$ink) +
      geom_text(data = keys, aes(gx0 + 2.0, y, label = lab), hjust = 0, parse = TRUE, size = pt(6),
                family = FONT, colour = COL$ink2) +
      annotate("text", x = gx0 - 0.2, y = y0 + 7.75, hjust = 0, size = pt(6), family = FONT, colour = COL$ink2,
               label = "M, men; F, women. Colors clipped at \u00b14;\nvalues beyond \u00b14 are printed.", vjust = 1, lineheight = 1.05)
  }
  p
}
pD1 <- make_block(1) + labs(tag = "D") + TAGTH +
  theme(plot.tag.position = c(-0.04, 1), plot.margin = margin(0, 1.5, 0, 3))
pD2 <- make_block(2, legend = TRUE)

# =============================================================================
# assemble
# =============================================================================
left <- free(pB) / pCbar / pCmat + plot_layout(heights = c(60, 46, 44))
right <- (pD1 | pD2) + plot_layout(widths = c(1, 1))
bottom <- (left | right) + plot_layout(widths = c(87, 96))
fig <- free(pA) / bottom + plot_layout(heights = c(AH, 166)) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 3, 1)))
HEIGHT <- 226
save_fig(fig, "Figure1", width_mm = 183, height_mm = HEIGHT, outdir = "out")
if (dev.cur() > 1) invisible(dev.off()); invisible(gc())
writeLines(unlist(CHECKS), "out/Figure1_checks.txt")
cat(unlist(CHECKS), sep = "\n")
