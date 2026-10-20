# Supplementary Figure S4 - gseGO of sex-shared genes (R / ggplot2 redraw of build_supp_figs_S1-S5.py, S4 block).
# The 12 GO-BP (A HPC, C PFC) and GO-MF (B HPC, D PFC) terms with the lowest nominal P; bars = NES,
# labels = nominal fgsea P. Source: Supplementary Table S6 (gseGO rows, sex-shared lists).
# Run from a fresh R session in the working directory that contains inputs/ and theme_scm.R.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
if (!file.exists("inputs/supp_csv/S6.csv")) stop("inputs/supp_csv/S6.csv not found - unzip the figure inputs first")
source("theme_scm.R")
tmp_dev <- tempfile(fileext = ".png"); ragg::agg_png(tmp_dev)   # font-metric device for grob construction (avoids the Arial-less default pdf device)

s6  <- fread("inputs/supp_csv/S6.csv")
gse <- s6[analysis == "gseGO" & category_list == "sex-shared"]

# same term abbreviations as Figure 3, then truncate with an ellipsis
ABBR <- list(c("RNA polymerase II-specific", "RNA Pol II-specific"), c("DNA-binding transcription factor", "DNA-binding TF"),
             c("transcription factor", "TF"), c("positive regulation of", "pos. reg. of"),
             c("negative regulation of", "neg. reg. of"), c("regulation of", "reg. of"),
             c("involved in", "in"), c("establishment of", "estab. of"), c(" signaling pathway", " signaling"),
             c("structural constituent", "struct. constituent"), c("transmembrane transporter", "transm. transporter"))
short <- function(s, n = 46) vapply(s, function(x) {
  for (a in ABBR) x <- gsub(a[1], a[2], x, fixed = TRUE)
  x <- trimws(x)
  if (nchar(x) <= n) x else paste0(sub("\\s+$", "", substr(x, 1, n - 1)), "\u2026")
}, character(1), USE.NAMES = FALSE)
# P label, always a x 10^b (2 significant digits), as a plotmath string
p_lab <- function(p) vapply(p, function(v) {
  e <- floor(log10(v)); m <- v / 10^e
  if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("italic(P)~'='~'%s' %%*%% 10^%d", formatC(m, format = "f", digits = 1), e)
}, character(1))

# ---- legend checks: term counts at nominal P < 0.05 and positive NES ------------------------
expect <- data.table(region = c("HPC", "HPC", "PFC", "PFC"), database = c("GO-BP", "GO-MF", "GO-BP", "GO-MF"),
                     n_exp = c(308, 30, 732, 60), pos_exp = c(301, 29, 663, 53))
cnt <- gse[, .(n = .N, pos = sum(NES > 0), maxP = max(P_value)), by = .(region, database)]
cnt <- merge(expect, cnt, by = c("region", "database"))
cnt[, pass := n == n_exp & pos == pos_exp & maxP < 0.05]
print(cnt)

XMAX <- max(gse[, .SD[order(P_value)][1:12], by = .(region, database)]$NES) + 1.25
panel <- function(reg, db) {
  all <- gse[region == reg & database == db]
  d <- all[order(P_value)][1:12]
  d[, y := rev(seq_len(.N))]
  d[, fill := ifelse(NES > 0, COL$shared, SEQ_COLS[2])]
  ggplot(d) +
    geom_vline(xintercept = 0, colour = COL$muted, linewidth = 0.3) +
    geom_col(aes(NES, y, fill = I(fill)), orientation = "y", width = 0.66) +
    geom_text(aes(NES + ifelse(NES > 0, 0.05, -0.05), y, label = p_lab(P_value), hjust = ifelse(NES > 0, 0, 1)),
              parse = TRUE, size = pt(FS$tick), colour = COL$ink2, family = FONT) +
    scale_y_continuous(breaks = d$y, labels = short(d$Description), expand = expansion(add = 0.6)) +
    scale_x_continuous(breaks = 0:3, limits = c(min(-0.3, min(d$NES) - 0.9), XMAX), expand = c(0, 0)) +
    labs(x = "Normalized enrichment score (NES)", y = NULL,
         title = sprintf("%s sex-shared genes, %s", reg, db),
         subtitle = parse(text = sprintf("'%d terms at nominal'~italic(P)~'< 0.05 (%d with NES > 0)'", nrow(all), sum(all$NES > 0)))) +
    theme_scm() +
    theme(axis.text.y = element_text(size = FS$small, colour = COL$ink, margin = margin(r = 2)),
          axis.ticks.y = element_blank(), axis.title.x = element_text(size = FS$small, margin = margin(t = 2)),
          plot.title = element_text(margin = margin(b = 1.5, l = 11)),
          plot.subtitle = element_text(margin = margin(b = 4, l = 11)),
          plot.margin = margin(2, 4, 6, 2))
}
pA <- panel("HPC", "GO-BP"); pB <- panel("HPC", "GO-MF"); pC <- panel("PFC", "GO-BP"); pD <- panel("PFC", "GO-MF")
fig <- add_tags(wrap_plots(pA, pB, pC, pD, ncol = 2, byrow = FALSE))
save_fig(fig, "Supplementary_Figure_S4", width_mm = 183, height_mm = 122, outdir = "out")
invisible(gc())
invisible(dev.off()); unlink(tmp_dev)
