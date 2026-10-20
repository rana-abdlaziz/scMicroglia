# Supplementary Figure S3 - sex bias of a priori gene sets and of genes with a sex difference > 1 log2 unit
# Run from the working directory containing theme_scm.R and inputs/ (unzipped scm_figure_inputs.zip).
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggtext) })
source("theme_scm.R")
SUP <- "inputs/supp_csv"

fmt_p_md <- function(p) {                      # 0.49 / 0.049 / 1.51 x 10^-4
  if (p >= 0.1) return(sprintf("%.2f", p))
  if (p >= 0.01) return(sprintf("%.3f", p))
  e <- floor(log10(p)); sprintf("%.2f \u00d7 10<sup>%s</sup>", p / 10^e, fmt_minus(as.character(e)))
}
BIAS_LAB <- expression("|" * log[2] * FC[F] * "| \u2212 |" * log[2] * FC[M] * "|")
th_md_title <- theme(plot.title = element_markdown(size = FS$title, face = "bold", hjust = 0, lineheight = 1.15,
                                                   margin = margin(b = 4, l = 11)))

# =========================================================================================
# A  a priori gene sets (Block F genes, summary statistics from Table S5 Block A)
# =========================================================================================
f <- fread(file.path(SUP, "S5_F_program_genes.csv"))
s5 <- fread(file.path(SUP, "S5.csv"))
groups <- data.table(program = c("chemokine/cytokine", "IEG/AP-1/EGR/NR4A", "chemokine/cytokine", "IEG/AP-1/EGR/NR4A"),
                     region = c("HPC", "HPC", "PFC", "PFC"),
                     set_key = c("chemokine_cytokine_HPC", "IEG_AP1_EGR_NR4A_HPC", "chemokine_cytokine_PFC", "IEG_AP1_EGR_NR4A_PFC"),
                     xlab = c("chemokine/\ncytokine\nHPC", "IEG/AP-1/\nEGR/NR4A\nHPC", "chemokine/\ncytokine\nPFC",
                              "IEG/AP-1/\nEGR/NR4A\nPFC"))
groups[, i := .I]
pts <- rbindlist(lapply(seq_len(nrow(groups)), function(k) {
  gk <- groups[k]; col <- paste0("sexbias_", gk$region, "_absF_minus_absM")
  d <- f[program == gk$program & !is.na(get(col)), .(gene, bias = get(col))]
  d[, i := gk$i]
}))
set.seed(2); pts[, xj := i + runif(.N, -0.15, 0.15)]
chk <- pts[, .(n_pts = .N, median_pts = median(bias), P_recomputed = wilcox.test(bias, exact = TRUE)$p.value), by = i]
hdr <- merge(groups, chk, by = "i")
tab <- rbindlist(lapply(hdr$set_key, function(kk) s5[startsWith(s5$set_id, kk), .(n_tab = n, med_tab = median_absF_minus_absM, P_tab = Wilcoxon_P)]))
stopifnot(nrow(tab) == nrow(hdr)); hdr <- cbind(hdr, tab)
stopifnot(all(hdr$n_pts == hdr$n_tab), all(abs(hdr$median_pts - hdr$med_tab) < 1e-9),
          all(abs(hdr$P_recomputed - hdr$P_tab) < 1e-6))
hdr[, lab := sprintf("n = %d<br>median %s<br>*P* = %s", n_tab, fmt_minus(sprintf("%+.2f", med_tab)),
                     sapply(P_tab, fmt_p_md))]
YA <- c(-2.6, 4.95)
pA <- ggplot(pts) +
  geom_hline(yintercept = 0, colour = COL$muted, linewidth = 0.3) +
  geom_point(aes(xj, bias), shape = 16, size = 1.3, colour = COL$ink2, alpha = 0.85) +
  geom_segment(data = hdr, aes(x = i - 0.28, xend = i + 0.28, y = med_tab, yend = med_tab), colour = COL$ink,
               linewidth = 0.8, lineend = "butt") +
  geom_richtext(data = hdr, aes(x = i, y = YA[2] - 0.05, label = lab), vjust = 1, size = pt(FS$small), family = FONT,
                fill = NA, label.colour = NA, label.padding = unit(0, "pt"), lineheight = 1.15, colour = COL$ink) +
  annotate("text", x = 4.45, y = 2.3, label = "larger in women", hjust = 1, vjust = 1, size = pt(FS$small),
           family = FONT, colour = COL$women) +
  annotate("text", x = 4.45, y = -2.5, label = "larger in men", hjust = 1, vjust = 0, size = pt(FS$small),
           family = FONT, colour = COL$men) +
  scale_x_continuous(breaks = groups$i, labels = groups$xlab) +
  scale_y_continuous(breaks = -2:4, labels = lab_minus) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = YA, expand = FALSE) +
  labs(x = NULL, y = BIAS_LAB,
       title = "A priori gene sets: only the PFC<br>chemokine/cytokine set is sex-biased") +
  theme_scm() + th_md_title +
  theme(axis.text.x = element_text(size = FS$tick, colour = COL$ink, lineheight = 0.95), axis.ticks.x = element_blank())

# =========================================================================================
# B  genes with |log2FC_F| - |log2FC_M| > 1 (Table S5 Block G)
# =========================================================================================
g <- fread(file.path(SUP, "S5_G_sexdiff_gt1.csv"))
stopifnot(all(abs(g$sexbias_absF_minus_absM) > 1), all(g$effect_stratum == "high"),
          all(abs(abs(g$log2FC_women) - abs(g$log2FC_men) - g$sexbias_absF_minus_absM) < 1e-6))
setorder(g, region, -sexbias_absF_minus_absM)
g[, `:=`(y = rev(seq_len(.N)), col = fifelse(sexbias_absF_minus_absM > 0, COL$women, COL$men),
         lab = sprintf("<i>%s</i> (%s)", gene, region),
         vals = fmt_minus(sprintf("%.2f / %.2f", log2FC_men, log2FC_women)))]
nB <- g[, .N, by = region]
VX <- 3.65
pB <- ggplot(g) +
  geom_vline(xintercept = 0, colour = COL$muted, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = sexbias_absF_minus_absM, y = y, yend = y, colour = col), linewidth = 0.6) +
  geom_point(aes(sexbias_absF_minus_absM, y, colour = col), size = 1.7) +
  geom_text(aes(VX, y, label = vals), hjust = 0, size = pt(FS$small), family = FONT, colour = COL$ink2) +
  annotate("richtext", x = VX, y = nrow(g) + 0.75, label = "**log<sub>2</sub>FC men / women**", hjust = 0, vjust = 0.5,
           size = pt(FS$small), family = FONT, fill = NA, label.colour = NA, label.padding = unit(0, "pt")) +
  scale_colour_identity() +
  scale_y_continuous(breaks = g$y, labels = g$lab) +
  scale_x_continuous(breaks = -2:5, labels = lab_minus) +
  coord_cartesian(xlim = c(-2.6, 6.0), ylim = c(0.4, nrow(g) + 1.2), expand = FALSE, clip = "off") +
  labs(x = BIAS_LAB, y = NULL,
       title = "Genes with a sex difference &gt; 1 log<sub>2</sub> unit<br>(high-effect shared genes)") +
  theme_scm() + th_md_title +
  theme(axis.text.y = element_markdown(size = FS$tick, colour = COL$ink), plot.margin = margin(3, 22, 3, 3))

fig <- add_tags(pA + free(pB, type = "label", side = "b") + plot_layout(widths = c(1, 0.92)))
dir.create("out", showWarnings = FALSE)
save_fig(fig, "Supplementary_Figure_S3", width_mm = 183, height_mm = 88)
invisible(gc())

# ---- numeric checks against the SI legend ------------------------------------------------
print(hdr[, .(region, program, n_pts, n_tab, median = round(med_tab, 3), P_tab = signif(P_tab, 3),
              P_recomputed = signif(P_recomputed, 3))])
print(nB)
