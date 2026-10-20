# Supplementary Figure S1 - sensitivity of the cross-region core to the PFC effect-size filter.
# Self-contained: run from the working directory containing inputs/ and theme_scm.R
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr) })
source("theme_scm.R")
CHECKS <- character()
chk <- function(name, got, expected) {
  ok <- isTRUE(all.equal(as.numeric(got), as.numeric(expected)))
  CHECKS <<- c(CHECKS, sprintf("%s: got %s, legend %s -> %s", name, paste(got, collapse = "/"),
                               paste(expected, collapse = "/"), if (ok) "PASS" else "MISMATCH"))
}
RET <- COL$ink2   # retained by the PFC filter (filtered core)
REM <- COL$pfc    # removed by the PFC filter (orange)

s4 <- read.csv("inputs/supp_csv/S4.csv", stringsAsFactors = FALSE) %>%
  mutate(HPC_mean = (log2FC_HPC_men + log2FC_HPC_women) / 2,
         PFC_mean = (log2FC_PFC_men + log2FC_PFC_women) / 2,
         status = ifelse(in_filtered_core_109 == "yes", "retained", "removed"))
n_ret <- sum(s4$status == "retained"); n_rem <- sum(s4$status == "removed")
chk("genes in panel A", nrow(s4), 163)
chk("retained (filtered core)", n_ret, 109)
chk("removed (intermediate PFC effect)", n_rem, 54)
conc <- sign(s4$HPC_mean) == sign(s4$PFC_mean)
chk("direction concordant HPC vs PFC (sign of means)", sum(conc), 162)
chk("cross_region_direction == concordant (table column)", sum(s4$cross_region_direction == "concordant"), 162)
CHECKS <- c(CHECKS, sprintf("only discordant gene = %s (legend: C22orf34) -> %s", paste(s4$gene[!conc], collapse = ","),
                            if (identical(s4$gene[!conc], "C22orf34")) "PASS" else "MISMATCH"))
lab_genes <- data.frame(gene = c("C22orf34", "RUNX1", "FOXP1", "CDKN1A", "RELB"),
                        tx = c(-1.0, 5.2, 5.2, 6.4, 6.4), ty = c(2.6, 0.4, 1.3, 2.6, 3.6),
                        hj = c(1, 0, 0, 0, 0)) %>%
  left_join(s4[, c("gene", "HPC_mean", "PFC_mean", "status")], by = "gene")
stopifnot(!anyNA(lab_genes$HPC_mean))
CHECKS <- c(CHECKS, sprintf("labelled RUNX1/FOXP1/CDKN1A/RELB are removed genes -> %s",
                            if (all(lab_genes$status[-1] == "removed")) "PASS" else "MISMATCH"))

leg_lab <- c(retained = sprintf("retained by PFC filter (filtered core, n = %d)", n_ret),
             removed = sprintf("removed by PFC filter (intermediate PFC effect, n = %d)", n_rem))
s4$status <- factor(s4$status, c("retained", "removed"))
rx <- range(c(s4$HPC_mean, -3, 10.5)); ry <- range(c(s4$PFC_mean, -3, 10.5))
pA <- ggplot(s4, aes(HPC_mean, PFC_mean)) +
  geom_hline(yintercept = 0, colour = "#cccccc", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#cccccc", linewidth = 0.3) +
  annotate("segment", x = -3, y = -3, xend = 10.5, yend = 10.5, colour = COL$muted, linewidth = 0.3,
           linetype = "22") +
  geom_point(aes(shape = status, colour = status, alpha = status), size = 1.15, stroke = 0.45) +
  geom_segment(data = lab_genes, aes(x = HPC_mean, y = PFC_mean, xend = tx - 0.12 * (1 - 2 * hj), yend = ty),
               colour = "#555555", linewidth = 0.25) +
  geom_text(data = lab_genes, aes(tx, ty, label = gene, hjust = hj), fontface = "italic", size = pt(6),
            family = FONT, colour = COL$ink) +
  scale_shape_manual(values = c(retained = 16, removed = 1), labels = leg_lab, name = NULL) +
  scale_colour_manual(values = c(retained = RET, removed = REM), labels = leg_lab, name = NULL) +
  scale_alpha_manual(values = c(retained = 0.8, removed = 1), labels = leg_lab, name = NULL) +
  scale_x_continuous(breaks = seq(-2, 10, 2), labels = lab_minus, expand = expansion(mult = 0.04)) +
  scale_y_continuous(breaks = seq(-2, 10, 2), labels = lab_minus, expand = expansion(mult = 0.04)) +
  coord_cartesian(xlim = c(-3.4, 10.9), ylim = c(-3.4, 10.9)) +
  labs(x = expression("mean " * log[2] * "FC in HPC (men, women)"),
       y = expression("mean " * log[2] * "FC in PFC (men, women)"),
       title = "Direction concordant for 162 of 163 genes\nwhen the PFC filter is relaxed") +
  guides(shape = guide_legend(override.aes = list(size = 1.4))) +
  theme_scm() +
  theme(legend.position = "inside", legend.position.inside = c(0.0, 1.0), legend.justification = c(0, 1),
        legend.key.size = unit(6, "pt"), legend.key.spacing.y = unit(0.5, "pt"),
        legend.text = element_text(size = FS$tick, margin = margin(l = 1)),
        plot.title = element_text(lineheight = 1.05))

# B: |log2FC| in PFC per sex, removed vs retained genes
bl <- s4 %>% select(gene, status, log2FC_PFC_men, log2FC_PFC_women) %>%
  pivot_longer(c(log2FC_PFC_men, log2FC_PFC_women), names_to = "sex", values_to = "lfc") %>%
  mutate(sex = ifelse(grepl("women", sex), "women", "men"), v = abs(lfc),
         grp = factor(paste(status, sex, sep = "\n"),
                      c("removed\nmen", "removed\nwomen", "retained\nmen", "retained\nwomen")),
         gx = as.integer(grp))
set.seed(1)
bl$x <- bl$gx + runif(nrow(bl), -0.15, 0.15)
med <- bl %>% group_by(grp, gx, status) %>% summarise(m = median(v), n = n(), .groups = "drop")
CHECKS <- c(CHECKS, sprintf("B medians |log2FC| PFC: %s (no legend values; reported)",
                            paste(sprintf("%s=%.3f (n=%d)", gsub("\n", " ", med$grp), med$m, med$n), collapse = "; ")))
rem_in_band <- with(filter(bl, status == "removed"), all(v > 0))
pB <- ggplot(bl, aes(x, v)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = "#f0efec") +
  geom_point(aes(colour = status), size = 0.75, alpha = 0.75, stroke = 0) +
  geom_segment(data = med, aes(x = gx - 0.25, xend = gx + 0.25, y = m, yend = m), linewidth = 0.7,
               colour = COL$ink) +
  annotate("text", x = 4.45, y = 1.0, label = "intermediate\nstratum\n(0.5\u20131.5)", hjust = 0, size = pt(6),
           family = FONT, colour = COL$ink2, lineheight = 1.05) +
  scale_colour_manual(values = c(retained = RET, removed = REM), guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = levels(bl$grp), limits = c(0.5, 5.3), expand = c(0, 0)) +
  scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10), labels = c("0.1", "0.3", "1", "3", "10"),
                minor_breaks = c(seq(0.1, 0.9, 0.1), 2:9), guide = guide_axis(minor.ticks = TRUE)) +
  coord_cartesian(ylim = c(0.09, 10), clip = "off") +
  labs(x = NULL, y = expression("|" * log[2] * "FC| in PFC"),
       title = "Removed genes have intermediate\nPFC effects (0.5\u20131.5) by construction") +
  theme_scm() +
  theme(axis.text.x = element_text(lineheight = 1.0, colour = COL$ink2), axis.minor.ticks.length.y = rel(0.55),

        plot.title = element_text(lineheight = 1.05))

fig <- add_tags((pA | pB) + plot_layout(widths = c(1.05, 1)))
save_fig(fig, "Supplementary_Figure_S1", width_mm = 183, height_mm = 88, outdir = "out")
if (dev.cur() > 1) invisible(dev.off()); invisible(gc())
writeLines(CHECKS, "out/Supplementary_Figure_S1_checks.txt")
cat(CHECKS, sep = "\n")
