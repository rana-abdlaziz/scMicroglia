# Supplementary Figure S2 - overlap of DEGs with the curated state-marker lists (circularity check).
# Self-contained: run from the working directory containing inputs/ and theme_scm.R
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr) })
source("theme_scm.R")
CHECKS <- character()
chk <- function(name, got, expected) {
  ok <- isTRUE(all.equal(as.numeric(got), as.numeric(expected)))
  CHECKS <<- c(CHECKS, sprintf("%s: got %s, legend %s -> %s", name, paste(got, collapse = "/"),
                               paste(expected, collapse = "/"), if (ok) "PASS" else "MISMATCH"))
}
s2 <- read.csv("inputs/supp_csv/S2.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))
s1 <- read.csv("inputs/supp_csv/S1.csv", stringsAsFactors = FALSE)
g  <- read.csv("inputs/supp_csv/S5_G_sexdiff_gt1.csv", stringsAsFactors = FALSE)

# ---- A: marker overlap of sex-shared lists ----------------------------------------------------
sets <- data.frame(region = c("HPC", "PFC", "PFC"), list = c("filtered", "filtered", "unfiltered"),
                   lab0 = c("HPC\nsex-shared", "PFC filtered\nsex-shared", "PFC unfiltered\nsex-shared"))
cnt <- sets %>% rowwise() %>% mutate(
  d = list(s2[s2$region == region & s2$list == list & s2$category == "sex-shared", ]),
  n = nrow(d),
  infl = sum(d$marker_of_contrast_group == "Inflammatory", na.rm = TRUE),
  homeo = sum(d$marker_of_contrast_group == "Homeostatic", na.rm = TRUE),
  other = sum(!is.na(d$state_marker_lists)) - infl - homeo) %>% ungroup() %>% select(-d) %>%
  mutate(lab = sprintf("%s\n(%s)", lab0, format(n, big.mark = ",", trim = TRUE)), x = row_number())
chk("A list sizes HPC / PFC filtered / PFC unfiltered", cnt$n, c(253, 1128, 1935))
chk("A HPC infl/homeo/other", unlist(cnt[1, c("infl", "homeo", "other")]), c(11, 1, 7))
chk("A PFC filtered infl/homeo/other", unlist(cnt[2, c("infl", "homeo", "other")]), c(9, 0, 17))
chk("A PFC unfiltered infl/homeo/other", unlist(cnt[3, c("infl", "homeo", "other")]), c(23, 7, 41))
nonshared_markers <- sum(!is.na(s2$state_marker_lists[s2$list == "filtered" & s2$category != "sex-shared"]))
chk("A sex-specific/dimorphic genes in filtered lists that are markers", nonshared_markers, 0)
# marker lists in S2 should come from the S1 curated lists
mk_states <- unique(trimws(unlist(strsplit(na.omit(s2$state_marker_lists), "[;,]"))))
CHECKS <- c(CHECKS, sprintf("A marker-list labels in S2 all present in S1 state lists -> %s",
                            if (all(mk_states %in% s1$state)) "PASS" else "MISMATCH"))
stk <- cnt %>% select(x, lab, infl, homeo, other) %>% pivot_longer(c(infl, homeo, other), names_to = "k", values_to = "v") %>%
  mutate(k = factor(k, c("other", "homeo", "infl")))
KEY <- c(infl = "marker of MG2/MG8/MG10 (inflammatory, contrast group)",
         homeo = "marker of MG0 (homeostatic, reference group)",
         other = "marker of another state list")
tot <- cnt %>% mutate(t = infl + homeo + other, txt = sprintf("%d / %d / %d", infl, homeo, other))
ymax <- max(tot$t) * 1.35
pA <- ggplot(stk, aes(x, v, fill = k)) +
  geom_col(width = 0.6) +
  geom_text(data = tot, aes(x, t + 1.2, label = txt), inherit.aes = FALSE, vjust = 0, size = pt(6),
            family = FONT, colour = COL$ink) +
  scale_fill_manual(values = c(infl = COL$act, homeo = COL$inh, other = COL$uns), breaks = c("infl", "homeo", "other"),
                    labels = KEY, name = NULL) +
  scale_x_continuous(breaks = cnt$x, labels = cnt$lab, expand = expansion(add = 0.45)) +
  scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 80, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "sex-shared DEGs that are curated\nstate markers (n)",
       title = "Few shared DEGs are annotation markers") +
  theme_scm() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 1), legend.justification = c(0, 1),
        legend.key.size = unit(6, "pt"), legend.key.spacing.y = unit(1, "pt"),
        legend.text = element_text(size = FS$tick, margin = margin(l = 2)),
        axis.text.x = element_text(lineheight = 1.0, colour = COL$ink2, size = FS$tick),
        axis.title.y = element_text(lineheight = 1.05))

# ---- B: HPC high-effect sex-shared genes with sex difference > 1 log2 unit --------------------
gh <- g %>% filter(region == "HPC") %>% arrange(desc(sexbias_absF_minus_absM)) %>%
  mutate(y = rev(seq_len(n())), side = ifelse(sexbias_absF_minus_absM > 0, "women", "men"))
mk <- s2 %>% filter(region == "HPC") %>% distinct(gene, state_marker_lists)
mk <- mk[!duplicated(mk$gene), ]
gh <- gh %>% left_join(mk, by = "gene") %>%
  mutate(mlab = ifelse(is.na(state_marker_lists), "\u2013", paste(state_marker_lists, "list")))
chk("B genes shown (HPC, |dF-dM| > 1)", nrow(gh), 11)
chk("B larger in women / men", c(sum(gh$side == "women"), sum(gh$side == "men")), c(9, 2))
chk("B all high-effect sex-shared", sum(gh$effect_stratum == "high"), 11)
chk("B all |sex bias| > 1", sum(abs(gh$sexbias_absF_minus_absM) > 1), 11)
fem <- filter(gh, side == "women")
chk("B female-biased genes that are MG8/MG10 markers", sum(grepl("MG8|MG10", fem$state_marker_lists)), 7)
SEXC <- c(women = COL$women, men = COL$men)
XMK <- 3.95
pB <- ggplot(gh, aes(y = y)) +
  geom_vline(xintercept = 0, colour = COL$muted, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = sexbias_absF_minus_absM, yend = y, colour = side), linewidth = 0.6) +
  geom_point(aes(x = sexbias_absF_minus_absM, colour = side), size = 1.7) +
  geom_text(aes(x = XMK, label = mlab), hjust = 0, size = pt(6), family = FONT, colour = "#333333") +
  annotate("text", x = XMK, y = nrow(gh) + 0.55, label = "marker list", hjust = 0, vjust = 0, fontface = "bold",
           size = pt(6), family = FONT, colour = COL$ink) +
  annotate("text", x = -2.2, y = 0.45, label = "\u2190 larger in men", hjust = 0, size = pt(6), family = FONT,
           colour = COL$men) +
  annotate("text", x = 3.8, y = 0.45, label = "larger in women \u2192", hjust = 1, size = pt(6), family = FONT,
           colour = COL$women) +
  scale_colour_manual(values = SEXC, guide = "none") +
  scale_y_continuous(breaks = gh$y, labels = gh$gene, expand = c(0, 0)) +
  scale_x_continuous(breaks = -2:5, labels = lab_minus, expand = c(0, 0)) +
  coord_cartesian(xlim = c(-2.3, 5.2), ylim = c(0.2, nrow(gh) + 1.15), clip = "off") +
  labs(x = expression("|" * log[2] * FC[F] * "| \u2212 |" * log[2] * FC[M] * "| in HPC"), y = NULL,
       title = "7 of 9 female-biased HPC genes are MG8/MG10 markers") +
  theme_scm() +
  theme(axis.text.y = element_text(face = "italic", colour = COL$ink, size = FS$tick))

fig <- add_tags((pA | free(pB, type = "label", side = "b")) + plot_layout(widths = c(0.9, 1)))
save_fig(fig, "Supplementary_Figure_S2", width_mm = 183, height_mm = 90, outdir = "out")
if (dev.cur() > 1) invisible(dev.off()); invisible(gc())
writeLines(CHECKS, "out/Supplementary_Figure_S2_checks.txt")
cat(CHECKS, sep = "\n")
