# Supplementary Figure S5 - evidence behind the prior-knowledge networks (R / ggplot2 redraw)
# Inputs: inputs/supp_csv/S7.csv, S7_B_regulator_stats.csv, S8.csv ; style: theme_scm.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source("theme_scm.R")

s7 <- read.csv("inputs/supp_csv/S7.csv", stringsAsFactors = FALSE, check.names = FALSE)
rs <- read.csv("inputs/supp_csv/S7_B_regulator_stats.csv", stringsAsFactors = FALSE)
s8 <- read.csv("inputs/supp_csv/S8.csv", stringsAsFactors = FALSE)
has_member <- function(src, r) vapply(strsplit(src, "/", fixed = TRUE), function(x) r %in% x, logical(1))

# ------------------------------------------------------------------ A: sign composition
rs <- rs %>% arrange(desc(n_unique_core_targets)) %>% slice_head(n = 12)   # stable sort as pandas
nt <- s7 %>% filter(effect != "Technical", in_core_target == "yes")
ev <- bind_rows(lapply(rs$regulator_HGNC, function(r) {
  sub <- nt[has_member(nt$source_HGNC, r), ]
  data.frame(reg = r, effect = c("Activation", "Inhibition", "Unspecified"),
             n = c(sum(sub$effect == "Activation"), sum(sub$effect == "Inhibition"), sum(sub$effect == "Unspecified")))
}))
ev$effect <- factor(ev$effect, levels = c("Unspecified", "Inhibition", "Activation"))  # activation leftmost
ev$reg <- factor(ev$reg, levels = rev(rs$regulator_HGNC))
write.csv(pivot_wider(ev, names_from = effect, values_from = n), "out_S5A_counts.csv", row.names = FALSE)
lab_reg <- function(r) {
  core <- rs$regulator_in_core[match(r, rs$regulator_HGNC)] == "yes"
  tag <- ifelse(core, "'(core)'", "'(not core)'")
  is_gene <- !grepl("complex", r)
  nm <- ifelse(is_gene, sprintf("italic('%s')", r), sprintf("'%s'", sub("\\(", " (", r)))
  parse(text = paste0(nm, "~", tag))
}
pA <- ggplot(ev, aes(x = n, y = reg, fill = effect)) +
  geom_col(width = 0.7, position = position_stack()) +
  scale_fill_manual(values = c(Activation = COL$act, Inhibition = COL$inh, Unspecified = COL$uns),
                    breaks = c("Activation", "Inhibition", "Unspecified"),
                    labels = c("activation", "inhibition", "unspecified"), name = NULL) +
  scale_y_discrete(labels = lab_reg) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)), breaks = seq(0, 80, 10)) +
  labs(x = "non-technical edges to core genes\n(n; composite sources counted for each member)", y = NULL,
       title = "Most JunD/JunB edges to the core\nare of unspecified sign") +
  theme_scm() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.04), legend.justification = c(1, 0),
        axis.text.y = element_text(colour = COL$ink, size = FS$tick), axis.ticks.y = element_blank(),
        legend.key.size = unit(6, "pt"), legend.key.spacing.y = unit(2.5, "pt"))

# ------------------------------------------------------------------ B: references per edge
allnt <- s7 %>% filter(effect != "Technical")
jd <- allnt %>% filter(source_HGNC == "JUND")
jb <- allnt %>% filter(source_HGNC == "JUNB")
oth <- allnt[!(has_member(allnt$source_HGNC, "JUND") | has_member(allnt$source_HGNC, "JUNB")), ]
nref <- function(d) { x <- suppressWarnings(as.integer(d$n_references)); x[is.na(x)] <- 0L; x }
blab <- c("1", "2", "3\u20134", "5\u20139", "\u226510")
binit <- function(x) cut(x, c(0.5, 1.5, 2.5, 4.5, 9.5, 1000), labels = blab)
grp <- c(sprintf("JunD edges (n = %d)", nrow(jd)), sprintf("JunB edges (n = %d)", nrow(jb)),
         sprintf("all other regulators (n = %d)", nrow(oth)))
hb <- bind_rows(lapply(1:3, function(i) {
  d <- list(jd, jb, oth)[[i]]; b <- binit(nref(d))
  tb <- table(b); data.frame(group = grp[i], bin = names(tb), pct = as.numeric(tb) / sum(tb) * 100)
}))
hb$group <- factor(hb$group, levels = grp); hb$bin <- factor(hb$bin, levels = blab)
write.csv(hb, "out_S5B_hist.csv", row.names = FALSE)
pm_n <- sum(grepl("40455848", ifelse(is.na(jd$references_PMID), "", jd$references_PMID), fixed = TRUE))
pB <- ggplot(hb, aes(bin, pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.81), width = 0.78) +
  scale_fill_manual(values = setNames(c("#4d4d4d", "#969696", "#cccccc"), grp), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA), breaks = seq(0, 80, 20)) +
  labs(x = "supporting references per edge (n)", y = "edges (%)",
       title = sprintf("JunD edges rest on one study:\n%d of %d single-source JunD edges cite PMID 40455848", pm_n, nrow(jd))) +
  theme_scm() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1), plot.margin = margin(3, 10, 3, 3),
        legend.key.size = unit(6, "pt"), legend.key.spacing.y = unit(2.5, "pt"))

# ------------------------------------------------------------------ C: decoupleR / CollecTRI dot matrix
tfs <- c("SP1", "STAT3", "JUN", "NFKB1", "RELA", "HIF1A", "NFKB", "CREB1", "AP1")
tfs <- c(tfs[tfs %in% s8$TF_source], sort(setdiff(unique(s8$TF_source), tfs)))
tg_m <- unique(s8$target[s8$sex == "men"]); tg_w <- unique(s8$target[s8$sex == "women"])
targets <- c(tg_m, setdiff(tg_w, tg_m))
d8 <- s8 %>% mutate(x = match(target, targets), y = match(TF_source, rev(tfs)),
                    col = ifelse(sign == "activating", COL$act, COL$inh),
                    fill = ifelse(target_in_final_HPC_list == "yes", col, "white"),
                    shp = ifelse(sex == "men", 21, 22))
nx_ <- length(targets); ny_ <- length(tfs)
# in-panel key (placed in the empty right-hand strip, as in the reference)
kx <- nx_ + 1.4; ky <- c(6.4, 5.7, 5.0, 4.3)
key <- data.frame(x = kx, y = ky, shp = c(21, 21, 22, 21),
                  col = c(COL$act, COL$inh, COL$act, COL$ink2), fill = c(COL$act, COL$inh, COL$act, "white"),
                  lab = c("men, activating", "men, repressing", "women, activating", "target not in final HPC list"))
nm <- sum(s8$sex == "men"); nw <- sum(s8$sex == "women")
pC <- ggplot(d8, aes(x, y)) +
  geom_hline(yintercept = 1:ny_, colour = "#eeeeee", linewidth = 0.25) +
  geom_vline(xintercept = 1:nx_, colour = "#eeeeee", linewidth = 0.25) +
  geom_point(aes(shape = I(shp), colour = I(col), fill = I(fill)), size = 2.3, stroke = 0.7) +
  geom_point(data = key, aes(shape = I(shp), colour = I(col), fill = I(fill)), size = 2.3, stroke = 0.7) +
  geom_text(data = key, aes(x = x + 0.35, label = lab), hjust = 0, size = pt(FS$small), family = FONT, colour = COL$ink) +
  scale_x_continuous(breaks = seq_len(nx_), labels = targets, limits = c(0.4, nx_ + 5.6), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_len(ny_), labels = parse(text = ifelse(rev(tfs) %in% c("NFKB", "AP1"), sprintf("'%s'", rev(tfs)), sprintf("italic('%s')", rev(tfs)))), limits = c(0.4, ny_ + 0.6), expand = c(0, 0)) +
  labs(x = "sex-specific target gene", y = "transcription factor",
       title = sprintf("CollecTRI edges from the 10 top-ranked decoupleR TFs to HPC sex-specific targets\n(men %d edges, women %d; run on an earlier, unfiltered HPC list)", nm, nw)) +
  theme_scm() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, face = "italic", colour = COL$ink),
        axis.text.y = element_text(colour = COL$ink))

# ------------------------------------------------------------------ layout + export
fig <- (pA | plot_spacer() | pB) + plot_layout(widths = c(1, 0.06, 1))
fig <- fig / free(pC, type = "label", side = "l") + plot_layout(heights = c(1, 1.05))
fig <- add_tags(fig)
dir.create("out", showWarnings = FALSE)
save_fig(fig, "Supplementary_Figure_S5", width_mm = 183, height_mm = 130)
cat("A counts:\n"); print(as.data.frame(pivot_wider(ev, names_from = effect, values_from = n)))
cat(sprintf("B: JunD n=%d, JunB n=%d, others n=%d, JunD citing 40455848=%d, JunD single-ref=%d\n",
            nrow(jd), nrow(jb), nrow(oth), pm_n, sum(nref(jd) == 1)))
cat(sprintf("C: men %d edges, women %d edges, TFs %d, targets %d, open symbols %d\n", nm, nw, ny_, nx_,
            sum(s8$target_in_final_HPC_list != "yes")))
invisible(gc())
