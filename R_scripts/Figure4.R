# Figure 4 - prior-knowledge regulatory networks (R / ggplot2 redraw of fig4_submission.py)
#   A  simplified MetaCore map "substance P-mediated inflammation" with log2FC mini-heatmaps
#   B  curated MetaCore network of the 109 core genes (+ JunB, JunD)
#   C  TF-TF subnetwork among 14 regulators
#   D  MetaCore expansions: male (vector redraw) and female (MetaCore raster + call-outs)
#   E  genetic-algorithm Boolean networks (male, female)
# All coordinates are in inches of the final figure (1 data unit = 1 inch, coord_fixed), so
# node radii / arrow shrink given in points are converted exactly (pt / 72).
suppressPackageStartupMessages({ library(dplyr); library(png); library(grid) })
source("theme_scm.R")

SD <- "inputs/source_data/"
INK <- COL$ink; INK2 <- COL$ink2; MUTED <- COL$muted; HL <- COL$highlight
ECOL <- c(Activation = COL$act, Inhibition = COL$inh, Unspecified = COL$uns, Both = COL$both)
NODE_TF <- "#8f8d88"; NODE_T <- "#e6e5e1"
FILL <- c("1" = "#f2a3a2", "-1" = "#86b6ef")
FEEDBACK <- c("NR4A1", "NR4A3", "BCL6", "BCL3", "SOCS3", "DUSP1", "ZFP36L1")
NONCORE <- c("JUNB", "JUND", "KDM3A")
R14 <- c("JUN", "JUNB", "JUND", "FOSL2", "EGR1", "EGR3", "NR4A1", "NR4A3", "BCL6", "BCL3", "SMARCA4", "MKL1",
         "ZFP36L1", "KDM3A")
EORDER <- c(Unspecified = 0, Activation = 1, Inhibition = 2, Both = 2)
FS6 <- pt(6)                                    # 6-pt text (geom_text size units)
lwp <- function(p) p * 96 / 72 / .pt            # line width in pt -> ggplot linewidth
P2I <- 1 / 72                                   # pt -> inch
DASH <- "22"

# ============================================================== geometry helpers (inches)
circ <- function(x, y, r, n = 48) { t <- seq(0, 2 * pi, length.out = n + 1)[-1]; data.frame(x = x + r * cos(t), y = y + r * sin(t)) }
rrect <- function(x0, y0, w, h, r = 0.04, n = 8) {
  th <- seq(0, pi / 2, length.out = n)
  cx <- c(x0 + w - r, x0 + r, x0 + r, x0 + w - r); cy <- c(y0 + h - r, y0 + h - r, y0 + r, y0 + r)
  do.call(rbind, lapply(1:4, function(k) data.frame(x = cx[k] + r * cos((k - 1) * pi / 2 + th),
                                                     y = cy[k] + r * sin((k - 1) * pi / 2 + th))))
}
# quadratic Bezier equivalent to matplotlib connectionstyle arc3(rad), trimmed by sa / sb (inches) at the ends
bez <- function(x1, y1, x2, y2, rad, sa, sb, n = 120) {
  dx <- x2 - x1; dy <- y2 - y1; cx <- (x1 + x2) / 2 + rad * dy; cy <- (y1 + y2) / 2 - rad * dx
  t <- seq(0, 1, length.out = n)
  x <- (1 - t)^2 * x1 + 2 * (1 - t) * t * cx + t^2 * x2; y <- (1 - t)^2 * y1 + 2 * (1 - t) * t * cy + t^2 * y2
  keep <- sqrt((x - x1)^2 + (y - y1)^2) >= sa & sqrt((x - x2)^2 + (y - y2)^2) >= sb
  if (sum(keep) < 2) keep[c(n %/% 2, n %/% 2 + 1)] <- TRUE
  data.frame(x = x[keep], y = y[keep])
}
# edges: data.frame(u, v, effect, hl, rad, lw (pt), alpha); pos: data.frame(name, x, y); r: named radii (pt)
# returns list(path, bars, hl) ready for geom_path / geom_segment
build_edges <- function(ed, pos, r, sa_add = 0.5, sb_add = 1.0, bar_pt = 2.4) {
  P <- list(); B <- list(); H <- list()
  for (i in seq_len(nrow(ed))) {
    e <- ed[i, ]; a <- pos[pos$name == e$u, ]; b <- pos[pos$name == e$v, ]
    p <- bez(a$x, a$y, b$x, b$y, e$rad, (r[[e$u]] + sa_add) * P2I, (r[[e$v]] + sb_add) * P2I)
    p$id <- i; p$col <- ECOL[[e$effect]]; p$lw <- e$lw; p$alpha <- e$alpha; p$z <- EORDER[[e$effect]]
    p$head <- e$effect %in% c("Activation", "Both")
    P[[i]] <- p
    if (e$effect == "Inhibition") {
      m <- nrow(p); d <- c(p$x[m] - p$x[m - 1], p$y[m] - p$y[m - 1]); d <- d / sqrt(sum(d^2)); nrm <- c(-d[2], d[1]) * bar_pt * P2I
      B[[i]] <- data.frame(x = p$x[m] - nrm[1], y = p$y[m] - nrm[2], xend = p$x[m] + nrm[1], yend = p$y[m] + nrm[2],
                           col = p$col[1], lw = e$lw, alpha = e$alpha)
    }
    if (isTRUE(e$hl)) H[[i]] <- transform(p, lw = e$lw + 3.2)
  }
  list(path = do.call(rbind, P), bars = do.call(rbind, B), hl = do.call(rbind, H))
}
edge_layers <- function(E, head_len = 1.2) {
  pa <- E$path
  L <- list()
  if (!is.null(E$hl)) L <- c(L, list(geom_path(data = E$hl, aes(x, y, group = id), colour = HL, alpha = 0.55,
                                               linewidth = lwp(E$hl$lw[1]), lineend = "round")))
  for (z in sort(unique(pa$z))) {
    s <- pa[pa$z == z & !pa$head, ]; h <- pa[pa$z == z & pa$head, ]
    if (nrow(s)) L <- c(L, list(geom_path(data = s, aes(x, y, group = id, colour = col, alpha = alpha, linewidth = lwp(lw)))))
    if (nrow(h)) L <- c(L, list(geom_path(data = h, aes(x, y, group = id, colour = col, alpha = alpha, linewidth = lwp(lw)),
                                          arrow = arrow(angle = 24, length = unit(head_len, "mm"), type = "closed"),
                                          linejoin = "mitre")))
    if (!is.null(E$bars) && z == 2) L <- c(L, list(geom_segment(data = E$bars, aes(x, y, xend = xend, yend = yend, colour = col,
                                                                                    alpha = alpha, linewidth = lwp(lw)))))
  }
  L
}
self_loop_df <- function(x, y, ang_deg, size, id) {
  a <- ang_deg * pi / 180; d <- circ(x + cos(a) * size, y + sin(a) * size, size * 0.55, n = 60)
  rbind(d, d[1, ])  |> transform(id = id)
}
node_df <- function(nodes) {          # nodes: name, x, y, r (inch), fill, col, lw (pt), lt
  do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) { n <- nodes[i, ]
    transform(circ(n$x, n$y, n$r), id = i, fill = n$fill, col = n$col, lw = n$lw, lt = n$lt) }))
}
node_layer <- function(nd) geom_polygon(data = nd, aes(x, y, group = id, fill = fill, colour = col,
                                                         linewidth = lwp(lw), linetype = lt))
canvas <- function(p, xlim, ylim) p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
  scale_fill_identity() + scale_colour_identity() + scale_alpha_identity() + scale_linewidth_identity() +
  scale_linetype_identity() + theme_void_scm()
# split "SYMBOL\n(alias)" into rows: line 1 italic symbol, line 2 upright alias
bg_text <- function(data, mapping, pad = 0.4, alpha_bg = 0.75, r = 0.5, ...) list(
  geom_label(data = data, mapping = mapping, text.colour = NA, fill = alpha("white", alpha_bg), label.size = 0, border.colour = NA,
             label.padding = unit(pad, "pt"), label.r = unit(r, "pt"), size = FS6, family = FONT, ...),
  geom_text(data = data, mapping = mapping, colour = INK, size = FS6, family = FONT, ...))
two_line <- function(x, y, l1, l2, vj, lh = 6.6 * P2I, face1 = "italic", face2 = "plain") {
  if (vj == 0)   { y1 <- y + lh; y2 <- y } else if (vj == 1) { y1 <- y; y2 <- y - lh } else { y1 <- y + lh / 2; y2 <- y - lh / 2 }
  data.frame(x = x, y = c(y1, y2), lab = c(l1, l2), face = c(face1, face2), vj = vj)
}

# ============================================================== A: MetaCore map redraw
mp <- read.csv(file.path(SD, "Fig4A_map_nodes_log2FC.csv"), stringsAsFactors = FALSE)
lfc_pal <- scales::gradient_n_pal(DIV_COLS)                     # = scale_fill_lfc colours, limits -4..4
lfc_col <- function(v) lfc_pal(scales::rescale(pmin(pmax(v, -4), 4), from = c(-4, 4)))
lfc <- function(g) {
  out <- matrix(NA_real_, 2, 2, dimnames = list(c("HPC", "PFC"), c("M", "F")))
  for (reg in c("HPC", "PFC")) { r <- mp[mp$gene == g & mp$region == reg & mp$in_filtered_list == "True", ]
    if (nrow(r)) out[reg, ] <- c(r$log2FC_men[1], r$log2FC_women[1]) }
  out
}
A_box <- list(); A_txt <- list(); A_heat <- list(); A_seg <- list(); A_arr <- list()
addbox <- function(d, fill, col, lw, lt) { k <- length(A_box) + 1; A_box[[k]] <<- transform(d, id = k, fill = fill, col = col, lw = lw, lt = lt) }
addtxt <- function(x, y, lab, col = INK, hj = 0.5, vj = 0.5, size = FS6, parse = TRUE) A_txt[[length(A_txt) + 1]] <<-
  data.frame(x = x, y = y, lab = lab, col = col, hj = hj, vj = vj, size = size, parse = parse)
seg <- function(x, y, xend, yend, col, lw) A_seg[[length(A_seg) + 1]] <<- data.frame(x = x, y = y, xend = xend, yend = yend, col = col, lw = lw)
arr <- function(x, y, xend, yend, col, lw) A_arr[[length(A_arr) + 1]] <<- data.frame(x = x, y = y, xend = xend, yend = yend, col = col, lw = lw)
measured_genes <- character(0)
heat <- function(x, y, g, cw = 0.215, ch = 0.115) {
  v <- lfc(g)
  for (i in 1:2) for (j in 1:2) {
    val <- v[i, j]; xx <- x + (j - 1) * cw; yy <- y - i * ch
    A_heat[[length(A_heat) + 1]] <<- data.frame(xmin = xx, xmax = xx + cw, ymin = yy, ymax = yy + ch,
      fill = if (is.na(val)) "#f4f3f0" else lfc_col(val),
      lab = if (is.na(val)) "\u2013" else fmt_minus(sprintf("%.1f", val)),
      tcol = if (is.na(val)) MUTED else if (abs(val) > 2.2) "white" else INK)
  }
}
node <- function(x, y, lab, gene = NULL, w = 0.56, h = 0.44, lab2 = NULL) {
  measured <- !is.null(gene) && any(!is.na(lfc(gene)))
  addbox(rrect(x - w / 2, y - h / 2, w, h, 0.04), "white", if (measured) INK else MUTED, if (measured) 0.8 else 0.5,
         if (measured) "solid" else DASH)
  if (measured) {
    measured_genes <<- c(measured_genes, gene)
    addtxt(x, y + h / 2 - 0.035, lab, INK, vj = 1)
    heat(x - 0.215, y + h / 2 - 0.165, gene)
  } else if (is.null(lab2)) addtxt(x, y, lab, INK2) else { addtxt(x, y + 0.054, lab, INK2); addtxt(x, y - 0.054, lab2, INK2) }
}
q <- function(s) sprintf("'%s'", s)                      # quoted plotmath string
it <- function(s) sprintf("italic('%s')", s)
par_it <- function(s) sprintf("'('*italic('%s')*')'", s) # (SYMBOL) with italic symbol
# upstream signalling (map context)
addbox(rrect(0.18, 2.72, 2.86, 0.38, 0.05), "#f7f6f3", MUTED, 0.5, "solid")
addtxt(1.61, 2.965, paste0(q("Substance P \u2192 NK1 receptor \u2192 G\u03b1q / PLC-\u03b2 \u2192 PKC, Ca"), "^'2+'"), INK2)
addtxt(1.61, 2.855, q("\u2192 ERK1/2, JNK, p38 MAPK and IKK (map context)"), INK2)
xj <- 0.95; xn <- 2.27; yt <- 2.25
for (xx in c(xj, xn)) arr(xx, 2.72, xx, yt + 0.23, INK2, 0.6)
node(xj, yt, "bold('c-Jun (')*bolditalic('JUN')*bold(')')", "JUN", w = 0.70)
node(xn, yt, q("NF-\u03baB"), NULL, w = 0.70, lab2 = q("p50/p65"))
ybus <- 1.86
for (xx in c(xj, xn)) seg(xx, yt - 0.22, xx, ybus, COL$act, 0.8)
seg(0.32, ybus, 3.17, ybus, COL$act, 0.8)
A_lab_tr <- data.frame(x = 1.61, y = ybus + 0.03, lab = "transcriptional regulation")
row1 <- list(list("bolditalic('CCL3L1')", "CCL3L1"), list("bold('MIP-1\u03b2/')*bolditalic('CCL4')", "CCL4"),
             list(it("CCL13"), NULL), list(q("COX-2"), NULL, par_it("PTGS2")), list(q("TNF-\u03b1"), NULL))
row2 <- list(list("bolditalic('CCL2')", "CCL2"), list(q("IL-8"), NULL, par_it("CXCL8")), list(q("IL-6"), NULL),
             list(q("GRO-2"), NULL, par_it("CXCL2")), list(q("GRO-3"), NULL, par_it("CXCL3")))
xs1 <- 0.32 + (0:4) * 0.645; y1 <- 1.40; y2 <- 0.80
for (i in 1:5) { arr(xs1[i], ybus, xs1[i], y1 + 0.225, COL$act, 0.6)
  r <- row1[[i]]; node(xs1[i], y1, r[[1]], r[[2]], lab2 = if (length(r) > 2) r[[3]] else NULL) }
ybus2 <- y2 + 0.33
seg(3.17, ybus, 3.17, ybus2, COL$act, 0.8); seg(xs1[1], ybus2, 3.17, ybus2, COL$act, 0.8)
for (i in 1:5) { arr(xs1[i], ybus2, xs1[i], y2 + 0.225, COL$act, 0.6)
  r <- row2[[i]]; node(xs1[i], y2, r[[1]], r[[2]], lab2 = if (length(r) > 2) r[[3]] else NULL) }
seg(0.32, 0.50, 2.90, 0.50, MUTED, 0.6)
for (xx in xs1) seg(xx, y2 - 0.22, xx, 0.50, MUTED, 0.6)
arr(1.61, 0.50, 1.61, 0.36, MUTED, 0.6)
addbox(rrect(0.95, 0.14, 1.32, 0.22, 0.04), "#f7f6f3", MUTED, 0.5, "solid")
addtxt(1.61, 0.25, q("Inflammatory response"), INK2)
# key for the mini-heatmaps
kx <- 2.45; ky <- 0.30; kw <- 0.13; kh <- 0.09
A_key <- expand.grid(i = 1:2, j = 1:2) |> transform(xmin = kx + (j - 1) * kw, ymin = ky + 0.08 - i * kh) |>
  transform(xmax = xmin + kw, ymax = ymin + kh)
addtxt(kx + kw / 2, ky + 0.09, q("M"), INK2, vj = 0); addtxt(kx + 1.5 * kw, ky + 0.09, q("F"), INK2, vj = 0)
addtxt(kx + 2 * kw + 0.04, ky + 0.08 - 0.5 * kh, q("HPC"), INK2, hj = 0); addtxt(kx + 2 * kw + 0.04, ky + 0.08 - 1.5 * kh, q("PFC"), INK2, hj = 0)

A_box <- do.call(rbind, A_box); A_txt <- do.call(rbind, A_txt); A_heat <- do.call(rbind, A_heat)
A_seg <- do.call(rbind, A_seg); A_arr <- do.call(rbind, A_arr)
pA <- ggplot() +
  geom_segment(data = A_seg, aes(x, y, xend = xend, yend = yend, colour = col, linewidth = lwp(lw))) +
  geom_segment(data = A_arr, aes(x, y, xend = xend, yend = yend, colour = col, linewidth = lwp(lw)),
               arrow = arrow(angle = 24, length = unit(1.4, "mm"), type = "closed"), linejoin = "mitre") +
  geom_polygon(data = A_box, aes(x, y, group = id, fill = fill, colour = col, linewidth = lwp(lw), linetype = lt)) +
  geom_label(data = A_lab_tr, aes(x, y, label = lab), colour = COL$act, fill = "white", size = FS6, family = FONT,
             vjust = 0, label.size = 0, border.colour = NA, label.padding = unit(0.6, "pt")) +
  geom_rect(data = A_heat, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "white", linewidth = lwp(0.5)) +
  geom_text(data = A_heat, aes((xmin + xmax) / 2, (ymin + ymax) / 2, label = lab, colour = tcol), size = FS6, family = FONT) +
  geom_rect(data = A_key, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#f4f3f0", colour = INK2, linewidth = lwp(0.4)) +
  geom_text(data = A_txt, aes(x, y, label = lab, colour = col, hjust = hj, vjust = vj), size = FS6, family = FONT, parse = TRUE) +
  labs(title = "MetaCore map: substance P-mediated inflammation",
       subtitle = expression("simplified redraw; DEG nodes show measured " * log[2] * "FC"))
pA <- canvas(pA, c(0, 3.24), c(0.08, 3.14))

# ============================================================== B: curated core network
EB <- read.csv(file.path(SD, "Fig4B_core_network_edges.csv"), stringsAsFactors = FALSE)
NB <- read.csv(file.path(SD, "Fig4B_core_network_nodes.csv"), stringsAsFactors = FALSE)
selfB <- EB[EB$src == EB$tgt, ]; EBd <- EB[EB$src != EB$tgt, ]
nodesB <- sort(NB$node, method = "radix")
outdeg <- setNames(as.integer(table(factor(EBd$src, levels = nodesB))), nodesB)
stopifnot(all(outdeg[NB$node] == NB$out_degree))
predB <- split(EBd$src, factor(EBd$tgt, levels = nodesB))
inner <- c("JUND", "JUNB", "FOSL2", "JUN", "SMARCA4", "BCL6", "EGR1", "MKL1")
TFset <- c(R14, "KDM2A", "ZHX2", "POU2F2")
regs <- nodesB[outdeg > 0]
middle <- regs[!regs %in% inner & regs %in% TFset]
outer_src <- regs[!regs %in% inner & !regs %in% TFset]
targetsB <- nodesB[outdeg == 0]
dominant <- function(n) { p <- predB[[n]]; k <- which(inner %in% p); if (length(k)) k[1] - 1 else length(inner) }
middle <- middle[order(sapply(middle, dominant), middle, method = "radix")]
outer <- c(targetsB, outer_src)
k1 <- sapply(outer, function(n) { p <- predB[[n]]; m <- match(p[p %in% inner], inner) - 1; min(c(m, length(inner))) })
k2 <- sapply(outer, function(n) { p <- sort(predB[[n]], method = "radix"); if (length(p)) p[1] else "" })
outer <- outer[order(k1, k2, outer, method = "radix")]
sB <- min(3.45 / 2.84, 3.12 / 2.96)                 # inch per layout unit in the reference axes
posB <- rbind(
  data.frame(name = inner, a = 2 * pi * (seq_along(inner) - 1) / length(inner) + 0.2, rr = 0.30),
  data.frame(name = middle, a = 2 * pi * (seq_along(middle) - 1) / length(middle) + 0.15, rr = 0.62),
  data.frame(name = outer, a = 2 * pi * (seq_along(outer) - 1) / length(outer), rr = 1.0))
posB <- transform(posB, x = rr * cos(a) * sB, y = rr * sin(a) * sB)
stopifnot(nrow(posB) == length(nodesB), setequal(posB$name, nodesB))
rB <- setNames(2.2 + 1.35 * sqrt(outdeg[posB$name]), posB$name)          # radius in pt
edB <- EBd |> transmute(u = src, v = tgt, effect, hl = FALSE, rad = 0.06,
                        lw = ifelse(effect == "Unspecified", 0.35, 0.5), alpha = ifelse(effect == "Unspecified", 0.55, 0.85))
EBg <- build_edges(edB, posB, rB, bar_pt = 2.2)
ndB <- data.frame(name = nodesB) |> left_join(posB, by = "name") |>
  mutate(r = rB[name] * P2I, fill = ifelse(name %in% TFset, NODE_TF, NODE_T), col = INK2, lw = 0.4, lt = "solid")
ndB <- ndB |> mutate(col = ifelse(name %in% FEEDBACK, HL, col), lw = ifelse(name %in% FEEDBACK, 1.4, lw)) |>
  mutate(fill = ifelse(name %in% NONCORE, "white", fill), col = ifelse(name %in% NONCORE, INK, col),
         lw = ifelse(name %in% NONCORE, 0.8, lw), lt = ifelse(name %in% NONCORE, DASH, lt))
LAB <- intersect(nodesB, unique(c(inner, middle, FEEDBACK, "CCL2", "CCL3L1", "IL1B", "SOCS3", "BCL2", "CXCR4", "CD83", "SRGN", "ACSL1")))
ALIAS <- c(SMARCA4 = "(BRG1)", FOSL2 = "(Fra-2)", KDM3A = "(JMJD1A)")
labIn <- do.call(rbind, lapply(inner, function(n) { p <- posB[posB$name == n, ]
  if (n %in% names(ALIAS)) two_line(p$x, p$y, n, ALIAS[[n]], 0.5, lh = 6.0 * P2I, "bold.italic", "bold")
  else data.frame(x = p$x, y = p$y, lab = n, face = "bold.italic", vj = 0.5) }))
# small label-only nudges (inches) so neighbouring inner-hub labels do not overprint (node positions unchanged)
NUDGE <- list(FOSL2 = c(-0.07, 0.01), JUNB = c(0.07, 0.03), BCL6 = c(-0.05, 0), EGR1 = c(0.05, -0.01))
for (n in names(NUDGE)) { k <- labIn$lab == n | (labIn$lab %in% ALIAS[n] & !is.na(ALIAS[n])); labIn$x[k] <- labIn$x[k] + NUDGE[[n]][1]; labIn$y[k] <- labIn$y[k] + NUDGE[[n]][2] }
labOut <- posB[posB$name %in% setdiff(LAB, inner), ] |>
  mutate(ux = cos(a), uy = sin(a), lx = x + ux * 0.075 * sB, ly = y + uy * 0.075 * sB,
         hj = ifelse(ux > 0.25, 0, ifelse(ux < -0.25, 1, 0.5)), vj = ifelse(uy > 0.25, 0, ifelse(uy < -0.25, 1, 0.5)))
# legend (top-left, as in the reference) and footnote
lgx <- -1.534 * sB; lgy <- 1.569 * sB - 0.055; dy <- 0.118
lgB_l <- data.frame(x = lgx, xend = lgx + 0.13, y = lgy - (0:3) * dy, col = unname(ECOL[c("Activation", "Inhibition", "Unspecified", "Both")]))
lgB_n <- data.frame(name = c("tf", "fb", "nc"), x = lgx + 0.065, y = lgy - (4:6) * dy, r = 2.6 * P2I,
                    fill = c(NODE_TF, "white", "white"), col = c(INK2, HL, INK), lw = c(0.4, 1.4, 0.8), lt = c("solid", "solid", DASH))
lgB_t <- data.frame(x = lgx + 0.19, y = lgy - (0:6) * dy,
                    lab = c("activation", "inhibition", "unspecified", "both signs", "TF / regulator", "feedback regulator", "not a core gene"))
nB_nodes <- length(nodesB); nB_pairs <- nrow(EBd)
footB <- sprintf("%d nodes, %d directed pairs\n(%d self-loop not drawn); node size proportional to out-degree", nB_nodes, nB_pairs, nrow(selfB))
pB <- ggplot() + edge_layers(EBg, head_len = 1.1) + node_layer(node_df(ndB)) +
  bg_text(labOut, aes(lx, ly, label = name, hjust = hj, vjust = vj), pad = 0.4, alpha_bg = 0.7, fontface = "italic") +
  bg_text(labIn, aes(x, y, label = lab, vjust = vj, fontface = face), pad = 0.6, alpha_bg = 0.75, r = 0.8) +
  geom_segment(data = lgB_l, aes(x, y, xend = xend, yend = y, colour = col), linewidth = lwp(1)) +
  node_layer(node_df(lgB_n)) +
  geom_text(data = lgB_t, aes(x, y, label = lab), hjust = 0, size = FS6, family = FONT, colour = INK) +
  annotate("text", x = 1.534 * sB, y = -1.273 * sB, label = footB, hjust = 1, vjust = 0, size = FS6, family = FONT,
           colour = INK2, lineheight = 0.95) +
  labs(title = "Curated MetaCore network of the 109 core genes")
pB <- canvas(pB, c(-1.66, 1.66), c(-1.52, 1.66))

# ============================================================== C: TF-TF subnetwork
EC <- read.csv(file.path(SD, "Fig4C_TF_TF_edges.csv"), stringsAsFactors = FALSE)
s7 <- read.csv("inputs/supp_csv/S7.csv", stringsAsFactors = FALSE)
anyR14 <- function(v) vapply(strsplit(v, "/", fixed = TRUE), function(x) any(x %in% R14), logical(1))
s7nt <- s7[s7$effect != "Technical", ]
nC_records <- sum(anyR14(s7nt$source_HGNC) & anyR14(s7nt$target_HGNC))
famorder <- c("JUN", "JUNB", "JUND", "FOSL2", "BCL3", "BCL6", "ZFP36L1", "NR4A1", "NR4A3", "EGR1", "EGR3", "MKL1", "SMARCA4", "KDM3A")
sC <- 2.15 / 2.76
posC <- data.frame(name = famorder, a = pi / 2 + 2 * pi * (0:13) / 14) |> transform(x = cos(a) * sC, y = sin(a) * sC)
odC <- setNames(as.integer(table(factor(EC$src, levels = famorder))), famorder)   # out-degree incl. self-loops (networkx)
rC <- setNames(3.2 + sqrt(odC), famorder)
ECd <- EC[EC$src != EC$tgt, ]; ECs <- EC[EC$src == EC$tgt, ]
edC <- ECd |> transmute(u = src, v = tgt, effect, hl = FALSE, rad = 0.12,
                        lw = ifelse(effect == "Unspecified", 0.45, 0.6), alpha = ifelse(effect == "Unspecified", 0.7, 0.9))
ECg <- build_edges(edC, posC, rC, bar_pt = 2.6)
slC <- do.call(rbind, lapply(seq_len(nrow(ECs)), function(i) { p <- posC[posC$name == ECs$src[i], ]
  transform(self_loop_df(p$x, p$y, p$a * 180 / pi, 0.11 * sC, i), col = ECOL[[ECs$effect[i]]]) }))
ndC <- posC |> mutate(r = rC[name] * P2I, fill = NODE_TF, col = INK2, lw = 0.4, lt = "solid") |>
  mutate(col = ifelse(name %in% FEEDBACK, HL, col), lw = ifelse(name %in% FEEDBACK, 1.4, lw),
         fill = ifelse(name %in% NONCORE, "white", fill), col = ifelse(name %in% NONCORE, INK, col),
         lw = ifelse(name %in% NONCORE, 0.8, lw), lt = ifelse(name %in% NONCORE, DASH, lt))
labC <- posC |> mutate(lx = x * 1.2, ly = y * 1.2, hj = ifelse(cos(a) > 0.3, 0, ifelse(cos(a) < -0.3, 1, 0.5)),
                       vj = ifelse(sin(a) > 0.3, 0, ifelse(sin(a) < -0.3, 1, 0.5)),
                       lab = ifelse(name %in% c("FOSL2", "KDM3A"), sprintf("italic('%s')~'%s'", name, ALIAS[name]), it(name)))
labC_s <- labC[labC$name != "SMARCA4", ]
pS <- labC[labC$name == "SMARCA4", ]
labC2 <- two_line(pS$lx, pS$ly, "SMARCA4", "(BRG1)", pS$vj, lh = 6.4 * P2I) |> transform(hj = pS$hj)
pC <- ggplot() + edge_layers(ECg, head_len = 1.2) +
  geom_path(data = slC, aes(x, y, group = id, colour = col), linewidth = lwp(0.6)) +
  node_layer(node_df(ndC)) +
  geom_text(data = labC_s, aes(lx, ly, label = lab, hjust = hj, vjust = vj), size = FS6, family = FONT, colour = INK, parse = TRUE) +
  geom_text(data = labC2, aes(x, y, label = lab, hjust = hj, vjust = vj, fontface = face), size = FS6, family = FONT, colour = INK) +
  labs(title = "TF\u2013TF edges among 14 regulators",
       subtitle = sprintf("%d interaction records, %d directed pairs", nC_records, nrow(EC)))
pC <- canvas(pC, c(-1.40, 1.26), c(-1.02, 1.00))

# ============================================================== D: MetaCore expansions
DM <- read.csv(file.path(SD, "Fig4D_male_expansion_edges_transcribed.csv"), stringsAsFactors = FALSE)
PMpx <- data.frame(name = c("NODAL", "FGF4", "Esrrb", "CDX2", "MED12L", "SDF-1 (CXCL12)", "DBL (MCF2)", "SMAD3", "ALK-4",
                            "Oct-3/4", "KCNH5", "SMAD2", "PRDM11", "DHRS9"),
                   px = c(385, 805, 1040, 1175, 80, 760, 1330, 510, 225, 955, 1340, 590, 880, 1140),
                   py = c(70, 105, 85, 215, 305, 300, 455, 485, 555, 540, 800, 785, 965, 905))
posD <- transform(PMpx, x = px / 1500 * 1.85, y = (1060 - py) / 1060 * 1.32)
stopifnot(all(c(DM$source, DM$target) %in% posD$name))
MALE_DIR <- c("MED12L" = -1, "SDF-1 (CXCL12)" = -1, "Esrrb" = -1, "DBL (MCF2)" = 1, "KCNH5" = 1, "PRDM11" = -1, "DHRS9" = -1)
rD <- setNames(rep(5, nrow(posD)), posD$name)
bidir <- mapply(function(u, v) any(DM$source == v & DM$target == u), DM$source, DM$target)
edD <- DM |> transmute(u = source, v = target, effect, hl = FALSE, rad = ifelse(bidir, 0.18, 0.04),
                       lw = ifelse(effect == "Unspecified", 0.45, 0.6), alpha = ifelse(effect == "Unspecified", 0.75, 0.95))
EDg <- build_edges(edD, posD, rD, bar_pt = 2.6)
ndD <- posD |> mutate(r = 5 * P2I, d = MALE_DIR[name],
                      fill = ifelse(is.na(d), "white", FILL[as.character(d)]), col = ifelse(is.na(d), INK, INK2),
                      lw = ifelse(is.na(d), 0.7, 0.4), lt = ifelse(is.na(d), DASH, "solid"))
HGNC_D <- c("NODAL", "FGF4", "CDX2", "MED12L", "SMAD3", "SMAD2", "KCNH5", "PRDM11", "DHRS9")
LOFF <- list("SDF-1 (CXCL12)" = c(-0.075, 0, 1, 0.5), "Esrrb" = c(0, 0.075, 0.5, 0), "KCNH5" = c(0, -0.075, 0.5, 1),
             "CDX2" = c(0.075, 0, 0, 0.5))
labD <- do.call(rbind, lapply(seq_len(nrow(posD)), function(i) {
  n <- posD$name[i]; o <- if (!is.null(LOFF[[n]])) LOFF[[n]] else c(0, -0.075, 0.5, 1)
  x <- posD$x[i] + o[1]; y <- posD$y[i] + o[2]
  if (grepl(" (", n, fixed = TRUE)) { pp <- strsplit(n, " ", fixed = TRUE)[[1]]
    d <- two_line(x, y, pp[1], pp[2], o[4], lh = 6.4 * P2I, "plain", "italic"); d$hj <- o[3]; d }
  else data.frame(x = x, y = y, lab = n, face = ifelse(n %in% HGNC_D, "italic", "plain"), vj = o[4], hj = o[3])
}))
# female expansion: MetaCore export, cropped as in the reference script and composited on white
im <- png::readPNG("inputs/images/female_network.png")
crop <- im[1:2140, 1360:4453, , drop = FALSE]; rm(im); invisible(gc())
a <- crop[, , 4]; rgbc <- crop[, , 1:3] * as.vector(a) + (1 - as.vector(a)); rm(crop, a); invisible(gc())
h2 <- floor(dim(rgbc)[1] / 2) * 2; w2 <- floor(dim(rgbc)[2] / 2) * 2         # 2x2 block average (~650 dpi at print size)
ds <- array(0, c(h2 / 2, w2 / 2, 3))
for (ch in 1:3) { m <- rgbc[1:h2, 1:w2, ch]
  ds[, , ch] <- (m[seq(1, h2, 2), seq(1, w2, 2)] + m[seq(2, h2, 2), seq(1, w2, 2)] + m[seq(1, h2, 2), seq(2, w2, 2)] + m[seq(2, h2, 2), seq(2, w2, 2)]) / 4 }
rm(rgbc, m); invisible(gc())
cw <- 3094; chh <- 2140; Wf <- 2.40; Hf <- Wf * chh / cw
fx0 <- 1.92; ftop <- 1.42; fy0 <- ftop - Hf
CALL <- data.frame(lab = c("E2F1/E2F4\u2013p53\nhub", "CDK1", "CDC25A", "p21", "CDK2", "Rb", "HDAC1", "PPAR-\u03b3", "E2F5"),
                   px = c(595, 400, 483, 535, 395, 880, 848, 1045, 660), py = c(560, 322, 285, 252, 92, 752, 990, 690, 890),
                   tx = c(1130, 300, 560, 720, 230, 1010, 1000, 1150, 480), ty = c(440, 230, 170, 210, 60, 840, 1000, 600, 960))
CALL <- CALL |> mutate(x = fx0 + px / 1500 * Wf, y = ftop - py / 1500 * Wf, lx = fx0 + tx / 1500 * Wf, ly = ftop - ty / 1500 * Wf)
yD_sub <- 1.60
pD <- ggplot() +
  annotation_raster(ds, xmin = fx0, xmax = fx0 + Wf, ymin = fy0, ymax = ftop, interpolate = TRUE) +
  edge_layers(EDg, head_len = 1.2) + node_layer(node_df(ndD)) +
  bg_text(labD, aes(x, y, label = lab, hjust = hj, vjust = vj, fontface = face), pad = 0.3, alpha_bg = 0.8) +
  geom_segment(data = CALL, aes(x = lx, y = ly, xend = x, yend = y), colour = INK, linewidth = lwp(0.45)) +
  geom_label(data = CALL, aes(lx, ly, label = lab), size = FS6, family = FONT, fontface = "bold", lineheight = 0.9, colour = INK,
             fill = alpha("white", 0.92), label.size = lwp(0.3), border.colour = INK2, linewidth = lwp(0.3),
             label.padding = unit(1, "pt"), label.r = unit(0.8, "pt")) +
  annotate("text", x = 0.08, y = yD_sub, label = "Male-specific (Oct-3/4\u2013SMAD2/3 hub)", hjust = 0, vjust = 1,
           size = pt(6.5), fontface = "bold", family = FONT, colour = INK) +
  annotate("text", x = 2.28, y = yD_sub, label = "Female-specific (E2F1/E2F4\u2013p53 hub)", hjust = 0, vjust = 1,
           size = pt(6.5), fontface = "bold", family = FONT, colour = INK) +
  annotate("text", x = 0.08, y = -0.17, label = "'node fill: HPC male '*log[2]*'FC (pink up, blue down);'", parse = TRUE,
           hjust = 0, vjust = 1, size = FS6, family = FONT, colour = INK2) +
  annotate("text", x = 0.08, y = -0.285, label = "dashed: added by MetaCore, not measured",
           hjust = 0, vjust = 1, size = FS6, family = FONT, colour = INK2) +
  annotate("text", x = 2.28, y = fy0 - 0.03, label = "MetaCore-added objects were not measured",
           hjust = 0, vjust = 1, size = FS6, family = FONT, colour = INK2) +
  labs(title = "MetaCore network expansions, HPC sex-specific genes")
pD <- canvas(pD, c(-0.08, 4.34), c(-0.40, 1.62))

# ============================================================== E: GA Boolean networks
gm <- read.csv(file.path(SD, "Fig4E_GA_male_edges.csv"), stringsAsFactors = FALSE)
gf <- read.csv(file.path(SD, "Fig4E_GA_female_edges.csv"), stringsAsFactors = FALSE)
for (nm in c("gm", "gf")) { d <- get(nm); d$differential <- d$differential %in% c("True", "TRUE", TRUE); assign(nm, d) }
CODE <- c("1" = "Activation", "-1" = "Inhibition", "2" = "Unspecified")
expr_of <- function(d) { e <- c(setNames(d$expr_src, d$src), setNames(d$expr_tgt, d$tgt)); e[!duplicated(names(e))] }
ga_nodes <- function(pos, ex, r_pt) pos |> mutate(r = r_pt[name] * P2I, s = ex[name],
  fill = ifelse(is.na(s), "white", FILL[as.character(s)]), col = ifelse(is.na(s), INK, INK2),
  lw = ifelse(is.na(s), 0.8, 0.4), lt = ifelse(is.na(s), DASH, "solid"))
# male: hand layout (reference data units -> inches)
sm <- min(1.85 / 2.7, 2.0 / 2.93); mX <- function(x) 1.025 + (x - 0.05) * sm; mY <- function(y) -1.30 + (y - 0.015) * sm
pm <- data.frame(name = c("RUNX1", "FOXP1", "ESRRB", "PRDM11", "CXCL12", "TATDN3", "DHRS9", "PFKFB2"),
                 x0 = c(-0.25, 0.65, 0.2, 0.2, 1.05, 1.05, -0.95, -0.95), y0 = c(0.2, 0.2, -0.45, 0.85, -0.5, 0.85, 0.85, -0.45))
pm <- transform(pm, x = mX(x0), y = mY(y0))
stopifnot(setequal(pm$name, unique(c(gm$src, gm$tgt))))
exm <- expr_of(gm); rm_ <- setNames(rep(13, nrow(pm)), pm$name)
gm$eff <- CODE[as.character(gm$code)]; gm$eff[gm$differential & gm$diff_type == "--|"] <- "Inhibition"
gmd <- gm[gm$src != gm$tgt, ]; gms <- gm[gm$src == gm$tgt, ]
edEm <- gmd |> transmute(u = src, v = tgt, effect = eff, hl = differential,
                         rad = mapply(function(s, t) if (any(gm$src == t & gm$tgt == s)) 0.18 else 0.05, src, tgt), lw = 0.7, alpha = 1)
EEm <- build_edges(edEm, pm, rm_, bar_pt = 3.0)
angm <- c(RUNX1 = 200, FOXP1 = 20, ESRRB = 225)
slm <- do.call(rbind, lapply(seq_len(nrow(gms)), function(i) { p <- pm[pm$name == gms$src[i], ]
  transform(self_loop_df(p$x, p$y, angm[[gms$src[i]]], 0.3 * sm, i), col = ECOL[[gms$eff[i]]]) }))
ndEm <- ga_nodes(pm, exm, rm_)
# female: RUNX1 star
sf <- 2.30 / 3.24; fX <- function(x) 5.85 + x * sf; fY <- function(y) -1.37 + y * sf
nodes_f <- sort(unique(c(gf$src, gf$tgt)), method = "radix")
others <- c("OTUD1", "SMC1A", setdiff(nodes_f, c("RUNX1", "OTUD1", "SMC1A")))
af <- pi / 2 - 2 * pi * ((seq_along(others) - 1) + 0.5) / length(others)
pf <- rbind(data.frame(name = "RUNX1", a = NA, x = fX(0), y = fY(0)),
            data.frame(name = others, a = af, x = fX(0.8 * cos(af)), y = fY(0.8 * sin(af))))
exf <- expr_of(gf); rf <- setNames(rep(3.4, nrow(pf)), pf$name); rf["RUNX1"] <- 8
gf$eff <- CODE[as.character(gf$code)]
gf <- gf[order(c("2" = 0, "1" = 1, "-1" = 2)[as.character(gf$code)]), ]
gfd <- gf[gf$src != gf$tgt, ]; gfs <- gf[gf$src == gf$tgt, ]
edEf <- gfd |> transmute(u = src, v = tgt, effect = eff, hl = differential,
                         rad = ifelse(src == "FZR1", 0.25, ifelse(src == "SMC1A", 0.35, 0)),
                         lw = ifelse(eff == "Unspecified", 0.45, 0.7), alpha = 1)
EEf <- build_edges(edEf, pf, rf, bar_pt = 2.4)
slf <- do.call(rbind, lapply(seq_len(nrow(gfs)), function(i) { p <- pf[pf$name == gfs$src[i], ]
  ang <- if (gfs$src[i] == "RUNX1") 200 else p$a * 180 / pi; size <- if (gfs$src[i] == "RUNX1") 0.14 else 0.1
  transform(self_loop_df(p$x, p$y, ang, size * sf, i + 100), col = ECOL[[gfs$eff[i]]], hl = gfs$differential[i]) }))
ndEf <- ga_nodes(pf, exf, rf)
labEf <- pf[pf$name != "RUNX1", ] |> mutate(deg = a * 180 / pi, deg = ifelse(deg < -180, deg + 360, deg),
  right = deg >= -90 & deg <= 90, rot = ifelse(right, deg, deg - 180), hj = ifelse(right, 0, 1),
  lx = fX(0.8 * 1.11 * cos(a)), ly = fY(0.8 * 1.11 * sin(a)), lab = ifelse(name == "C5ORF15", "C5orf15", name))
nfm <- length(unique(c(gm$src, gm$tgt))); nfe <- nrow(gm); nff <- length(nodes_f); nfef <- nrow(gf)
lgE_n <- data.frame(name = 1:3, x = 2.36, y = -1.06 - (0:2) * 0.118, r = 2.6 * P2I, fill = c(FILL[["1"]], FILL[["-1"]], "white"),
                    col = c(INK2, INK2, INK), lw = c(0.4, 0.4, 0.8), lt = c("solid", "solid", DASH))
lgE_t <- data.frame(x = 2.46, y = -1.06 - (0:3) * 0.118,
                    lab = c("input +1 (up)", "input \u22121 (down)", "not in the expression input", "differential edge"))
pE <- ggplot() +
  # male
  edge_layers(EEm, head_len = 1.5) +
  geom_path(data = slm, aes(x, y, group = id, colour = col), linewidth = lwp(0.6)) +
  node_layer(node_df(ndEm)) +
  geom_text(data = pm, aes(x, y, label = name), size = FS6, family = FONT, fontface = "italic", colour = INK) +
  annotate("text", x = mX(-1.3), y = -0.24, label = sprintf("Male (%d nodes, %d edges; best fitness 6.0)", nfm, nfe),
           hjust = 0, vjust = 1, size = pt(6.5), fontface = "bold", family = FONT, colour = INK) +
  annotate("segment", x = 0.66, y = -2.04, xend = mX(0.60), yend = mY(-0.49), colour = MUTED, linewidth = lwp(0.4)) +
  annotate("text", x = mX(-1.3), y = -2.06 - (0:3) * 0.098, hjust = 0, vjust = 1, size = FS6, family = FONT, colour = INK2,
           label = c("'differential:'", "italic('ESRRB')~'\u2014|'~italic('CXCL12')", "'(unspecified in DB;'", "'GA-resolved)'"), parse = TRUE) +
  # female
  edge_layers(EEf, head_len = 1.1) +
  geom_path(data = slf[slf$hl, ], aes(x, y, group = id), colour = HL, alpha = 0.55, linewidth = lwp(3.6)) +
  geom_path(data = slf, aes(x, y, group = id, colour = col), linewidth = lwp(0.6)) +
  node_layer(node_df(ndEf)) +
  geom_text(data = labEf, aes(lx, ly, label = lab, angle = rot, hjust = hj), vjust = 0.5, size = FS6, family = FONT,
            fontface = "italic", colour = INK) +
  annotate("text", x = fX(0), y = fY(-0.2), label = "RUNX1", vjust = 1, size = FS6, family = FONT, fontface = "italic", colour = INK) +
  annotate("text", x = 4.70, y = -0.24, label = sprintf("Female (%d nodes, %d edges; best fitness 4.0)", nff, nfef),
           hjust = 0, vjust = 1, size = pt(6.5), fontface = "bold", family = FONT, colour = INK) +
  annotate("text", x = 4.78, y = -2.52, hjust = 0, vjust = 1, size = FS6, family = FONT, colour = INK2, parse = TRUE,
           label = "'differential: '*italic('SMC1A')~'\u2192'~italic('OTUD1')*'; '*italic('BCL10')~'self-loop'") +
  # shared legend
  node_layer(node_df(lgE_n)) +
  annotate("segment", x = 2.30, xend = 2.42, y = lgE_t$y[4], yend = lgE_t$y[4], colour = HL, alpha = 0.6, linewidth = lwp(3)) +
  geom_text(data = lgE_t, aes(x, y, label = lab), hjust = 0, size = FS6, family = FONT, colour = INK) +
  labs(title = "Genetic-algorithm Boolean networks, HPC sex-specific genes")
pE <- canvas(pE, c(0, 7.12), c(-2.66, -0.18))

# ============================================================== layout + export
row1 <- (pA | pB) + plot_layout(widths = c(3.24, 3.32))
row2 <- (pC | pD) + plot_layout(widths = c(2.66, 4.42))
fig <- add_tags(row1 / row2 / pE + plot_layout(heights = c(3.18, 2.02, 2.48)))
FIG_H <- 228
save_fig(fig, "Figure4", width_mm = 183, height_mm = FIG_H)

# ---- numeric checks (printed for the report)
cat(sprintf("B: %d nodes, %d directed pairs, self-loops not drawn: %s; effects: %s\n", nB_nodes, nB_pairs,
            paste(selfB$src, collapse = ","), paste(names(table(EBd$effect)), table(EBd$effect), collapse = " ")))
cat(sprintf("C: %d records, %d pairs (%d self-loops)\n", nC_records, nrow(EC), nrow(ECs)))
cat(sprintf("D: male edges %d, nodes %d\n", nrow(DM), nrow(posD)))
cat(sprintf("E: male %d nodes %d edges; female %d nodes %d edges; RUNX1 partners %d\n", nfm, nfe, nff, nfef,
            length(setdiff(unique(c(gf$tgt[gf$src == "RUNX1"], gf$src[gf$tgt == "RUNX1"])), "RUNX1"))))
cat("A measured nodes:", unique(measured_genes), "\n")
rm(ds); invisible(gc())
