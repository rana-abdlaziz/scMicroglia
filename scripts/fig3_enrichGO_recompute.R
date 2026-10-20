## Re-implementation of clusterProfiler::enrichGO (v4 algorithm; enricher_internal) because
## clusterProfiler is not installable here (no Bioconductor access; not in Ubuntu apt).
## Uses Ubuntu-packaged org.Hs.eg.db 3.18.0 + GO.db 3.18.0.
## Defaults as in pathway_script.R: keyType ENTREZID, pAdjustMethod "fdr" (=BH), pvalueCutoff 0.05,
## minGSSize 10, maxGSSize 500, universe = all genes annotated to the ontology (GOALL, propagated),
## readable = TRUE. qvalueCutoff 0.2 is not applied (Storey q <= BH, so non-binding at FDR<0.05).
suppressMessages({library(AnnotationDbi); library(org.Hs.eg.db); library(GO.db)})
out <- "/home/claude/scm/figures/out/source_data"
up  <- "/mnt/user-data/uploads/scMicroglia/microglia-subtypes"
hpc <- read.csv(file.path(up, "hipo_analysis/dge_filterd.csv"))
pfc <- read.csv(file.path(up, "pfc_result/dge_filterd.csv"))

kk <- keys(org.Hs.eg.db, keytype = "ENTREZID")
goAnno <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = kk, keytype = "ENTREZID",
                                                 columns = c("GOALL", "ONTOLOGYALL")))
goAnno <- unique(goAnno[!is.na(goAnno$GOALL), c("ENTREZID", "GOALL", "ONTOLOGYALL")])
goterm <- suppressMessages(AnnotationDbi::select(GO.db, keys = unique(goAnno$GOALL), columns = "TERM", keytype = "GOID"))
term_name <- setNames(goterm$TERM, goterm$GOID)
e2s <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = kk, keytype = "ENTREZID", columns = "SYMBOL"))
e2s <- setNames(e2s$SYMBOL, e2s$ENTREZID)

bitr_sym <- function(sym) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(sym), keytype = "SYMBOL", columns = "ENTREZID"))
  unique(m$ENTREZID[!is.na(m$ENTREZID)])
}

enrichGO2 <- function(gene, ont, minGS = 10, maxGS = 500, cut = 0.05, keep_all = FALSE) {
  a <- goAnno[goAnno$ONTOLOGYALL == ont, ]
  t2g <- split(a$ENTREZID, a$GOALL)
  universe <- unique(a$ENTREZID)
  g <- intersect(unique(gene), universe)
  n <- length(g); N <- length(universe)
  qa <- a[a$ENTREZID %in% g, ]
  qT <- unique(qa$GOALL)
  sz <- lengths(t2g[qT])
  qT <- qT[sz >= minGS & sz <= maxGS]
  if (!length(qT)) return(NULL)
  hits <- lapply(qT, function(t) intersect(g, t2g[[t]]))
  k <- lengths(hits); M <- lengths(t2g[qT])
  p <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  res <- data.frame(ID = qT, Description = term_name[qT], GeneRatio = paste0(k, "/", n),
                    BgRatio = paste0(M, "/", N), Count = k, pvalue = p,
                    p.adjust = p.adjust(p, method = "fdr"),
                    geneID = sapply(hits, function(h) paste(e2s[h], collapse = "/")),
                    stringsAsFactors = FALSE)
  res <- res[order(res$pvalue), ]
  if (!keep_all) res <- res[res$pvalue < cut & res$p.adjust < cut, ]
  rownames(res) <- NULL
  res
}

sets <- list(
  core109 = intersect(hpc$Gene_Symbols[hpc$DEG_Type == "gender-shared"], pfc$Gene_Symbols[pfc$DEG_Type == "gender-shared"]),
  HPC_dimorphic = hpc$Gene_Symbols[hpc$DEG_Type == "gender-dimorphic"],
  PFC_dimorphic = pfc$Gene_Symbols[pfc$DEG_Type == "gender-dimorphic"],
  HPC_shared_validation = hpc$Gene_Symbols[hpc$DEG_Type == "gender-shared"])
cat("set sizes:", sapply(sets, length), "\n")
for (s in names(sets)) {
  ent <- bitr_sym(sets[[s]])
  cat(s, ": symbols", length(sets[[s]]), "-> entrez", length(ent), "\n")
  for (o in c("BP", "MF")) {
    if (s == "HPC_shared_validation" && o == "MF") next
    r <- enrichGO2(ent, o, keep_all = (s == "HPC_shared_validation"))
    nsig <- if (is.null(r)) 0 else sum(r$p.adjust < 0.05)
    cat(sprintf("  %s: %d terms FDR<0.05\n", o, nsig))
    if (!is.null(r)) write.csv(r, file.path(out, sprintf("Fig3_recomputed_enrichGO_%s_%s.csv", s, o)), row.names = FALSE)
  }
}
writeLines(c(paste("org.Hs.eg.db", packageVersion("org.Hs.eg.db")), paste("GO.db", packageVersion("GO.db")),
             capture.output(org.Hs.eg.db)[grep("GOSOURCEDATE|EGSOURCEDATE", capture.output(org.Hs.eg.db))]),
           file.path(out, "Fig3_recompute_annotation_versions.txt"))
