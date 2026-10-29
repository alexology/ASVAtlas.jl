# Validation, step 5 (docs/VALIDATION.md): TCS with the R package haplotypes
# (parsimnet), on the alignments ASVAtlas used; at most 300 s per species.
# A separate process: once loaded, haplotypes makes pegas::msn hang.
#   Rscript tools/validation/tcs.R WORKDIR
suppressPackageStartupMessages({ library(Biostrings); library(haplotypes) })
setwd(file.path(commandArgs(trailingOnly = TRUE)[1], "networks"))
summ <- read.delim("summary_julia.tsv", stringsAsFactors = FALSE)
sizes <- tapply(summ$haplotypes, summ$file, max)
for (f in names(sort(sizes))) {
  s <- readDNAStringSet(paste0(f, ".fasta"))
  m <- do.call(rbind, strsplit(as.character(s), "")); rownames(m) <- names(s)
  m <- m[, apply(m, 2, function(col) all(col %in% c("A", "C", "G", "T"))), drop = FALSE]
  mem <- read.delim(paste0(f, ".members.tsv"), stringsAsFactors = FALSE)
  nh <- max(mem$haplotype)
  hm <- m[mem$id[match(seq_len(nh), mem$haplotype)], , drop = FALSE]
  if (nh < 2) next
  j <- summ[summ$file == f & summ$method == "tcs", ]
  t0 <- Sys.time()
  r <- tryCatch({
    setTimeLimit(elapsed = 300, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
    p <- parsimnet(as.dna(hm), indels = "missing")
    sprintf("%d components, limit %d steps", length(p@d), p@conlimit[1])
  }, error = function(e) paste("stopped:", conditionMessage(e)))
  setTimeLimit(elapsed = Inf)
  cat(sprintf("%s (%d haplotypes) | Julia: %d components, %.2f s | haplotypes: %s, %.2f s
", j$species, nh,
              j$components, j$seconds, r, as.numeric(Sys.time() - t0, units = "secs")))
  flush.console()
}
