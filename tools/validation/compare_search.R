# Validation, step 3 (docs/VALIDATION.md): similarity, ASVAtlas
# (NucleotideSearch.jl) against pwalign (R) and BLAST.
#   Rscript tools/validation/compare_search.R WORKDIR
# pwalign aligns each query globally into a local stretch of the reference
# ("global-local"), with the scoring of NucleotideSearch/VSEARCH: match 2,
# mismatch -4, gap of length n = 20 + 2n. Identity = matches / alignment
# columns (VSEARCH iddef 2), coverage = aligned query bases / query length.
suppressPackageStartupMessages({ library(Biostrings); library(pwalign) })
dir <- commandArgs(trailingOnly = TRUE)[1]
setwd(dir)
ref <- readDNAStringSet("blastdb/ref.fasta")
names(ref) <- sub(" .*", "", names(ref))
qs <- readDNAStringSet("queries.fasta")
names(qs) <- sub(" .*", "", names(qs))

julia_app <- read.delim("julia_hits_app.tsv", stringsAsFactors = FALSE)
julia_ex <- read.delim("julia_hits_exhaustive.tsv", stringsAsFactors = FALSE)
blast <- read.delim("blast_hits.tsv", header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("query", "reference", "pident", "length", "mismatch", "gapopen",
                                  "qstart", "qend", "sstart", "send", "sstrand", "qlen"))

# Candidates: every pair found by Julia, and BLAST pairs that could reach 90%.
cand <- unique(rbind(
  julia_ex[, c("query", "reference")],
  blast[blast$pident >= 85 & blast$length >= 0.8 * blast$qlen, c("query", "reference")]))
cat("candidate pairs:", nrow(cand), "\n")

mat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -4, baseOnly = TRUE)
align_pairs <- function(q, r) {
  # Both strands; the better score wins (as in a search on both strands).
  best <- NULL
  for (s in c("+", "-")) {
    qq <- if (s == "+") q else reverseComplement(q)
    a <- pairwiseAlignment(qq, r, type = "global-local", substitutionMatrix = mat,
                           gapOpening = 20, gapExtension = 2)
    if (is.null(best) || score(a) > best$score) {
      m <- nmatch(a); mm <- nmismatch(a)
      best <- list(score = score(a), identity = 100 * m / nchar(a),
                   coverage = (m + mm) / length(q), strand = s)
    }
  }
  best
}
t0 <- Sys.time()
res <- do.call(rbind, lapply(seq_len(nrow(cand)), function(k) {
  a <- align_pairs(qs[[cand$query[k]]], ref[[cand$reference[k]]])
  data.frame(query = cand$query[k], reference = cand$reference[k], identity = a$identity,
             coverage = a$coverage, strand = a$strand)
}))
cat("pwalign:", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
write.table(res, "pwalign_hits.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

key <- function(d) paste(d$query, d$reference)
report <- NULL
for (t in c(100, 99, 98, 97, 95, 90)) {
  truth <- res[res$identity >= t - 1e-9 & res$coverage >= 0.9, ]
  for (lab in c("app", "exhaustive")) {
    j <- if (lab == "app") julia_app else julia_ex
    j <- j[j$identity >= t - 1e-9, ]
    report <- rbind(report, data.frame(
      threshold = t, julia = lab, pwalign_pairs = nrow(truth), julia_pairs = nrow(j),
      both = length(intersect(key(truth), key(j))),
      only_pwalign = length(setdiff(key(truth), key(j))),
      only_julia = length(setdiff(key(j), key(truth))),
      queries_pwalign = length(unique(truth$query)), queries_julia = length(unique(j$query))))
  }
}
print(report, row.names = FALSE)
write.table(report, "report_search.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Identity and strand of the pairs found by both.
m <- merge(julia_ex, res, by = c("query", "reference"), suffixes = c(".julia", ".pwalign"))
d <- m$identity.julia - m$identity.pwalign
cat("\npairs compared:", nrow(m), "\n")
cat("identical identity (|diff| < 1e-6):", sum(abs(d) < 1e-6), "\n")
cat("max |identity difference|:", max(abs(d)), "\n")
cat("same strand:", sum(m$strand.julia == m$strand.pwalign), "\n")
cat("max |coverage difference|:", max(abs(m$coverage.julia - m$coverage.pwalign)), "\n")
if (any(abs(d) >= 1e-6)) {
  cat("\nlargest differences:\n")
  print(head(m[order(-abs(d)), c("query", "reference", "identity.julia", "identity.pwalign",
                                  "coverage.julia", "coverage.pwalign")], 10), row.names = FALSE)
}

# BLAST (local) identity for the same pairs: comparable only when BLAST
# aligns the whole query.
mb <- merge(julia_ex, blast, by = c("query", "reference"))
full <- mb$length == mb$qlen & mb$qstart == 1 & mb$qend == mb$qlen
cat("\nBLAST pairs aligning the whole query:", sum(full), "of", nrow(mb), "\n")
cat("  same identity (rounded to 2 decimals as BLAST):",
    sum(abs(round(mb$identity[full], 2) - mb$pident[full]) < 0.006), "\n")

# Exhaustive check outside the candidates: 10 queries against all references.
set.seed(1)
pick <- sample(names(qs), 10)
t0 <- Sys.time()
miss <- 0
for (q in pick) {
  for (s in c("+", "-")) {
    qq <- if (s == "+") qs[[q]] else reverseComplement(qs[[q]])
    a <- pairwiseAlignment(rep(DNAStringSet(qq), length(ref)), ref, type = "global-local",
                           substitutionMatrix = mat, gapOpening = 20, gapExtension = 2)
    idn <- 100 * nmatch(a) / nchar(a)
    cov <- (nmatch(a) + nmismatch(a)) / length(qq)
    hit <- names(ref)[idn >= 90 & cov >= 0.9]
    if (length(hit)) miss <- miss + length(setdiff(paste(q, hit), key(cand)))
  }
}
cat("\nall-against-all for", length(pick), "queries:", round(as.numeric(Sys.time() - t0, units = "secs")),
    "s; pairs >= 90% outside the candidates:", miss, "\n")
