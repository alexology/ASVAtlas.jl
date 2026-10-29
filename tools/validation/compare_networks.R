# Validation, step 4 (docs/VALIDATION.md): haplotype networks, ASVAtlas
# (ProfileAlignments.jl + HaplotypeNetworks.jl) against R (ape, pegas,
# DECIPHER), species by species; TCS is compared in tcs.R.
#   Rscript tools/validation/compare_networks.R WORKDIR
suppressPackageStartupMessages({
  library(ape); library(pegas); library(Biostrings); library(DECIPHER)
})
dir <- commandArgs(trailingOnly = TRUE)[1]
setwd(file.path(dir, "networks"))
summ <- read.delim("summary_julia.tsv", stringsAsFactors = FALSE)
files <- unique(summ$file)

# Columns with a gap or an ambiguous base in any sequence are left out, as in
# HaplotypeNetworks.jl and PopART.
mask <- function(m) {
  ok <- apply(m, 2, function(col) all(col %in% c("A", "C", "G", "T")))
  m[, ok, drop = FALSE]
}
read_matrix <- function(path) {
  s <- readDNAStringSet(path)
  m <- do.call(rbind, strsplit(as.character(s), ""))
  rownames(m) <- names(s)
  m
}
edges_key <- function(a, b, w) paste(pmin(a, b), pmax(a, b), w)

# Union of all minimum spanning trees: an edge u-v of weight w belongs to one
# iff u and v are not connected by edges lighter than w (Kruskal by levels).
msn_union <- function(D) {
  n <- nrow(D)
  comp <- seq_len(n)
  find <- function(x) { while (comp[x] != x) x <- comp[x]; x }
  out <- NULL
  for (w in sort(unique(D[upper.tri(D)]))) {
    ij <- which(D == w & upper.tri(D), arr.ind = TRUE)
    keep <- ij[mapply(function(i, j) find(i) != find(j), ij[, 1], ij[, 2]), , drop = FALSE]
    out <- rbind(out, cbind(keep, w))
    for (k in seq_len(nrow(keep))) {
      a <- find(keep[k, 1]); b <- find(keep[k, 2])
      if (a != b) comp[a] <- b
    }
  }
  out
}

rows <- list()
for (f in files) {
  sp <- summ$species[summ$file == f][1]
  cat("start", sp, format(Sys.time(), "%H:%M:%S"), "
")
  m <- read_matrix(paste0(f, ".fasta"))
  mem <- read.delim(paste0(f, ".members.tsv"), stringsAsFactors = FALSE)
  if (max(mem$haplotype) < 3) { cat("  skipped: fewer than 3 haplotypes
"); next }
  DJ <- as.matrix(read.delim(paste0(f, ".dist.tsv"), header = FALSE))
  mm <- mask(m)
  # Haplotypes in R: identical masked sequences.
  strs <- apply(mm, 1, paste, collapse = "")
  hapR <- match(strs, unique(strs))
  same_partition <- all(table(paste(hapR, mem$haplotype[match(rownames(mm), mem$id)])) > 0) &&
    length(unique(hapR)) == length(unique(mem$haplotype)) &&
    length(unique(paste(hapR, mem$haplotype[match(rownames(mm), mem$id)]))) == length(unique(hapR))
  # Distances between haplotypes, in the Julia numbering (one member each).
  nh <- max(mem$haplotype)
  repr <- mem$id[match(seq_len(nh), mem$haplotype)]
  hm <- mm[repr, , drop = FALSE]
  dn <- as.DNAbin(tolower(hm))
  DR <- as.matrix(dist.dna(dn, model = "N", pairwise.deletion = FALSE))
  dimnames(DR) <- NULL
  dist_equal <- all(DR == DJ)

  get_j <- function(method) {
    e <- read.delim(paste0(f, ".", method, ".edges.tsv"), stringsAsFactors = FALSE)
    e
  }
  hk <- function(x) as.integer(sub("^H", "", x))

  # MST: total length (unique even if trees differ).
  t0 <- Sys.time(); mstR <- pegas::mst(as.dist(DR)); t_mst <- as.numeric(Sys.time() - t0, units = "secs")
  mstJ <- get_j("mst")
  # MSN: union of all MSTs computed here, pegas::msn and Julia.
  U <- msn_union(DR)
  t0 <- Sys.time(); msnR <- pegas::msn(as.dist(DR)); t_msn <- as.numeric(Sys.time() - t0, units = "secs")
  msnJ <- get_j("msn")
  keyU <- edges_key(U[, 1], U[, 2], U[, 3])
  keyJ <- edges_key(hk(msnJ$from), hk(msnJ$to), msnJ$steps)
  keyP <- edges_key(msnR[, 1], msnR[, 2], msnR[, 3])

  row <- data.frame(species = sp, n = nrow(m), sites_julia = summ$sites[summ$file == f][1],
                    sites_R = ncol(mm), haplotypes_julia = nh, haplotypes_R = length(unique(hapR)),
                    same_haplotypes = same_partition, same_distances = dist_equal,
                    mst_steps_julia = sum(mstJ$steps), mst_steps_pegas = sum(mstR[, 3]),
                    msn_edges_julia = nrow(msnJ), msn_edges_union = nrow(U), msn_edges_pegas = nrow(msnR),
                    msn_julia_eq_union = setequal(keyJ, keyU),
                    msn_pegas_extra = length(setdiff(keyP, keyU)), msn_pegas_missing = length(setdiff(keyU, keyP)),
                    pegas_mst_s = round(t_mst, 3), pegas_msn_s = round(t_msn, 3))

  # MJN with pegas (slow on many haplotypes): nodes, median vectors, edges, length.
  mjnJ <- get_j("mjn")
  row$mjn_julia <- sprintf("%d mv, %d edges, %d steps",
                           summ$inferred[summ$file == f & summ$method == "mjn"], nrow(mjnJ), sum(mjnJ$steps))
  row$mjn_pegas <- NA
  if (nh <= 40) {
    t0 <- Sys.time()
    r <- tryCatch({
      net <- pegas::mjn(dn, epsilon = 0, quiet = TRUE)
      nmv <- length(attr(net, "labels")) - nh
      sprintf("%d mv, %d edges, %d steps (%.0f s)", nmv, nrow(net), sum(net[, 3]),
              as.numeric(Sys.time() - t0, units = "secs"))
    }, error = function(e) paste("error:", conditionMessage(e)))
    row$mjn_pegas <- r
    cat("  mjn pegas:", r, "
")
  }

  # TCS with haplotypes::parsimnet: networks (components) and connection limit.
  tcsJ <- summ[summ$file == f & summ$method == "tcs", ]
  row$tcs_julia <- sprintf("%d components, %d intermediates, %d edges", tcsJ$components, tcsJ$inferred, tcsJ$edges)
  # TCS of the R package haplotypes: in tcs.R, a separate process (once
  # loaded, haplotypes makes pegas::msn hang).

  # Independent alignment: DECIPHER on the same oriented sequences.
  raw <- DNAStringSet(gsub("-", "", apply(m, 1, paste, collapse = "")))
  names(raw) <- rownames(m)
  if (length(unique(width(raw))) == 1 && all(!grepl("-", apply(m, 1, paste, collapse = "")))) {
    row$decipher <- "no gaps: alignment trivial"
  } else {
    al <- AlignSeqs(raw, verbose = FALSE, processors = 4)
    dm <- do.call(rbind, strsplit(as.character(al), ""))
    rownames(dm) <- names(al)
    dm <- mask(dm)[rownames(m), , drop = FALSE]
    hd <- as.matrix(dist.dna(as.DNAbin(tolower(dm[repr, , drop = FALSE])), model = "N"))
    dimnames(hd) <- NULL
    s2 <- apply(dm, 1, paste, collapse = "")
    row$decipher <- sprintf("%d sites, %d haplotypes, distances equal %s, max |diff| %d",
                            ncol(dm), length(unique(s2)), all(hd == DJ), max(abs(hd - DJ)))
  }
  rows[[f]] <- row
  cat("done", sp, "\n")
}
out <- do.call(rbind, rows)
write.table(out, "report_networks.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250)
print(out[, c(1:16, 17, 18)], row.names = FALSE)
print(out[, c("species", "mjn_julia", "mjn_pegas")], row.names = FALSE)
print(out[, c("species", "decipher")], row.names = FALSE)
