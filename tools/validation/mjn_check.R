# Validation, step 6 (docs/VALIDATION.md): MJN, pegas against HaplotypeNetworks.jl.
#   Rscript tools/validation/mjn_check.R WORKDIR For the observed haplotypes, the
# shortest path in the network is compared with the number of differences
# (in a median-joining network the two should agree as far as possible:
# every extra step is a distortion). Also times pegas::mjn.
suppressPackageStartupMessages({ library(ape); library(pegas); library(Biostrings) })
setwd(file.path(commandArgs(trailingOnly = TRUE)[1], "networks"))
summ <- read.delim("summary_julia.tsv", stringsAsFactors = FALSE)
path_lengths <- function(from, to, w, nodes) {
  n <- length(nodes); D <- matrix(Inf, n, n, dimnames = list(nodes, nodes)); diag(D) <- 0
  for (k in seq_along(from)) D[from[k], to[k]] <- D[to[k], from[k]] <- min(D[from[k], to[k]], w[k])
  for (k in seq_len(n)) D <- pmin(D, outer(D[, k], D[k, ], "+"))  # Floyd-Warshall
  D
}
# pegas::mjn is run on the species with at most 40 haplotypes.
sizes <- tapply(summ$haplotypes, summ$file, max)
for (f in names(sizes)[sizes >= 3 & sizes <= 40]) {
  sp <- summ$species[summ$file == f][1]
  s <- readDNAStringSet(paste0(f, ".fasta"))
  m <- do.call(rbind, strsplit(as.character(s), "")); rownames(m) <- names(s)
  m <- m[, apply(m, 2, function(col) all(col %in% c("A", "C", "G", "T"))), drop = FALSE]
  mem <- read.delim(paste0(f, ".members.tsv"), stringsAsFactors = FALSE)
  nh <- max(mem$haplotype)
  hm <- m[mem$id[match(seq_len(nh), mem$haplotype)], , drop = FALSE]
  rownames(hm) <- paste0("H", seq_len(nh))
  dn <- as.DNAbin(tolower(hm))
  H <- as.matrix(dist.dna(dn, model = "N"))
  # Julia
  e <- read.delim(paste0(f, ".mjn.edges.tsv"), stringsAsFactors = FALSE)
  nodes <- unique(c(rownames(hm), e$from, e$to))
  P <- path_lengths(e$from, e$to, e$steps, nodes)[rownames(hm), rownames(hm)]
  u <- upper.tri(H)
  cat(sprintf("%-26s Julia: %3d mv, %3d steps; pairs with path = distance %d/%d, extra steps %d\n",
              sp, length(nodes) - nh, sum(e$steps), sum(P[u] == H[u]), sum(u), sum(P[u] - H[u])))
  # pegas
  t0 <- Sys.time()
  net <- tryCatch(pegas::mjn(dn, epsilon = 0, quiet = TRUE, max.n.cost = 1e6), error = function(e) NULL)
  tt <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(net)) { cat(sprintf("%-26s pegas: error (%.1f s)\n", "", tt)); next }
  lab <- attr(net, "labels")
  Pp <- path_lengths(lab[net[, 1]], lab[net[, 2]], net[, 3], lab)[rownames(hm), rownames(hm)]
  cat(sprintf("%-26s pegas: %3d mv, %3d steps; pairs with path = distance %d/%d, extra steps %d (%.1f s)\n",
              "", length(lab) - nh, sum(net[, 3]), sum(Pp[u] == H[u]), sum(u), sum(Pp[u] - H[u]), tt))
  flush.console()
}
