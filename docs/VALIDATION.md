# Validation of ASVAtlas against R and BLAST

ASVAtlas was compared with independent tools on a real reference: 8,198 COI
ASVs (418 bp) of freshwater invertebrates with 21,383 georeferenced
occurrences. The reference itself is private and not part of the package;
the scripts in [`tools/validation/`](../tools/validation) repeat the
comparison on any reference in the ASVAtlas format (see "Reproducing" at
the end). The comparison was run on 2026-10-28.

In short:

- **Similarity search:** ASVAtlas found exactly the same reference ASVs as an
  exhaustive alignment in R at every similarity level, with identical
  identity, coverage and strand. No pair was missed, also with the early stop
  used by the web page.
- **Haplotype networks:** haplotypes, distances, minimum spanning trees,
  minimum spanning networks and the separate TCS networks agree with R. The
  median-joining networks differ from `pegas::mjn` (a known difference);
  those of ASVAtlas keep the distances between observed haplotypes better.
- **Speed:** the network algorithms of ASVAtlas take milliseconds where R
  takes seconds to minutes; in the web page most of the time of a network
  goes to the multiple alignment.

## Software

| Tool | Version | Used for |
|---|---|---|
| Julia | 1.12.2 | ASVAtlas 0.1.0, NucleotideSearch.jl 0.1.0, ProfileAlignments.jl 0.3.0, HaplotypeNetworks.jl 0.1.0 |
| BLAST+ | 2.17.0 | independent candidate pairs (`blastn`) |
| R | 4.5.0 | |
| pwalign / Biostrings | 1.4.0 / 2.76.0 | reference pairwise alignments |
| ape / pegas | 5.8.1 / 1.3 | distances, MST, MSN, MJN |
| haplotypes | 1.1.3.1 | TCS (statistical parsimony, `parsimnet`) |
| DECIPHER | 3.4.0 | independent multiple alignment |

## Similarity search

### Queries

215 query ASVs, generated with a fixed random seed:

- 5 test queries (identical to a reference ASV, one mutation, reverse
  strand, about 97%, unknown);
- 200 reference ASVs with 0–40 random substitutions (100% to about 90%
  identity); every 7th with 1–3 indels and every 6th reverse-complemented
  (33);
- 10 shuffled reference sequences, which should match nothing.

### Method

ASVAtlas searched at 90% identity, query coverage 0.9, both strands, up to
2,000 matches per query, in two ways: as the web page does (the search of a
query stops after 64 rejected candidates) and exhaustively (no early stop).

R (pwalign) realigned every candidate pair: those found by ASVAtlas and those
found by BLAST with at least 85% local identity over at least 80% of the
query (16,920 pairs). The query was aligned whole into a stretch of the
reference (`type = "global-local"`), on both strands, with the scoring of
NucleotideSearch and VSEARCH (match 2, mismatch −4, gap of length *n*
20 + 2*n*). Identity is matches / alignment columns, coverage aligned query
bases / query length, as in ASVAtlas.

To check that no pair lies outside the candidates, 10 random queries were
also aligned in R against all 8,198 reference ASVs on both strands.

### Results

| Similarity | Pairs (R) | Pairs (ASVAtlas, page) | Pairs (ASVAtlas, exhaustive) | Queries with matches | Differences |
|---|---|---|---|---|---|
| 100% | 7 | 7 | 7 | 7 | 0 |
| ≥ 99% | 72 | 72 | 72 | 23 | 0 |
| ≥ 98% | 237 | 237 | 237 | 37 | 0 |
| ≥ 97% | 532 | 532 | 532 | 55 | 0 |
| ≥ 95% | 1,477 | 1,477 | 1,477 | 92 | 0 |
| ≥ 90% | 5,333 | 5,333 | 5,333 | 204 | 0 |

- Identity, coverage and strand are identical for all 5,333 pairs (largest
  difference 0).
- BLAST gives the same identity (rounded to two decimals) for the 3,521 pairs
  where it aligns the whole query; for the others it clips the ends of the
  query (local alignment), so its identity is not comparable.
- The all-against-all alignment of 10 queries found no pair at ≥ 90%
  that ASVAtlas had not found.
- The 10 shuffled sequences matched nothing; the 204 queries with matches are
  all the others but the unknown test query.

| Search | Time (215 queries) |
|---|---|
| ASVAtlas, as in the page | **7 s** |
| ASVAtlas, exhaustive | 271 s |
| BLAST (`blastn`, 4 threads) | 340 s |
| pwalign on the 16,920 candidate pairs | 343 s |

## Haplotype networks

### Species and method

11 species: the four species of the test queries (*Cricotopus bicinctus*,
*Gammarus pulex*, *Potamophylax latipennis*, *Simulium tuberosum*), the three
largest other species with at most 300 reference ASVs (*Stenacron
interpunctatum*, *Serratella ignita*, *Simulium venustum complex*) and four
species with 15 reference ASVs (*Amphinemura sulcicollis*, *Baetis
buceratus*, *Chironomus harpi*, *Euglesa casertana*); 3 to 236 haplotypes.

For each species ASVAtlas built the MST, MSN, MJN and TCS networks as the web
page does (orientation, alignment with ProfileAlignments.jl, haplotypes and
network with HaplotypeNetworks.jl) and wrote its alignment. R started from
the same alignment, with the same rule: alignment columns with a gap or an
ambiguous base in any sequence are left out (as in PopART).

### Results

| Check | Result |
|---|---|
| Haplotypes (which sequences are identical) | identical, 11/11 species |
| Distances between haplotypes (`ape::dist.dna`, `model = "N"`) | identical, 11/11 |
| Length of the minimum spanning tree (`pegas::mst`) | identical, 11/11 (e.g. 532 steps for *C. bicinctus*) |
| Minimum spanning network | identical to the union of all minimum spanning trees computed in R, 11/11 |
| TCS: number of separate networks (`haplotypes::parsimnet`, limit 8 steps) | identical, 10/10 species finished in R |
| Median-joining network (`pegas::mjn`, ε = 0) | different (see below) |

**Minimum spanning network.** The MSN is the union of all minimum spanning
trees. ASVAtlas returns exactly that set of links in every species.
`pegas::msn` 1.3 returned fewer links in 10 of 11 species (1 to 188 missing,
e.g. 235 instead of 423 for *C. bicinctus*), and none outside the union.

**Median-joining network.** `pegas::mjn` was run on the 6 species with at most
40 haplotypes; with its default settings it stops on *G. pulex* ("too many
connection costs computed"), so `max.n.cost` was raised. The two programs
infer different median vectors (also on the worked examples of Bandelt et
al. 1999, see HaplotypeNetworks.jl). As an objective measure, the shortest
path in the network between two observed haplotypes was compared with their
number of differences: every extra step is a distortion.

| Species (haplotypes) | ASVAtlas: median vectors, steps | pairs with path = distance | extra steps | pegas: median vectors, steps | pairs with path = distance | extra steps |
|---|---|---|---|---|---|---|
| *Gammarus pulex* (34) | 8, 108 | 198/561 | 2,393 | 26, 122 | 113/561 | 3,128 |
| *Potamophylax latipennis* (3) | 3, 27 | 3/3 | 0 | 3, 26 | 2/3 | 1 |
| *Amphinemura sulcicollis* (15) | 4, 49 | 92/105 | 22 | 7, 56 | 47/105 | 196 |
| *Baetis buceratus* (15) | 4, 68 | 68/105 | 144 | 4, 66 | 45/105 | 288 |
| *Chironomus harpi* (15) | 0, 18 | 97/105 | 16 | 3, 19 | 60/105 | 140 |
| *Euglesa casertana* (15) | 10, 102 | 41/105 | 584 | 11, 112 | 19/105 | 2,116 |

In every species the ASVAtlas network keeps more distances exact, with fewer
extra steps. In *P. latipennis* one site has three different bases in the
three sequences: ASVAtlas joins three median vectors in a triangle (one step
each), so that every pair of haplotypes keeps its real distance; pegas uses
one fewer link and lengthens one path by a step.

### Speed of the network algorithms

Seconds for the algorithm alone, on the same aligned sequences (ASVAtlas: best
of 3 after compilation).

| Species (haplotypes) | MST ASVAtlas / pegas | MSN ASVAtlas / pegas | MJN ASVAtlas / pegas | TCS ASVAtlas / haplotypes |
|---|---|---|---|---|
| *Cricotopus bicinctus* (236) | 0.002 / 0.002 | 0.002 / 2.2 | 0.064 / not run | 0.42 / stopped after 300 |
| *Stenacron interpunctatum* (124) | 0.0006 / 0.001 | 0.0005 / 0.12 | 0.003 / not run | 0.057 / 16.8 |
| *Simulium tuberosum* (116) | 0.0004 / 0.001 | 0.0004 / 0.12 | 0.007 / not run | 0.026 / 7.2 |
| *Serratella ignita* (104) | 0.0004 / < 0.001 | 0.0004 / 0.10 | 0.010 / not run | 0.011 / 1.2 |
| *Simulium venustum complex* (101) | 0.0003 / < 0.001 | 0.0003 / 0.05 | 0.008 / not run | 0.037 / 93.0 |
| *Gammarus pulex* (34) | 0.0001 / < 0.001 | < 0.001 / 0.004 | 0.0004 / 5.1 | 0.008 / 0.11 |
| 15 haplotypes (4 species) | < 0.001 / < 0.001 | < 0.001 / 0.001 | 0.0002 / ≤ 0.1 | 0.008 / 0.03–0.37 |

In the web page a network takes a few seconds (about 6 s for *C. bicinctus*):
almost all of it is the multiple alignment of the sequences.

### Alignment

DECIPHER (`AlignSeqs`) realigned the sequences of every species. In 10 species
all sequences have 418 bases and both programs align them without gaps. In
*C. bicinctus* (236 ASVs, all of 418 bases) ProfileAlignments gives 10 ASVs of
a lineage about 11% divergent from the others a 3-base deletion (columns
290–292) and a 3-base insertion 14 columns later (304–306); DECIPHER aligns
them without gaps, with 6–7 more substitutions in that T-rich stretch. Both
alignments are plausible. Since columns with gaps are left out, ASVAtlas uses
415 sites instead of 418 there, and distances between the two lineages are up
to 6 steps shorter than with the DECIPHER alignment; distances within each
lineage are not affected.

### Note on R

Once the R package haplotypes is loaded, `pegas::msn` does not return (it
takes 0.3 s on the same data in a fresh session). The TCS comparison is
therefore run in a separate R process.

## Reproducing

The scripts write their files to a working folder (`WORKDIR`). These files
contain reference sequences: keep the folder private when the reference is.

```sh
# 1. Queries, ASVAtlas searches and networks (Julia)
julia -t 4 tools/validation/export_julia.jl reference.fasta metadata.tsv WORKDIR \
      [test_queries.fasta] [species ...]
# 2. BLAST candidates (makeblastdb and blastn on the PATH)
sh tools/validation/run_blast.sh reference.fasta WORKDIR
# 3-6. Comparisons in R (Biostrings, pwalign, ape, pegas, haplotypes, DECIPHER)
Rscript tools/validation/compare_search.R WORKDIR
Rscript tools/validation/compare_networks.R WORKDIR
Rscript tools/validation/tcs.R WORKDIR
Rscript tools/validation/mjn_check.R WORKDIR
# 7. Time of the network algorithms alone (Julia)
julia tools/validation/time_julia.jl WORKDIR
```

Without species names, step 1 takes the three largest species with at most
300 reference ASVs and the four with the number of ASVs nearest to 15. The
results above were obtained naming the 11 species listed in "Species and
method". The scripts also run on the synthetic example of the package
(`examples/data`), where every check agrees.
