# ASVAtlas.jl

Search amplicon sequence variants (ASVs) against a reference set of ASVs whose
taxonomy and sampling locations are already known, at a chosen percent
identity (100%, 99%, 97%, ...), and show where the matching reference ASVs
were found on a world map. A web application built on
[Genie.jl](https://github.com/GenieFramework/Genie.jl) and a Julia API do the
same search.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

ASVAtlas does **not** assign taxonomy: each hit carries the taxonomy already
attached to the reference ASV. When the best hits of a query carry different
names, all of them are shown. The search is done by
[NucleotideSearch.jl](https://github.com/alexology/NucleotideSearch.jl) (pure
Julia, no BLAST or VSEARCH).

Julia ≥ 1.10. Not registered; three dependencies are not registered either:
NucleotideSearch.jl, HaplotypeNetworks.jl and ProfileAlignments.jl. The
`[sources]` entries of `Project.toml` need Julia ≥ 1.11 and the three
folders next to `ASVAtlas.jl`. On Julia 1.10, `Pkg.develop` the three
packages first.

## Quick start

```julia
using ASVAtlas

ref = load_reference("reference.fasta"; metadata = "reference_metadata.tsv")
serve(ref)                       # open http://127.0.0.1:8080 (Ctrl-C to stop)

server = serve(ref; async = true) # or in the background...
shutdown(server)                  # ...and stop it
```

The example in `examples/` uses synthetic data (invented taxa, sequences and
sites, written by `tools/make_example_data.jl`):

```
julia -t auto examples/run_app.jl        # or include("examples/run_app.jl") in a REPL; then paste examples/data/queries.fasta
```

The first tab answers **"Where in the world is my ASV?"**. Before a search
it shows only a form: paste your ASVs (FASTA, or one sequence per line), read
or drop a FASTA file, or press **Try an example** (when the server offers
one), then **Find in the world**. The search starts with identical records
(100%); query coverage, number of matches and strands are under **Advanced
options**.

After the search the map comes first. Above it, the similarity (100, 99, 98,
97, 95 or 90%) widens the search at once, and a sentence gives the answer,
e.g. "4 of your 5 ASVs are found at 11 places in Canada, Germany and United
Kingdom: 76 records ≥ 97% similar". **Edit ASVs** goes back to the form with the
same ASVs and coordinates; **New search** empties the form, the coordinates,
the results and the haplotype network, and sets the similarity back to 100%
(in the form, **Clear** empties sequences and coordinates).
Each sequence must appear once in the FASTA: the same sequence under two
identifiers is refused, with the identifiers in the message.
Under **Sampling locations (optional)** a CSV file gives the places of the
ASVs: a header with exactly the columns `asv_id`, `latitude` and `longitude`,
then one row per place, so that an ASV found at several places has several
rows. `asv_id` is the identifier in the FASTA; coordinates are decimal
degrees, WGS84 (e.g. `my_ASV,47.393030,8.544811`; a decimal comma is accepted
in semicolon- or tab-separated files). Other columns, degrees-minutes-seconds,
out-of-range or swapped coordinates are refused with a message. The file is
read in the browser and never sent to the server. Each sampling place is a gold
star on the map (distinct by shape and colour from the records); the line of
an ASV gives the distance from its nearest sampling place to its nearest
matching record, the tooltip of a star the distance from that place; each
place of the selected ASV is joined to its nearest record by a dashed line;
the results per ASV get the columns `n_sampling_places`, `sampling_places`
and `nearest_record_km`. When a table is loaded, an ASV without a row in it
says "no sampling place in your table" (and the node of the network that holds
it names it as not on the map), so that a missing row is not taken for a
fault of the page. In the haplotype network tab each star is filled with the colour of the
group of its ASVs (orange within the chosen number of mutations from the
centre, blue farther; split in halves when the ASVs sampled at a place differ)
and keeps a gold border; the tooltip gives the haplotype and the distance from
the centre of each ASV. Below the
map there is one line per ASV: the species of its best matches, its places
and countries, records and best similarity, and a button to its haplotype
network; an ASV without matches offers the next, wider similarity ("Try
95%"). Selecting an ASV keeps only its records on the map and shows its
matches grouped by taxon (number of reference ASVs, similarity range,
records, lineage); coverage, strand and position on the reference are under
**Technical details**. The countries come from the bundled world map (none
for places at sea or on coasts left out by its outlines). The map never zooms
closer than a region, so that the place can be recognised. Markers are
coloured by the best identity of the hits at that place, on a single aqua
scale (brighter = more similar: 100%, ≥ 99%, ≥ 97%, lower), checked for
lightness order and contrast against the land and sea colours; hovering shows
the country, records, best similarity and your ASVs found there, clicking
lists the reference ASVs, the queries matching them and the occurrence
metadata. **Download** saves the results per ASV, all matches or all records
as TSV.

The interface has one dark theme (black background, aqua accents), whatever
the system setting. The tab bar links to the group behind the application,
"Powered by dnambientale.it" (https://dnambientale.it/en/index.html, new
tab), with its icon `src/assets/dnambientale.png`, a copy of
`immagini/loghi/dnambientale-header.png` of the dnambientale site.

The default base map is a free world map bundled with the package:
[Natural Earth](https://www.naturalearthdata.com/) 1:50m countries (public
domain), drawn with [Leaflet](https://leafletjs.com/) 1.9.4 (BSD-2-Clause,
licence in `src/assets/leaflet/LICENSE`), also bundled. The application
therefore works **without internet access**. The layer button of the map
switches to OpenStreetMap tiles (free, online, more detailed when zoomed
in). `src/assets/world.geojson` is rebuilt from the Natural Earth file with
`tools/make_world_map.jl` (country names only, coordinates rounded to
3 decimals).

`serve` listens on `127.0.0.1` by default; `host = "0.0.0.0"` serves the
network, without authentication. Genie keeps one set of routes per Julia
session, so one reference dataset is served at a time.

## Reference data

A reference is the ASVs of **one marker and one primer pair**: name them with
`load_reference(...; marker = "COI-5P", primers = "BF3-BR2")` (or `--marker`
and `--primers`) and the page shows them next to the name of the reference,
in the top bar and in the first sentence of the search. ASVs of another
marker or primer pair belong to another reference.

The tab **Reference data** documents the reference: menus for the marker and
the reference database (one entry each for now), a short description, the
studies (country, area, sample type, primers, link to the raw data, number of
ASVs and records, and where the coordinates come from) and a general account
of the processing. It is read from a JSON file given with
`load_reference(...; about = "reference_about.json")` (or `--about`; keys in
the docstring of `load_reference`); the counts per study are computed from
the occurrence column `dataset`. Records of studies whose coordinates are
approximate (random points in the study area, or one point for the whole
area) are marked as such in the popups of both maps. Everything in the file
is public: write a description to document the reference, not the scripts.

- **FASTA** with the reference ASVs (unique identifiers). Taxonomy may be
  given in SINTAX headers: `>REF_1;tax=k:Animalia,p:Arthropoda,...,s:Genus_species;`
  (`_` becomes a space).
- **Metadata table** (optional, TSV, or CSV for `.csv`; quoted fields
  allowed), **one row per occurrence**:

  | id | kingdom | ... | species | latitude | longitude | site | date |
  |---|---|---|---|---|---|---|---|
  | REF_1 | Animalia | ... | Genus species | 46.35 | 10.70 | Site A | 2025-08-26 |
  | REF_1 | Animalia | ... | Genus species | 46.06 | 11.38 | Site B | 2025-07-01 |

  - identifier: `id`, `asv`, `asv_id`, `reference`, `reference_id`,
    `sequence_id` (or `id_column`);
  - coordinates in decimal degrees (WGS84): `latitude`/`lat`,
    `longitude`/`lon`/`lng`, Darwin Core `decimalLatitude`/`decimalLongitude`
    (or `latitude_column`, `longitude_column`). Rows with both empty only add
    taxonomy;
  - taxonomy: columns named after the ranks in `KNOWN_RANKS` (domain ...
    species, subspecies) and/or `taxonomy`/`lineage`. They must agree among
    the rows of a reference; metadata values replace header values rank by
    rank;
  - any other column (site, date, habitat, sample, ...) is shown in the map
    popups.

  **Everything in the occurrences is public.** Any visitor can search at low
  identity and read the occurrence columns and coordinates of most
  references. Remove personal or internal data (collectors, e-mail addresses,
  notes) from the table, or publish only some columns with `fields`
  (`--fields site,habitat,date`; `none` for no columns). Round the
  coordinates of sensitive sites (e.g. protected species) with
  `coordinate_digits` (`--coordinate-digits 2` ≈ 1 km, `1` ≈ 11 km). The
  server logs the published columns when it starts.

  Rows of identifiers absent from the FASTA file are skipped with a warning;
  invalid or out-of-range coordinates are errors.

## Haplotype networks

The second tab of the page, **Haplotype network**, builds the haplotype
network of the reference ASVs of one taxon and shows where each haplotype was
found:

1. networks exist only for the **species matched by your ASVs** in the last
   search: the tab lists them (none before a search); the button under the
   matches of an ASV opens the network of the species of its best match;
2. your ASVs matching the species are placed in its network (an ASV
   identical to a reference falls in its haplotype; nodes containing them
   have an aqua ring); the server checks again that they match the species,
   at the similarity of the search, and refuses any other species;
3. choose the method: median-joining (MJN, default), statistical parsimony
   (TCS), minimum spanning network (MSN) or tree (MST), and epsilon for
   MJN/MSN.

The sequences are first put in the same orientation (8-mers shared on
either strand; e.g. an ASV of yours read on the other strand is turned
around, and the page lists the sequences reversed); a sequence sharing no
8-mer with the others cannot be compared and is left out, with a note. They
are aligned with ProfileAlignments.jl and the network is computed by
HaplotypeNetworks.jl. As in PopART, columns with a gap or an
ambiguous base in any sequence are left out: ASVs of different lengths are
compared over their common region and indels do not count. Haplotypes are
numbered by decreasing total number of occurrences, after identical
sequences are merged (`H1` is the most widespread);
node area grows with the occurrences, inferred nodes (median vectors,
intermediates) are small grey dots and edges with more than one mutation show
their number.

The tab shows the map, with the network in a small translucent panel at its
top right (**Group** frames the group, **Fit** the whole network, ⤢ enlarges
the panel; wheel to zoom, drag to move; circles keep their size on screen,
so zooming in spreads them apart). Under **Grouping**:

- **Centre**: by default the haplotype of your ASV (the one the tab was
  opened for, else the first of your ASVs, else `H1`); any haplotype can be
  chosen, also by clicking a node or a row of the table;
- **Mutations from the centre** (slider): the group holds the haplotypes at
  most that many differences away from the centre (0: the centre only; 1:
  also those one mutation away; ...). Distances are numbers of differences
  between the haplotypes over the sites used, whatever the network method.

Places and nodes have two colours: orange, with a glow, for the group and
blue for the other haplotypes; each place is a pie chart of its occurrences
in and out of the group. The pair was checked for the dark map (any pair:
ΔE 31.8 normal vision, 26.8 colour-blind, at least 3:1 against land and sea).
The haplotype table lists the distance from the centre, sorted by it, and
can be downloaded as TSV with the group of each haplotype.

In Julia:

```julia
tn = taxon_network(ref, "species", "Fictiferus alpinus"; method = :mjn,
                   queries = ["ASV_1" => "ACGT..."])
tn.network                     # HaplotypeNetworks.HaplotypeNetwork
network_sites(ref, tn)         # places with haplotype counts
taxon_list(ref)                # taxa of each rank with their number of references

import HaplotypeNetworks as HN
HN.write_svg("network.svg", tn.network)   # any HaplotypeNetworks function
```

At most `max_network_sequences` (300) sequences per network
(`--max-network-sequences`): median-joining grows quickly with the number of
haplotypes. Networks run one at a time, like searches.

## Running on a server

`bin/asvatlas.jl` starts the application from the command line; every
option can also be given as an environment variable (`--help` lists them):

```
julia -t auto bin/asvatlas.jl --fasta reference.fasta --metadata reference_metadata.tsv       --name "Reference ASVs" --marker COI-5P --primers BF3-BR2 --host 127.0.0.1 --port 8080
```

At start-up the reference is loaded and indexed and one search is run, so
that the first user does not wait for compilation (about 4 s on the
example; the following searches of the example take 0.02 s). Each search is
logged (`[ Info: Search: 6 queries at ≥ 97% identity, 4 with hits, 0.02 s`)
and written at once, also when the output is redirected to a file.
`GET /health` answers `{"status": "ok", "references": ...}` for monitoring.

The recommended set-up on a Linux server is the application on
`127.0.0.1` behind nginx, which adds HTTPS, an optional password and a limit
on the number of searches per client. The `deploy/` folder has examples to
adapt:

- `asvatlas.service`: systemd unit (dedicated user, read-only data in
  `/srv/asvatlas`, restart on failure, log in `journalctl -u asvatlas`).
  Copy `ASVAtlas.jl`, `NucleotideSearch.jl`, `HaplotypeNetworks.jl` and
  `ProfileAlignments.jl` side by side in
  `/opt/asvatlas`, then once, as the `asvatlas` user and with the unit's
  `JULIA_DEPOT_PATH`:
  `julia --project=/opt/asvatlas/ASVAtlas.jl -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`
  (delete a `Manifest.toml` copied from Windows first).
- `nginx.conf`: the site, with HTTPS via certbot, `auth_basic` (commented
  out), 30 searches per minute per address and the variant for a sub-path
  (`https://example.org/asvatlas/`): the page uses relative URLs, so it works
  under any path.
- `Dockerfile` and `docker-compose.yml`: a container image built from the
  folder that contains the four packages
  (`docker build -f ASVAtlas.jl/deploy/Dockerfile -t asvatlas .`), with the
  reference data mounted read-only in `/data`.

There is no authentication in the application itself: on a server reachable
from the internet, use the nginx password or a private network. Searches run
one at a time, so the limits below also bound how long one user can hold the
server.

### Container image

The workflow `.github/workflows/docker.yml` builds `deploy/Dockerfile` on
GitHub and publishes the image as `ghcr.io/alexology/asvatlas`: tag `latest`
(and `main`, `sha-…`) for every push to `main`, `1.2.3` and `1.2` for a
release tag `v1.2.3`; pull requests only build it. The workflow checks out the
three unregistered dependencies (public repositories) next to ASVAtlas.jl, as
the Dockerfile expects. The image holds the code only, compiled for generic
x86-64 CPUs: the reference files are mounted at run time (see
`deploy/docker-compose.yml`, which uses the published image). On first
publication the package on GitHub is private: make it public in its settings,
or log the server in to `ghcr.io`.

### Keeping the reference private

The page answers questions about the ASVs it is sent; it does not publish the
reference:

- no reference sequence is ever sent (matches give identifiers, similarity,
  coordinates on the reference and taxonomy of the reference ASVs found);
- there is no list of taxa or of reference ASVs (`GET /api/dataset` gives
  only the name and counts);
- matches need at least `min_identity` percent similarity (90 by default,
  `--min-identity`), so unrelated sequences reveal nothing;
- networks are built only for species matched by the ASVs sent.

Keep the reference files (FASTA and metadata) **outside the repository**: the
application reads them at start from any folder (`--fasta`, `--metadata`, or
the arguments of `load_reference`); only the synthetic example is part of the
package. What the page shows of each occurrence is chosen with `fields` and
`coordinate_digits`.

### Security measures

- Requests: bodies over `max_body` (16 MiB) are refused before being read
  (by `Content-Length`, or while reading a chunked body); at most
  `max_connections` (128) open connections; per search at most
  `--max-queries` sequences of `--max-length` bases (20,000), identity from
  50%, and at most 1,000,000 possible hit rows.
- `--allowed-hosts` / `ASVATLAS_ALLOWED_HOSTS` refuses requests for other
  host names (DNS rebinding); set it to the public name on a server.
- `--fields` / `ASVATLAS_FIELDS` and `--coordinate-digits` limit what the
  occurrences publish (see "Reference data").
- `--example` / `ASVATLAS_EXAMPLE` (`serve(...; example)`) is a FASTA file
  offered by "Try an example" and sent to every visitor: use sequences that
  may be public, not private reference sequences.
- `POST /api/search` accepts only `Content-Type: application/json`, which a
  cross-site page cannot send without a CORS preflight; no origin is allowed.
- Only a fixed list of files is served (no file system paths); Genie's static
  file handler is disabled.
- Every response carries a strict Content Security Policy (scripts, styles
  and data only from the server; OpenStreetMap tiles as the only external
  images), `X-Content-Type-Options: nosniff`, `X-Frame-Options: DENY`,
  `Referrer-Policy: strict-origin` (only the origin, which the OpenStreetMap
  tile servers require) and same-origin COOP/CORP. Error pages are a
  fixed JSON message (Genie's own pages echo the requested path).
- The page writes all data with `textContent`, never as HTML. Exported TSV
  cells starting with `=`, `+`, `-` or `@` get a leading apostrophe, so
  spreadsheet programs do not run them as formulas (the Julia `write_table`
  writes values unchanged).
- Internal errors are logged on the server and answered with a generic
  message.

Not covered by the application, to be handled by the reverse proxy: HTTPS,
authentication, rate limiting and slow-client timeouts (HTTP.jl's read
timeout would also cut long searches). For large references the index
takes memory (see NucleotideSearch: about 3.3 GB per million 650-bp
references) and is rebuilt at every start.

## Julia API

```julia
queries = read_queries("asvs.fasta")              # or parse_queries(text), or a vector of sequences
res = search_asvs(ref, queries; identity = 97, query_cov = 0.9, maxhits = 50)

res[1].hits[1].identity, ref.ids[res[1].hits[1].reference]
best_taxa(ref, res[1])                            # names of the best hits

write_table("summary.tsv", summary_table(ref, res))       # one row per query
write_table("hits.tsv", hit_table(ref, res))              # one row per hit, with lineage
write_table("occurrences.tsv", occurrence_table(ref, res)) # one row per hit occurrence
```

The tables are `NamedTuple`s of vectors, so `DataFrame(hit_table(ref, res))`
works with DataFrames.jl.

Search options (`search_asvs` and the web page):

- `identity` (percent): minimum identity of a hit. Identity is matches over
  aligned columns, end gaps excluded (`iddef = 2` of NucleotideSearch,
  VSEARCH's default). 100 keeps only references identical over the aligned
  region.
- `query_cov`: minimum fraction of the query bases aligned to reference
  bases (0.9 by default; 1 requires the whole query).
- `mode = :query` (default) aligns the whole query inside the reference, so
  references longer than the queries (e.g. complete barcodes) match; a
  query end is soft-clipped only where it diverges strongly. `:global`
  aligns whole sequences.
- `maxhits` (50) and `maxrejects` (64): candidates, ranked by shared 8-mers,
  are aligned until `maxhits` are accepted or `maxrejects` rejected. With
  many near-identical references above the threshold, raise `maxhits` to see
  them all. `ASVResult.exhaustive` tells whether every reference was
  examined.
- `strand = :both` (default) also searches the reverse complement.

`serve` options: `host` (an IP address or host name), `port`, `async`
(returns the server; stop it with `shutdown`), `max_queries` (10,000),
`max_hits` (500), `max_length` (20,000 bases), `max_body` (16 MiB),
`max_connections` (128), `allowed_hosts` (any), `threads`, `warmup`. The API
used by the page is `GET /api/dataset` and `POST /api/search` with
`Content-Type: application/json` and `{"sequences": "<FASTA>", "identity": 97,
"query_cov": 0.9, "max_hits": 50, "strand": "both"}` (identity from 50).

## Validation

`julia -t 4 --project=. test/runtests.jl` (409 tests): delimited files with
quotes and CRLF, SINTAX headers, metadata precedence and errors, query
parsing, search at 100/97/90% (nested hit sets, thresholds, reverse
complement, partial queries with their coordinates on the reference, ties
with different taxa), tables and TSV round trip, and the web server
(pages, bundled Leaflet and world map, API results equal to `search_asvs`,
400/404/405/413/415/421 errors, path traversal, Host check, oversized and
chunked bodies on a raw socket, security headers, `/health`, command-line
options and environment variables), and the haplotype networks (the four
methods, numbering by occurrences, queries joining identical references, a
trimmed query aligned to its reference, a single sequence, place counts,
limits and the `/api/network` endpoint, restricted to species matched by the ASVs sent; no `/api/taxa`). The page was also driven in
headless Chrome (desktop and phone width, all external requests
blocked to prove offline use, OpenStreetMap layer; network tab: building,
selection shared by network, map and table, opening from the search, zoom,
export), and through a reverse
proxy under a sub-path with `bin/asvatlas.jl` configured by environment
variables on `0.0.0.0`. The systemd, nginx and Docker files were not run
(no Linux server, nginx or Docker available when they were written).

**Comparison with R and BLAST** ([docs/VALIDATION.md](docs/VALIDATION.md)),
on a real reference of 8,198 COI ASVs: for 215 queries ASVAtlas found the
same reference ASVs as an exhaustive alignment in R (pwalign) at every
similarity level from 100 to 90%, with identical identity, coverage and
strand, also with the early stop of the web page (7 s against 271 s
exhaustive and 340 s for BLAST). For 11 species, haplotypes, distances, MST,
MSN (union of all minimum spanning trees) and TCS components agree with ape,
pegas and haplotypes; the median-joining networks differ from `pegas::mjn`
and keep the distances between observed haplotypes better. The scripts are in
`tools/validation/`.

Search speed and memory on larger references are those of NucleotideSearch
(see its README: about 3.3 GB of index per million 650-bp references).
