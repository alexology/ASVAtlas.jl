# Changelog

## 0.1.0 (2026-10-25)

- First version: reference ASVs from FASTA (optionally SINTAX headers) and an
  occurrence/taxonomy table; search at a chosen identity with
  NucleotideSearch; summary, hit and occurrence tables; web application on
  Genie.jl with a bundled offline world map (Natural Earth, public domain) and
  Leaflet, OpenStreetMap as optional online layer.
- Server use: `bin/asvatlas.jl` (`ASVAtlas.main`, options or `ASVATLAS_*`
  environment variables), warm-up search, one log line per search written at
  once, `GET /health`, relative URLs for reverse proxies under a sub-path,
  examples for systemd, nginx and Docker in `deploy/`.
- Interface redesigned: one dark theme (black, aqua accents), header with the
  dataset figures, identity as buttons, result figures, identity marks in the
  table, marker tooltips, dark map and popups, FASTA drag and drop.
- Security review: request bodies limited before reading, connection limit,
  optional Host check (`allowed_hosts`), JSON-only search API, files served
  from a fixed list, security headers and fixed error messages on every
  response, query length and identity limits, TSV formula neutralization.
- Privacy: `load_reference(...; fields, coordinate_digits)` (CLI
  `--fields`, `--coordinate-digits`) choose the occurrence columns published
  and round the coordinates; the server logs the published columns at start.
- Haplotype network tab: networks (MJN, TCS, MSN, MST; HaplotypeNetworks.jl)
  of the reference ASVs of a taxon after alignment (ProfileAlignments.jl),
  optionally with the user's ASVs; pie charts of the haplotypes at each place
  on the map; `taxon_network`, `network_sites`, `taxon_list`,
  `taxon_references`, `reference_sequence`; API `GET /api/taxa`,
  `POST /api/network`.
- Networks (review of 2026-10-26): sequences are oriented on a common strand
  before the alignment, unrelated ones left out and reported (AT1);
  haplotypes numbered by total frequency after merging (AT3); Docker image and
  instructions include HaplotypeNetworks.jl and ProfileAlignments.jl (AT2).
- Network tab redesigned: map with the network in a small translucent panel;
  slider grouping the haplotypes within a number of mutations of a centre (by
  default the haplotype of the user's ASV), places and nodes in two colours
  (group / others); distances between haplotypes in `POST /api/network`;
  circles of constant size on screen; scrolling haplotype table.
- "Powered by dnambientale.it" link with the group's icon in the tab bar.
- First tab reorganised around "Where in the world is my ASV?": similarity
  100% by default (identical records), 99–90% to widen the search (changing it
  repeats the search), summary above the map, advanced options folded.
- Reference kept private: networks only for species matched by the ASVs sent
  (checked by the server), no `GET /api/taxa`, `min_identity` (90% by
  default) for every search from the page.
- OpenStreetMap tiles were shown as "Access blocked": the page sent no
  Referer (`Referrer-Policy: no-referrer`), which the tile servers require.
  Now `strict-origin` (origin only, never the path).
- First tab redesigned: a single form before the search ("Try an example"
  when the server has one: `serve(...; example)`, `--example`); after it the
  map first, with the similarity above it (immediate) and a one-sentence
  answer naming the countries; one line per ASV (species, places, countries,
  records, best similarity, network button, "Try 95%" when nothing is found);
  matches grouped by taxon, alignment details folded; map never closer than
  a region and refitted when the window changes size; country in tooltips
  and popups; a single Download menu; the key-figure boxes are gone.
- While a new similarity is searched, "Searching at …" is aqua with a dot and
  pulses slowly (still for users who ask for reduced motion), and the map is
  dimmed until the answer arrives.
- Validation against R and BLAST on a real reference (docs/VALIDATION.md),
  with the scripts to repeat it (tools/validation/).
- Optional sampling locations of the query ASVs: a CSV/TSV table (id,
  latitude, longitude in decimal degrees WGS84) read in the page and never
  sent to the server; gold stars on both maps, distance to the nearest
  matching record, dashed lines to the records of the selected ASV, columns
  in the results per ASV. Files edited on Windows were back to LF endings.
- An ASV can have several sampling places (several rows of the CSV, whose
  header is exactly asv_id, latitude, longitude); each sequence must appear
  once in the FASTA, otherwise the search is refused.
- Network tab: the stars of the sampling places are filled with the colour of
  the group (split when the ASVs of a place differ), with a gold border.
- Marker and primer pair of the reference (`load_reference(...; marker,
  primers)`, `--marker`, `--primers`), shown in the top bar and in the first
  sentence of the search; COI-5P / BF3-BR2 for freshwater_invertebrates.
- Tab "Reference data": marker and database menus, description, studies
  (with counts computed from the column `dataset`), processing and origin of
  the coordinates, from `load_reference(...; about)` / `--about`; popups mark
  records with approximate coordinates. Primers of freshwater_invertebrates:
  BF3-BR2.
- "New search" above the map (empties ASVs, coordinates, results and network)
  and "Clear" in the form; the reference is shown as "Freshwater
  invertebrates" (display name only; files unchanged).
- With a table of coordinates, ASVs without a row in it are marked ("no
  sampling place in your table"), also in the tooltip of their network node.
- GitHub workflow building the container image and publishing it on the
  GitHub Container Registry (ghcr.io/alexology/asvatlas); image built for
  generic x86-64 CPUs; docker-compose uses the published image.
- Usage statistics: usage log of counts only (`serve(...; usage_log)`,
  `--usage-log`, `usage_report`, `bin/asvatlas-usage.jl`; on in the image);
  `deploy/caddy/` with Caddy (HTTPS, access log with truncated IP addresses,
  reports at /stats/ behind a password) and GoAccess (`update-stats.sh`);
  privacy note in the footer.
