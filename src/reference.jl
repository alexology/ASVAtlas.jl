"""
    KNOWN_RANKS

Metadata columns read as taxonomic ranks (case-insensitive), from the most
general to the most specific. A column named `taxonomy` or `lineage` is read
as a free-text lineage after them.
"""
const KNOWN_RANKS = ("domain", "superkingdom", "kingdom", "phylum", "division", "class",
                     "subclass", "order", "suborder", "family", "subfamily", "tribe", "genus",
                     "subgenus", "species", "subspecies")
const _LINEAGE_COLUMNS = ("taxonomy", "lineage")
const _SINTAX_RANKS = Dict('d' => "domain", 'k' => "kingdom", 'p' => "phylum", 'c' => "class",
                           'o' => "order", 'f' => "family", 'g' => "genus", 's' => "species")
const _ID_COLUMNS = ("id", "asv", "asv_id", "reference", "reference_id", "sequence_id")
const _LATITUDE_COLUMNS = ("latitude", "lat", "decimallatitude")
const _LONGITUDE_COLUMNS = ("longitude", "lon", "long", "lng", "decimallongitude")
const _MISSING_VALUES = Set(["", "na", "nan", "null", "none", "missing"])
const _IUPAC = Set("ACGTURYSWKMBDHVN")

"""
    Occurrence

A place where a reference ASV was found: `latitude` and `longitude` in
decimal degrees (WGS84) and `fields`, the other non-empty metadata of that
row (column name => value, in file order).
"""
struct Occurrence
    latitude::Float64
    longitude::Float64
    fields::Vector{Pair{String,String}}
end

"""
    ReferenceDataset

Reference ASVs with their taxonomy and occurrences, built by
[`load_reference`](@ref). Fields: `name`, `marker` and `primers` (as given to
`load_reference`, `""` when unknown), `about` (the description of the
reference, JSON text, `""` when none), `db` (the NucleotideSearch
`SequenceDB`), `ids`, `lengths`, `ranks` (taxonomy columns, general to
specific), `lineages` (one vector of names per reference, aligned with
`ranks`, `""` when unknown) and `occurrences` (one vector per reference).
"""
struct ReferenceDataset
    name::String
    marker::String
    primers::String
    about::String
    db::SequenceDB
    ids::Vector{String}
    lengths::Vector{Int}
    ranks::Vector{String}
    lineages::Vector{Vector{String}}
    occurrences::Vector{Vector{Occurrence}}
end

Base.length(ref::ReferenceDataset) = length(ref.ids)

function Base.show(io::IO, ref::ReferenceDataset)
    n_occ = sum(length, ref.occurrences; init = 0)
    n_geo = count(!isempty, ref.occurrences)
    print(io, "ReferenceDataset(\"", ref.name, "\", ",
          isempty(ref.marker) ? "" : ref.marker * ", ", isempty(ref.primers) ? "" : ref.primers * ", ",
          length(ref), " references, ", n_geo,
          " georeferenced, ", n_occ, " occurrences, ranks: ",
          isempty(ref.ranks) ? "none" : join(ref.ranks, ", "), ")")
end

"""
    taxon(ref, i) -> String

Deepest known name in the lineage of reference `i` (`""` if none).
"""
function taxon(ref::ReferenceDataset, i::Integer)
    for name in Iterators.reverse(ref.lineages[i])
        isempty(name) || return name
    end
    return ""
end

# Upper case, `U` to `T`, without whitespace and alignment gaps; only IUPAC
# nucleotide codes are accepted.
function _normalize_sequence(seq::AbstractString, id::AbstractString)
    io = IOBuffer()
    for c in seq
        (isspace(c) || c == '-' || c == '.') && continue
        u = uppercase(c)
        u in _IUPAC || throw(ArgumentError("invalid nucleotide '$c' in sequence $id"))
        print(io, u == 'U' ? 'T' : u)
    end
    s = String(take!(io))
    isempty(s) && throw(ArgumentError("sequence $id is empty"))
    return s
end

_value(s::AbstractString) = (v = String(strip(s)); lowercase(v) in _MISSING_VALUES ? "" : v)

# Identifier and SINTAX taxonomy (`id;tax=k:Animalia,p:Arthropoda,...;`) of a
# FASTA header. Names in SINTAX headers use `_` for spaces.
function _parse_header(identifier::AbstractString, desc::AbstractString)
    m = match(r";tax=([^;]*)", desc)
    m === nothing && return String(identifier), Dict{String,String}()
    id = String(strip(first(split(desc, ';'))))
    tax = Dict{String,String}()
    for item in split(m.captures[1], ',')
        p = split(item, ':'; limit = 2)
        length(p) == 2 && length(p[1]) == 1 || continue
        rank = get(_SINTAX_RANKS, only(p[1]), nothing)
        rank === nothing && continue
        v = _value(replace(p[2], '_' => ' '))
        isempty(v) || (tax[rank] = v)
    end
    return id, tax
end

function _read_reference_fasta(path::AbstractString)
    ids = String[]
    seqs = String[]
    taxa = Dict{String,String}[]
    FASTA.Reader(open(path)) do reader
        for rec in reader
            id, tax = _parse_header(identifier(rec), description(rec))
            isempty(id) && throw(ArgumentError("$path has a record without identifier"))
            push!(ids, id)
            push!(seqs, _normalize_sequence(sequence(String, rec), id))
            push!(taxa, tax)
        end
    end
    isempty(ids) && throw(ArgumentError("$path contains no sequences"))
    if !allunique(ids)
        seen = Set{String}()
        dup = first(id for id in ids if !(id in seen || (push!(seen, id); false)))
        throw(ArgumentError("$path: identifier $dup appears more than once"))
    end
    return ids, seqs, taxa
end

function _find_column(header, given, candidates, what)
    lower = lowercase.(header)
    if given !== nothing
        j = findfirst(==(lowercase(String(given))), lower)
        j === nothing && throw(ArgumentError("metadata has no column \"$given\" ($what)"))
        return j
    end
    for c in candidates
        j = findfirst(==(c), lower)
        j === nothing || return j
    end
    return nothing
end

function _coordinate(v, what, limit, id, row)
    x = tryparse(Float64, v)
    (x === nothing || !isfinite(x)) &&
        throw(ArgumentError("metadata row $row ($id): invalid $what \"$v\""))
    abs(x) <= limit ||
        throw(ArgumentError("metadata row $row ($id): $what $x outside [-$limit, $limit]"))
    return x
end

"""
    load_reference(fasta; metadata = nothing, delim = nothing, id_column = nothing,
                   latitude_column = nothing, longitude_column = nothing,
                   fields = nothing, coordinate_digits = nothing,
                   name = basename(fasta), marker = "", primers = "", about = nothing,
                   k = 8) -> ReferenceDataset

Read the reference ASVs from `fasta` and, optionally, their taxonomy and
occurrences from the delimited file `metadata` (tab-separated, or
comma-separated for `.csv`, unless `delim` is given).

The metadata table has one row per occurrence: the reference identifier
(column `id`, `asv`, `asv_id`, `reference`, `reference_id` or `sequence_id`
unless `id_column` is given), the coordinates in decimal degrees (`latitude`
or `lat`, `longitude`, `lon` or `lng`, also Darwin Core
`decimalLatitude`/`decimalLongitude`), the taxonomy (columns named after
[`KNOWN_RANKS`](@ref) and/or `taxonomy`) and any other column (site, date,
sample, ...), shown with the occurrence. The occurrence columns are public in
the web application: `fields` (column names, case-insensitive; empty for
none) keeps only those, `nothing` keeps them all, and `coordinate_digits`
rounds the coordinates to that many decimals (e.g. 2 ≈ 1 km, 1 ≈ 11 km), for
sites that should not be published exactly. A reference may have several rows;
its taxonomy must agree among them. Rows with both coordinates empty only add
taxonomy; references without rows have no occurrences. Rows of identifiers
not in `fasta` are skipped with a warning.

The taxonomy can also come from SINTAX headers
(`>ASV1;tax=k:Animalia,p:Arthropoda,...,s:Genus_species;`); metadata values
take precedence rank by rank. `k` is the k-mer length of the search index.

`marker` (e.g. `"COI-5P"`) and `primers` (e.g. `"BF3-BR2"`) describe the
amplicon of the reference ASVs; the web page shows them next to its name. A
reference of another marker or primer pair is another dataset.

`about` is a JSON file describing the reference, shown by the tab "Reference
data" of the web page (at most 256 KiB): an object with the optional keys
`description` (text), `processing` (list of texts), `coordinates` (text) and
`sources`, a list of the studies, each an object with `dataset` (the value of
the occurrence column `dataset`), `label`, `citation`, `country`, `area`,
`sample`, `primers`, `accession`, `url` (http or https) and `coordinates`
(`"sites"`, `"approximate"`, `"study area"` or `"none"`: where the
coordinates of its records come from). All of it is public.
"""
function load_reference(fasta::AbstractString; metadata::Union{Nothing,AbstractString} = nothing,
                        delim::Union{Nothing,Char} = nothing, id_column = nothing,
                        latitude_column = nothing, longitude_column = nothing,
                        fields = nothing, coordinate_digits::Union{Nothing,Integer} = nothing,
                        name::AbstractString = basename(fasta), marker::AbstractString = "",
                        primers::AbstractString = "", about::Union{Nothing,AbstractString} = nothing,
                        k::Integer = 8)
    about_text = about === nothing ? "" : _read_about(about)
    coordinate_digits === nothing || 0 <= coordinate_digits <= 8 ||
        throw(ArgumentError("coordinate_digits must be between 0 and 8"))
    fields isa AbstractString && throw(ArgumentError("fields must be a collection of column names"))
    ids, seqs, taxa = _read_reference_fasta(fasta)
    index = Dict(id => i for (i, id) in enumerate(ids))
    occurrences = [Occurrence[] for _ in ids]
    meta_taxa = [Dict{String,String}() for _ in ids]
    if metadata !== nothing
        _read_metadata!(occurrences, meta_taxa, index, metadata; delim, id_column,
                        latitude_column, longitude_column, fields, coordinate_digits)
    elseif fields !== nothing && !isempty(fields)
        throw(ArgumentError("fields needs a metadata table"))
    end
    present = Set{String}()
    for t in taxa, r in keys(t)
        push!(present, r)
    end
    for t in meta_taxa, r in keys(t)
        push!(present, r)
    end
    ranks = [r for r in (KNOWN_RANKS..., _LINEAGE_COLUMNS...) if r in present]
    lineages = [[get(meta_taxa[i], r, get(taxa[i], r, "")) for r in ranks] for i in eachindex(ids)]
    db = SequenceDB(ids, seqs; k)
    return ReferenceDataset(String(name), String(strip(marker)), String(strip(primers)), about_text, db, ids,
                            length.(seqs), ranks, lineages, occurrences)
end

function _read_metadata!(occurrences, meta_taxa, index, path; delim, id_column,
                         latitude_column, longitude_column, fields = nothing,
                         coordinate_digits = nothing)
    header, rows = read_delimited(path; delim)
    jid = _find_column(header, id_column, _ID_COLUMNS, "identifier")
    jid === nothing &&
        throw(ArgumentError("metadata has no identifier column; name it with `id_column`"))
    jlat = _find_column(header, latitude_column, _LATITUDE_COLUMNS, "latitude")
    jlon = _find_column(header, longitude_column, _LONGITUDE_COLUMNS, "longitude")
    (jlat === nothing) == (jlon === nothing) ||
        throw(ArgumentError("metadata has only one of the latitude and longitude columns"))
    lower = lowercase.(header)
    rank_cols = [(j, lower[j]) for j in eachindex(header)
                 if lower[j] in KNOWN_RANKS || lower[j] in _LINEAGE_COLUMNS]
    length(unique(last.(rank_cols))) == length(rank_cols) ||
        throw(ArgumentError("metadata repeats a taxonomy column"))
    used = Set([jid; first.(rank_cols)])
    jlat === nothing || push!(used, jlat, jlon)
    field_cols = [j for j in eachindex(header) if !(j in used)]
    if fields !== nothing
        wanted = Set(lowercase(strip(String(f))) for f in fields)
        missing_cols = setdiff(wanted, lower[field_cols])
        isempty(missing_cols) ||
            throw(ArgumentError("fields: metadata has no occurrence column " *
                                join(("\"$c\"" for c in sort!(collect(missing_cols))), ", ")))
        filter!(j -> lower[j] in wanted, field_cols)
    end
    _round(x) = coordinate_digits === nothing ? x : round(x; digits = coordinate_digits)
    unknown = Set{String}()
    for (row, r) in enumerate(rows)
        id = String(strip(r[jid]))
        isempty(id) && throw(ArgumentError("metadata row $row has an empty identifier"))
        i = get(index, id, 0)
        if i == 0
            push!(unknown, id)
            continue
        end
        for (j, rank) in rank_cols
            v = _value(r[j])
            isempty(v) && continue
            old = get(meta_taxa[i], rank, "")
            isempty(old) || old == v ||
                throw(ArgumentError("reference $id has two values for $(header[j]): \"$old\" and \"$v\""))
            meta_taxa[i][rank] = v
        end
        jlat === nothing && continue
        lat = _value(r[jlat])
        lon = _value(r[jlon])
        isempty(lat) && isempty(lon) && continue
        (isempty(lat) || isempty(lon)) &&
            throw(ArgumentError("metadata row $row ($id) has only one coordinate"))
        fields = [header[j] => String(strip(r[j])) for j in field_cols if !isempty(_value(r[j]))]
        values = [header[j] => String(strip(r[j])) for j in field_cols if !isempty(_value(r[j]))]
        push!(occurrences[i], Occurrence(_round(_coordinate(lat, "latitude", 90, id, row)),
                                         _round(_coordinate(lon, "longitude", 180, id, row)), values))
    end
    isempty(unknown) ||
        @warn "Skipped metadata rows of $(length(unknown)) identifier(s) not in the FASTA file, " *
              "e.g. $(first(sort!(collect(unknown))))"
    return nothing
end

# The description of a reference (see `load_reference`): a JSON object, read
# once and kept as text; errors name the file.
const _ABOUT_KINDS = ("sites", "approximate", "study area", "none")

function _read_about(path::AbstractString)
    filesize(path) <= 2^18 || throw(ArgumentError("$path: the description must be at most 256 KiB"))
    text = read(path, String)
    obj = try
        JSON3.read(text)
    catch
        throw(ArgumentError("$path: not valid JSON"))
    end
    obj isa JSON3.Object || throw(ArgumentError("$path: the description must be a JSON object"))
    for (k, src) in enumerate(get(obj, :sources, []))
        src isa JSON3.Object || throw(ArgumentError("$path: sources[$k] must be an object"))
        c = get(src, :coordinates, "sites")
        c in _ABOUT_KINDS ||
            throw(ArgumentError("$path: sources[$k].coordinates must be one of $(join(_ABOUT_KINDS, ", "))"))
        url = get(src, :url, "")
        (isempty(url) || occursin(r"^https?://", url)) ||
            throw(ArgumentError("$path: sources[$k].url must start with http:// or https://"))
    end
    return text
end
