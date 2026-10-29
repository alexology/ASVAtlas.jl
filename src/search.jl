"""
    IDENTITY_LEVELS

Percent identities offered by the web application (any value in (0, 100]
can be used with [`search_asvs`](@ref)).
"""
const IDENTITY_LEVELS = (100, 99, 98, 97, 95, 90, 85, 80)

"""
    ReferenceHit

One reference ASV matching a query: `reference` (index in the
[`ReferenceDataset`](@ref)), `identity` (percent, end gaps excluded),
`coverage` (fraction of the query bases aligned to reference bases),
`strand` (`'-'` if the reverse complement of the query matched), the
1-based inclusive coordinates `query_start`, `query_end`, `target_start`,
`target_end` and the `cigar` string of the alignment.
"""
struct ReferenceHit
    reference::Int
    identity::Float64
    coverage::Float64
    strand::Char
    query_start::Int
    query_end::Int
    target_start::Int
    target_end::Int
    cigar::String
end

"""
    ASVResult

Search result of one query: `id`, `sequence`, `hits` (by decreasing
identity) and `exhaustive`, true when every reference was examined, so that
no reference missing from `hits` can match.
"""
struct ASVResult
    id::String
    sequence::String
    hits::Vector{ReferenceHit}
    exhaustive::Bool
end

"""
    parse_queries(text) -> Vector{Pair{String,String}}

Query ASVs from FASTA text, or from plain text with one sequence per line
(named `query_1`, `query_2`, ...). Sequences are converted to upper case,
`U` to `T`, with whitespace and gaps removed. Identifiers must be unique.
"""
function parse_queries(text::AbstractString)
    lines = [rstrip(l) for l in split(text, '\n')]
    filter!(!isempty, lines)
    isempty(lines) && throw(ArgumentError("no query sequences"))
    queries = Pair{String,String}[]
    if startswith(lstrip(first(lines)), '>')
        io = IOBuffer(join(lines, '\n'))
        FASTA.Reader(io) do reader
            for rec in reader
                id = String(identifier(rec))
                isempty(id) && throw(ArgumentError("a query has an empty identifier"))
                push!(queries, id => _normalize_sequence(sequence(String, rec), id))
            end
        end
    else
        any(l -> startswith(lstrip(l), '>'), lines) &&
            throw(ArgumentError("FASTA text must start with a header line (\">id\")"))
        for (i, l) in enumerate(lines)
            id = "query_$i"
            push!(queries, id => _normalize_sequence(l, id))
        end
    end
    ids = first.(queries)
    allunique(ids) || throw(ArgumentError("query identifiers must be unique"))
    return queries
end

"""
    read_queries(path) -> Vector{Pair{String,String}}

[`parse_queries`](@ref) on the contents of a file.
"""
read_queries(path::AbstractString) = parse_queries(read(path, String))

_queries(q::AbstractVector{<:Pair}) =
    [String(first(p)) => _normalize_sequence(last(p), first(p)) for p in q]
_queries(q::AbstractVector{<:AbstractString}) =
    ["query_$i" => _normalize_sequence(s, "query_$i") for (i, s) in enumerate(q)]

"""
    search_asvs(ref, queries; identity = 97, query_cov = 0.9, maxhits = 50,
                maxrejects = 64, mode = :query, strand = :both,
                threads = Threads.nthreads()) -> Vector{ASVResult}

Compare query ASVs with the reference ASVs of `ref` and keep the references
with percent identity at least `identity` (e.g. 100 for identical sequences,
97, 95, ...) and at least `query_cov` of the query aligned. `queries` is a
vector of `id => sequence` pairs (see [`parse_queries`](@ref)) or of
sequences.

Up to `maxhits` hits are returned per query, by decreasing identity.
Candidate references are examined in order of shared k-mers until `maxhits`
are accepted or `maxrejects` rejected (0: no limit). In `:query` mode
(default) the whole query is aligned within the reference, so references
longer than the query (e.g. full barcodes) can match; `:global` aligns
whole sequences. Both strands are searched unless `strand = :plus`.

The taxonomy reported for a hit is the one of the reference ASV: no new
assignment is computed. See [`best_taxa`](@ref) for the names of the best
hits and [`summary_table`](@ref), [`hit_table`](@ref) and
[`occurrence_table`](@ref) for tables.
"""
function search_asvs(ref::ReferenceDataset, queries::AbstractVector; identity::Real = 97,
                     query_cov::Real = 0.9, maxhits::Integer = 50, maxrejects::Integer = 64,
                     mode::Symbol = :query, strand::Symbol = :both,
                     threads::Integer = Threads.nthreads())
    (isfinite(identity) && 0 < identity <= 100) ||
        throw(ArgumentError("identity must be a percentage in (0, 100], got $identity"))
    0 < query_cov <= 1 || throw(ArgumentError("query_cov must be in (0, 1], got $query_cov"))
    maxhits >= 1 || throw(ArgumentError("maxhits must be at least 1"))
    maxrejects >= 0 || throw(ArgumentError("maxrejects must not be negative"))
    qs = _queries(queries)
    isempty(qs) && return ASVResult[]
    allunique(first.(qs)) || throw(ArgumentError("query identifiers must be unique"))
    reports = Vector{SearchReport}(undef, length(qs))
    found = search(last.(qs), ref.db; mode, strand, min_identity = identity, query_cov,
                   maxhits, maxaccepts = maxhits, maxrejects, cigar = true, reports, threads)
    return [ASVResult(id, seq,
                      [ReferenceHit(h.target, h.identity, (h.matches + h.mismatches) / length(seq),
                                    h.strand, h.query_start, h.query_end, h.target_start,
                                    h.target_end, h.cigar) for h in found[q]],
                      reports[q].exhaustive)
            for (q, (id, seq)) in enumerate(qs)]
end

"""
    best_taxa(ref, result) -> Vector{String}

Distinct deepest known names of the references hit with the highest
identity, in hit order (empty without hits or taxonomy). Several names mean
that equally similar references carry different taxonomies.
"""
function best_taxa(ref::ReferenceDataset, r::ASVResult)
    isempty(r.hits) && return String[]
    best = first(r.hits).identity
    names = String[]
    for h in r.hits
        h.identity == best || continue
        t = taxon(ref, h.reference)
        isempty(t) || t in names || push!(names, t)
    end
    return names
end
