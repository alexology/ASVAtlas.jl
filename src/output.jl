# Column tables (NamedTuples of vectors, usable with Tables.jl, e.g.
# `DataFrame(hit_table(ref, results))`) and their tab-separated files.

_join_fields(fields) = join((string(k, "=", v) for (k, v) in fields), "; ")

"""
    summary_table(ref, results)

One row per query: `query`, `length`, `n_hits`, `best_identity`,
`best_references` (identifiers of the hits with the best identity, `;`
separated), `best_taxa` (their distinct names, see [`best_taxa`](@ref),
` | ` separated), `n_georeferenced` (hits with at least one occurrence),
and `n_occurrences` (occurrences of all hits).
"""
function summary_table(ref::ReferenceDataset, results::AbstractVector{ASVResult})
    best = [isempty(r.hits) ? missing : first(r.hits).identity for r in results]
    return (query = [r.id for r in results],
            length = [length(r.sequence) for r in results],
            n_hits = [length(r.hits) for r in results],
            best_identity = best,
            best_references = [join((ref.ids[h.reference] for h in r.hits
                                     if h.identity == first(r.hits).identity), ";")
                               for r in results],
            best_taxa = [join(best_taxa(ref, r), " | ") for r in results],
            n_georeferenced = [count(h -> !isempty(ref.occurrences[h.reference]), r.hits)
                               for r in results],
            n_occurrences = [sum((length(ref.occurrences[h.reference]) for h in r.hits); init = 0)
                             for r in results])
end

"""
    hit_table(ref, results)

One row per hit: `query`, `hit` (rank of the hit for its query),
`reference`, `identity`, `coverage`, `strand`, `query_start`, `query_end`,
`reference_start`, `reference_end`, `cigar`, `taxon` (deepest name of the
reference), one column per rank of `ref` and `n_occurrences`.
"""
function hit_table(ref::ReferenceDataset, results::AbstractVector{ASVResult})
    rows = [(r, j, h) for r in results for (j, h) in enumerate(r.hits)]
    fixed = (query = [r.id for (r, _, _) in rows],
             hit = [j for (_, j, _) in rows],
             reference = [ref.ids[h.reference] for (_, _, h) in rows],
             identity = [h.identity for (_, _, h) in rows],
             coverage = [h.coverage for (_, _, h) in rows],
             strand = [h.strand for (_, _, h) in rows],
             query_start = [h.query_start for (_, _, h) in rows],
             query_end = [h.query_end for (_, _, h) in rows],
             reference_start = [h.target_start for (_, _, h) in rows],
             reference_end = [h.target_end for (_, _, h) in rows],
             cigar = [h.cigar for (_, _, h) in rows],
             taxon = [taxon(ref, h.reference) for (_, _, h) in rows])
    ranks = NamedTuple{Tuple(Symbol.(ref.ranks))}(
        Tuple([ref.lineages[h.reference][k] for (_, _, h) in rows] for k in eachindex(ref.ranks)))
    occ = (n_occurrences = [length(ref.occurrences[h.reference]) for (_, _, h) in rows],)
    return merge(fixed, ranks, occ)
end

"""
    occurrence_table(ref, results)

One row per occurrence of each hit, for mapping software: `query`,
`reference`, `identity`, `taxon`, `latitude`, `longitude` and `fields` (the
other metadata of the occurrence, `name=value` pairs separated by `; `).
"""
function occurrence_table(ref::ReferenceDataset, results::AbstractVector{ASVResult})
    rows = [(r, h, o) for r in results for h in r.hits for o in ref.occurrences[h.reference]]
    return (query = [r.id for (r, _, _) in rows],
            reference = [ref.ids[h.reference] for (_, h, _) in rows],
            identity = [h.identity for (_, h, _) in rows],
            taxon = [taxon(ref, h.reference) for (_, h, _) in rows],
            latitude = [o.latitude for (_, _, o) in rows],
            longitude = [o.longitude for (_, _, o) in rows],
            fields = [_join_fields(o.fields) for (_, _, o) in rows])
end

"""
    write_table(path_or_io, table)

Write a column table (a `NamedTuple` of equally long vectors, such as
[`summary_table`](@ref)) as tab-separated text with a header. Tabs and line
breaks inside values become spaces, `missing` an empty field.
"""
function write_table(io::IO, table::NamedTuple)
    println(io, join(string.(keys(table)), '\t'))
    n = isempty(table) ? 0 : length(first(table))
    for i in 1:n
        println(io, join((_tsv_field(col[i]) for col in table), '\t'))
    end
    return nothing
end

write_table(path::AbstractString, table::NamedTuple) = open(io -> write_table(io, table), path, "w")
