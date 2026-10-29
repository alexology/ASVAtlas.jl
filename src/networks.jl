# Haplotype networks of the reference ASVs of one taxon (HaplotypeNetworks.jl),
# after a multiple alignment (ProfileAlignments.jl), with the places where
# each haplotype was found.

"""
    NETWORK_METHODS

Network methods of [`taxon_network`](@ref): `:mjn` (median-joining),
`:tcs` (statistical parsimony), `:msn` (minimum spanning network) and `:mst`
(minimum spanning tree), from HaplotypeNetworks.jl.
"""
const NETWORK_METHODS = (:mjn, :tcs, :msn, :mst)

"""
    TaxonNetwork

Haplotype network of the reference ASVs of a taxon, built by
[`taxon_network`](@ref). Fields: `rank`, `taxon`, `method`, `ids` (input
sequences: the reference ASVs, then the query ASVs), `reference` (index in
the [`ReferenceDataset`](@ref) of each input sequence, 0 for a query),
`alignment` (aligned sequences), `network` (a HaplotypeNetworks
`HaplotypeNetwork`; nodes `1:nhaplotypes` are the observed haplotypes),
`layout` (node coordinates, `nnodes × 2`, in mutational steps) and
`reversed` (whether each sequence was reverse-complemented to match the
orientation of the others) and `excluded` (identifiers of the sequences left
out because they share no 8-mer with the others on either strand). `ids`,
`reference`, `alignment` and `reversed` are in the order given to
HaplotypeNetworks: by decreasing frequency of their haplotype.
"""
struct TaxonNetwork
    rank::String
    taxon::String
    method::Symbol
    ids::Vector{String}
    reference::Vector{Int}
    alignment::Vector{String}
    network::HN.HaplotypeNetwork
    layout::Matrix{Float64}
    reversed::Vector{Bool}
    excluded::Vector{String}
end

function Base.show(io::IO, tn::TaxonNetwork)
    print(io, "TaxonNetwork(", tn.taxon, ", ", uppercase(String(tn.method)), ", ",
          HN.nhaplotypes(tn.network.hapset), " haplotypes from ", length(tn.ids), " sequences, ",
          HN.nnodes(tn.network), " nodes)")
end

"""
    reference_sequence(ref, i) -> String

Sequence of reference ASV `i` as stored in the search index (bases other than
A, C, G and T are read back as `N`).
"""
function reference_sequence(ref::ReferenceDataset, i::Integer)
    db = ref.db
    codes = view(db.data, db.starts[i]:db.starts[i + 1] - 1)
    return String([c < 0x04 ? _BASES[c + 1] : 'N' for c in codes])
end

"""
    taxon_references(ref, rank, taxon) -> Vector{Int}

Indices of the reference ASVs whose lineage has `taxon` at `rank`.
"""
function taxon_references(ref::ReferenceDataset, rank::AbstractString, taxon::AbstractString)
    k = findfirst(==(lowercase(rank)), ref.ranks)
    k === nothing && throw(ArgumentError("the reference has no rank \"$rank\""))
    return [i for i in eachindex(ref.ids) if ref.lineages[i][k] == taxon]
end

"""
    taxon_list(ref; min_references = 1) -> Dict{String,Vector{Pair{String,Int}}}

Taxa of each rank of `ref` with their number of reference ASVs (at least
`min_references`), sorted by name.
"""
function taxon_list(ref::ReferenceDataset; min_references::Integer = 1)
    out = Dict{String,Vector{Pair{String,Int}}}()
    for (k, rank) in enumerate(ref.ranks)
        counts = Dict{String,Int}()
        for l in ref.lineages
            isempty(l[k]) || (counts[l[k]] = get(counts, l[k], 0) + 1)
        end
        out[rank] = sort!([n => c for (n, c) in counts if c >= min_references]; by = first)
    end
    return out
end

"""
    taxon_network(ref, rank, taxon; queries = Pair{String,String}[], method = :mjn,
                  epsilon = 0, max_sequences = 300, threads = Threads.nthreads())
        -> TaxonNetwork

Haplotype network of the reference ASVs with `taxon` at `rank` (e.g.
`"species", "Fictiferus alpinus"`), optionally with query ASVs (`id =>
sequence`) placed among them.

The sequences are first put in the same orientation: starting from the one
linked to most others, each sequence is kept or reverse-complemented
according to which strand shares more 8-mers with the sequences already
oriented (an ASV hit on the minus strand by [`search_asvs`](@ref) is thus
turned around). A sequence sharing no 8-mer with the others on either strand
cannot be oriented or compared: it is left out, with a warning, and listed in
`excluded`. The others are then aligned with
ProfileAlignments.jl (DNA scoring) and collapsed into haplotypes by
HaplotypeNetworks.jl. As in PopART, alignment columns with a gap or an
ambiguous base in any sequence are left out, so ASVs of different lengths
are compared over their common region and indels do not count as
differences. Haplotypes are numbered by decreasing total frequency (`H1` is
the most frequent; ties keep the order of their most frequent member); the
frequency of a reference ASV is its number of occurrences (at least 1), a
query adds 0, so a haplotype seen only in queries has frequency 0.

`method` is one of [`NETWORK_METHODS`](@ref); `epsilon` is the tolerance of
`:mjn` and `:msn`. At most `max_sequences` sequences are accepted (the
median-joining network grows quickly with the number of haplotypes).
"""
function taxon_network(ref::ReferenceDataset, rank::AbstractString, taxon::AbstractString;
                       queries::AbstractVector = Pair{String,String}[], method::Symbol = :mjn,
                       epsilon::Integer = 0, max_sequences::Integer = 300,
                       threads::Integer = Threads.nthreads())
    method in NETWORK_METHODS ||
        throw(ArgumentError("method must be one of $(join(NETWORK_METHODS, ", "))"))
    0 <= epsilon <= 100 || throw(ArgumentError("epsilon must be between 0 and 100"))
    refs = taxon_references(ref, rank, taxon)
    isempty(refs) && throw(ArgumentError("no reference ASV has $rank \"$taxon\""))
    qs = _queries(queries)
    n = length(refs) + length(qs)
    n <= max_sequences ||
        throw(ArgumentError("$n sequences for $taxon; the limit is $max_sequences"))
    weight(i) = max(1, length(ref.occurrences[i]))
    sort!(refs; by = i -> (-weight(i), ref.ids[i]))
    ids = [ref.ids[refs]; first.(qs)]
    seqs = [reference_sequence.(Ref(ref), refs); last.(qs)]
    reference = [refs; zeros(Int, length(qs))]
    weights = Float64[weight.(refs); zeros(length(qs))]
    seqs, reversed, out = _orient(seqs)
    excluded = ids[out]
    keep = .!out
    ids, seqs, reference, weights, reversed =
        ids[keep], seqs[keep], reference[keep], weights[keep], reversed[keep]
    isempty(excluded) ||
        @warn "Left out of the network of $taxon (no 8-mer shared with the others): " *
              join(excluded, ", ")
    aligned = _align(seqs, threads)
    # HaplotypeNetworks numbers haplotypes in order of first appearance: put
    # the members of the most frequent haplotypes first, then collapse again.
    first_pass = HN.haplotypes(ids, aligned; weights)
    freq = HN.haplotype_frequencies(first_pass)
    p = sortperm(eachindex(ids); by = j -> (-freq[first_pass.assignment[j]],
                                            first_pass.assignment[j], j))
    ids, reference, aligned, weights, reversed =
        ids[p], reference[p], aligned[p], weights[p], reversed[p]
    hs = HN.haplotypes(ids, aligned; weights)
    net = method === :mjn ? HN.mjn(hs; epsilon) :
          method === :msn ? HN.msn(hs; epsilon) :
          method === :tcs ? HN.tcs(hs) : HN.mst(hs)
    layout = HN.layout_network(net; max_edge_length = 12)
    return TaxonNetwork(String(lowercase(rank)), String(taxon), method, ids, reference, aligned,
                        net, layout, reversed, excluded)
end

const _COMPLEMENT = Dict(zip("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))

_revcomp(s::AbstractString) = String([_COMPLEMENT[c] for c in Iterators.reverse(s)])

function _kmers(s::AbstractString, k::Int = 8)
    out = Set{String}()
    for i in 1:length(s) - k + 1
        w = s[i:i + k - 1]
        all(in("ACGT"), w) && push!(out, w)
    end
    return out
end

# Orientation of the sequences: starting from the one that shares 8-mers with
# most of the others (on either strand), the sequence sharing most 8-mers with
# those already oriented is added next, on the strand that shares more.
# Returns the oriented sequences, whether each was reversed, and whether each
# was left out because it shares no 8-mer with the others on either strand.
function _orient(seqs::Vector{String})
    n = length(seqs)
    fwd = [_kmers(s) for s in seqs]
    rev = [_kmers(_revcomp(s)) for s in seqs]
    shared(a, b) = count(in(b), a)
    linked = [count(i -> i != j && (shared(fwd[j], fwd[i]) > 0 || shared(fwd[j], rev[i]) > 0), 1:n)
              for j in 1:n]
    anchor = argmax(linked)            # the first one on ties
    out = copy(seqs)
    reversed = falses(n)
    excluded = falses(n)
    done = falses(n)
    done[anchor] = true
    seen = copy(fwd[anchor])
    for _ in 2:n
        best, jbest, revbest = -1, 0, false
        for j in 1:n
            done[j] && continue
            nf, nr = shared(fwd[j], seen), shared(rev[j], seen)
            if max(nf, nr) > best
                best, jbest, revbest = max(nf, nr), j, nr > nf
            end
        end
        if best == 0
            excluded .= .!done
            break
        end
        done[jbest] = true
        if revbest
            out[jbest] = _revcomp(seqs[jbest])
            reversed[jbest] = true
        end
        union!(seen, revbest ? rev[jbest] : fwd[jbest])
    end
    return out, collect(reversed), collect(excluded)
end

function _align(seqs::Vector{String}, threads)
    length(seqs) == 1 && return seqs
    recs = [PA.Sequence("s$i", s) for (i, s) in enumerate(seqs)]
    aln = PA.align_sequences(recs; model = PA.ScoringModel(alphabet = :dna), strategy = :auto,
                             threads = max(1, threads))
    byid = Dict(r.id => r.sequence for r in aln.records)
    return [byid["s$i"] for i in eachindex(seqs)]
end

"""
    network_sites(ref, tn::TaxonNetwork)

Places where the haplotypes of `tn` were found: one named tuple per distinct
coordinate pair, with `latitude`, `longitude`, `counts` (haplotype index =>
number of occurrences there) and `occurrences` (`(reference, fields)` for
each occurrence).
"""
function network_sites(ref::ReferenceDataset, tn::TaxonNetwork)
    hs = tn.network.hapset
    sites = Dict{Tuple{Float64,Float64},Int}()
    out = NamedTuple{(:latitude, :longitude, :counts, :occurrences),
                     Tuple{Float64,Float64,Dict{Int,Int},Vector{Tuple{String,Vector{Pair{String,String}}}}}}[]
    for (j, i) in enumerate(tn.reference)
        i == 0 && continue
        h = hs.assignment[j]
        for o in ref.occurrences[i]
            key = (o.latitude, o.longitude)
            k = get!(sites, key) do
                push!(out, (latitude = o.latitude, longitude = o.longitude, counts = Dict{Int,Int}(),
                            occurrences = Tuple{String,Vector{Pair{String,String}}}[]))
                length(out)
            end
            out[k].counts[h] = get(out[k].counts, h, 0) + 1
            push!(out[k].occurrences, (ref.ids[i], o.fields))
        end
    end
    return out
end

# Deepest known name of reference `i` with its rank.
function _deepest_taxon(ref::ReferenceDataset, i::Integer)
    k = findlast(!isempty, ref.lineages[i])
    return (name = ref.lineages[i][k], rank = ref.ranks[k])
end

function _network_json(ref::ReferenceDataset, tn::TaxonNetwork)
    net = tn.network
    hs = net.hapset
    freq = HN.haplotype_frequencies(hs)
    radii = HN.node_radii(net)
    nh = HN.nhaplotypes(hs)
    haplotypes = [(name = hs.hap_names[h], frequency = freq[h],
                   references = [tn.ids[j] for j in hs.members[h] if tn.reference[j] != 0],
                   queries = [tn.ids[j] for j in hs.members[h] if tn.reference[j] == 0],
                   taxa = unique!([_deepest_taxon(ref, tn.reference[j])
                                   for j in hs.members[h]
                                   if tn.reference[j] != 0 && !isempty(taxon(ref, tn.reference[j]))]),
                   occurrences = sum((length(ref.occurrences[tn.reference[j]])
                                      for j in hs.members[h] if tn.reference[j] != 0); init = 0))
                  for h in 1:nh]
    nodes = [(name = net.names[v], haplotype = net.haplotype[v], x = tn.layout[v, 1],
              y = tn.layout[v, 2], r = radii[v]) for v in 1:HN.nnodes(net)]
    edges = [(source = net.src[e], target = net.dst[e], steps = net.steps[e])
             for e in 1:HN.nedges(net)]
    sites = [(lat = s.latitude, lon = s.longitude,
              counts = [(haplotype = h, n = c) for (h, c) in sort!(collect(s.counts))],
              occurrences = [(reference = r, fields = [(name = k, value = v) for (k, v) in f])
                             for (r, f) in s.occurrences])
             for s in network_sites(ref, tn)]
    return (rank = tn.rank, taxon = tn.taxon, method = String(tn.method),
            sequences = length(tn.ids), alignment_length = hs.alignment_length,
            sites_used = length(hs.sites), variable_sites = HN.nvariable(hs),
            components = maximum(HN.network_components(net); init = 0),
            reversed = tn.ids[tn.reversed], excluded = tn.excluded,
            # Number of differences between observed haplotypes (rows of the matrix).
            distances = (D = HN.distance_matrix(hs); [D[h, :] for h in 1:nh]),
            haplotypes = haplotypes, nodes = nodes, edges = edges, sites = sites)
end
