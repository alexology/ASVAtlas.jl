# Validation, step 1: queries, ASVAtlas searches and networks for the R
# scripts of this folder (see docs/VALIDATION.md).
#
#   julia -t 4 tools/validation/export_julia.jl REFERENCE.fasta METADATA.tsv WORKDIR \
#         [TEST_QUERIES.fasta] [SPECIES ...]
#
# Writes to WORKDIR: queries.fasta, julia_hits_app.tsv (search as in the page,
# maxrejects 64), julia_hits_exhaustive.tsv (maxrejects 0) and, in
# WORKDIR/networks, for each species the alignment used by ASVAtlas, the
# haplotype of each sequence, the distances between haplotypes and the edges
# of the MST, MSN, MJN and TCS networks, with summary_julia.tsv. Without
# SPECIES, the 3 largest species with at most 300 reference ASVs and the 4
# with the number of ASVs nearest to 15 are used.
#
# The files written contain reference sequences: keep WORKDIR private.
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io = devnull)
using ASVAtlas, Random
import HaplotypeNetworks as HN

length(ARGS) >= 3 || error("usage: export_julia.jl REFERENCE.fasta METADATA.tsv WORKDIR [TEST_QUERIES.fasta] [SPECIES ...]")
fasta, metadata, out = ARGS[1], ARGS[2], ARGS[3]
tests = length(ARGS) >= 4 && endswith(lowercase(ARGS[4]), r"\.(fa|fasta|fas)") ? ARGS[4] : nothing
named = ARGS[(tests === nothing ? 4 : 5):end]
mkpath(joinpath(out, "networks"))

ref = load_reference(fasta; metadata, name = splitext(basename(fasta))[1])

# ---- Queries -------------------------------------------------------------------
# The test queries, then 200 reference ASVs with 0-40 substitutions (100 to
# about 90% identity), every 7th with 1-3 indels and every 6th
# reverse-complemented, and 10 shuffled sequences that should match nothing.
rng = MersenneTwister(2026)
comp = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')
revcomp(s) = String([comp[c] for c in reverse(s)])
function mutate(rng, s, nsub, nindel)
    v = collect(s)
    for p in randperm(rng, length(v))[1:nsub]
        v[p] = rand(rng, filter(!=(v[p]), ['A', 'C', 'G', 'T']))
    end
    for _ in 1:nindel
        p = rand(rng, 20:length(v) - 20)
        rand(rng, Bool) ? deleteat!(v, p) : insert!(v, p, rand(rng, ['A', 'C', 'G', 'T']))
    end
    return String(v)
end
queries = tests === nothing ? Pair{String,String}[] : read_queries(tests)
for k in 1:200
    i = rand(rng, 1:length(ref))
    nsub = rand(rng, 0:40)
    nindel = k % 7 == 0 ? rand(rng, 1:3) : 0
    s = mutate(rng, reference_sequence(ref, i), nsub, nindel)
    rc = k % 6 == 0
    push!(queries, "q$(lpad(k, 3, '0'))_$(ref.ids[i])_s$(nsub)_i$(nindel)$(rc ? "_rc" : "")" =>
                   rc ? revcomp(s) : s)
end
for k in 1:10
    push!(queries, "shuffled_$k" => String(shuffle(rng, collect(reference_sequence(ref, rand(rng, 1:length(ref)))))))
end
open(joinpath(out, "queries.fasta"), "w") do io
    for (id, s) in queries
        println(io, ">", id, "\n", s)
    end
end

# ---- Searches: as in the page (maxrejects 64) and exhaustive ------------------
for (label, maxrejects) in (("app", 64), ("exhaustive", 0))
    t = @elapsed res = search_asvs(ref, queries; identity = 90, query_cov = 0.9, maxhits = 2000,
                                   maxrejects, strand = :both)
    open(joinpath(out, "julia_hits_$label.tsv"), "w") do io
        println(io, "query\treference\tidentity\tcoverage\tstrand\tref_start\tref_end")
        for r in res, h in r.hits
            println(io, join((r.id, ref.ids[h.reference], h.identity, h.coverage, h.strand,
                              h.target_start, h.target_end), '\t'))
        end
    end
    println("search $label: $(round(t; digits = 1)) s, ",
            sum(length(r.hits) for r in res), " hits, max per query ",
            maximum(length(r.hits) for r in res))
end

# ---- Networks ----------------------------------------------------------------
species = taxon_list(ref)["species"]
counts = Dict(species)
chosen = if isempty(named)
    big = first([n for (n, c) in sort(species; by = x -> -x[2]) if c <= 300], 3)
    small = first([n for (n, c) in sort(species; by = x -> (abs(x[2] - 15), x[1])) if !(n in big)], 4)
    [big; small]
else
    collect(named)
end
open(joinpath(out, "networks", "summary_julia.tsv"), "w") do io
    println(io, "species\tfile\treferences\tsequences\thaplotypes\tsites\tvariable\tmethod\tnodes\tinferred\tedges\ttotal_steps\tcomponents\tseconds")
    for (k, sp) in enumerate(chosen)
        file = "sp$(lpad(k, 2, '0'))"
        for method in (:mst, :msn, :mjn, :tcs)
            # Whole taxon_network: orientation, alignment, haplotypes, network.
            t = @elapsed tn = taxon_network(ref, "species", sp; method, threads = 4)
            net, hs = tn.network, tn.network.hapset
            if method === :mst
                open(joinpath(out, "networks", "$file.fasta"), "w") do f
                    for (id, s) in zip(tn.ids, tn.alignment)
                        println(f, ">", id, "\n", s)
                    end
                end
                open(joinpath(out, "networks", "$file.members.tsv"), "w") do f
                    println(f, "id\tweight\thaplotype\treversed")
                    for (j, id) in enumerate(tn.ids)
                        println(f, id, '\t', hs.weights[j], '\t', hs.assignment[j], '\t', tn.reversed[j])
                    end
                end
                D = HN.distance_matrix(hs)
                open(joinpath(out, "networks", "$file.dist.tsv"), "w") do f
                    for i in axes(D, 1)
                        println(f, join(D[i, :], '\t'))
                    end
                end
            end
            et = HN.edges_table(net)
            open(joinpath(out, "networks", "$file.$method.edges.tsv"), "w") do f
                println(f, "from\tto\tsteps")
                for e in eachindex(et.steps)
                    println(f, et.from[e], '\t', et.to[e], '\t', et.steps[e])
                end
            end
            nh = HN.nhaplotypes(hs)
            println(io, join((sp, file, counts[sp], length(tn.ids), nh, HN.nsites(hs), HN.nvariable(hs),
                              method, HN.nnodes(net), HN.nnodes(net) - nh, HN.nedges(net), sum(net.steps),
                              maximum(HN.network_components(net)), round(t; digits = 2)), '\t'))
            flush(io)
        end
        println("network $sp done")
    end
end
