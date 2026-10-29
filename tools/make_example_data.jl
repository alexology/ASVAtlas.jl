# Writes the synthetic example data set in examples/data: reference ASVs of
# invented taxa with invented sampling sites, and query ASVs at various
# distances from them. Run with `julia --project=. tools/make_example_data.jl`.

using Random

const OUT = joinpath(@__DIR__, "..", "examples", "data")
const BASES = ('A', 'C', 'G', 'T')

rng = Xoshiro(20261025)
random_seq(n) = String([rand(rng, BASES) for _ in 1:n])

function mutate(s, n)
    c = collect(s)
    for i in randperm(rng, length(c))[1:n]
        c[i] = rand(rng, filter(!=(c[i]), BASES))
    end
    return String(c)
end

revcomp(s) = String([Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')[b] for b in reverse(s)])

# (family, genus, species, number of reference ASVs, centre of its range)
const TAXA = [
    ("Fictiferidae", "Fictiferus", "Fictiferus alpinus", 3, (46.5, 10.5)),
    ("Fictiferidae", "Fictiferus", "Fictiferus lacustris", 3, (45.7, 8.6)),
    ("Fictiferidae", "Fictiferus", "Fictiferus sp.", 1, (44.4, 11.3)),
    ("Fictiferidae", "Pseudofictus", "Pseudofictus minor", 2, (42.0, 13.5)),
    ("Exemplidae", "Exemplus", "Exemplus primus", 3, (40.8, 16.0)),
    ("Exemplidae", "Exemplus", "Exemplus secundus", 2, (38.1, 14.0)),
]

len = 420
family_seq = Dict(f => random_seq(len) for f in unique(first.(TAXA)))
genus_seq = Dict(g => mutate(family_seq[f], 40) for (f, g) in unique((t[1], t[2]) for t in TAXA))

refs = Pair{String,String}[]
meta = Vector{Vector{String}}()
habitats = ["stream", "lake", "pond", "spring"]
for (f, g, sp, n, (lat, lon)) in TAXA
    base = mutate(genus_seq[g], 25)
    for v in 1:n
        id = "REF_" * lpad(length(refs) + 1, 3, '0')
        push!(refs, id => (v == 1 ? base : mutate(base, rand(rng, 1:5))))
        taxonomy = ["Animalia", "Demonstrata", "Syntheticae", "Exemplaria", f, g, sp]
        # The last variant of each species has no known location.
        if v == n && n > 1
            push!(meta, [id; taxonomy; ""; ""; ""; ""; ""])
            continue
        end
        for k in 1:rand(rng, 1:4)
            site = "Site " * string(Char('A' + length(meta) % 26)) * lpad(k, 2, '0')
            push!(meta, [id; taxonomy;
                         string(round(lat + randn(rng) * 0.6; digits = 5));
                         string(round(lon + randn(rng) * 0.6; digits = 5));
                         site; rand(rng, habitats);
                         string("2025-0", rand(rng, 4:9), "-", lpad(rand(rng, 1:28), 2, '0'))])
        end
    end
end

queries = [
    "ASV_identical" => refs[1].second,
    "ASV_99" => mutate(refs[4].second, 3),
    "ASV_97" => mutate(refs[11].second[11:410], 10),
    "ASV_reverse" => revcomp(refs[8].second),
    "ASV_genus" => mutate(genus_seq["Exemplus"], 15),
    "ASV_unknown" => random_seq(len),
]

mkpath(OUT)
open(joinpath(OUT, "reference.fasta"), "w") do io
    for (id, s) in refs
        println(io, ">", id, "\n", s)
    end
end
open(joinpath(OUT, "reference_metadata.tsv"), "w") do io
    println(io, join(["id", "kingdom", "phylum", "class", "order", "family", "genus", "species",
                      "latitude", "longitude", "site", "habitat", "date"], '\t'))
    for r in meta
        println(io, join(r, '\t'))
    end
end
open(joinpath(OUT, "queries.fasta"), "w") do io
    for (id, s) in queries
        println(io, ">", id, "\n", s)
    end
end
println("Wrote ", length(refs), " references, ", length(meta), " metadata rows, ",
        length(queries), " queries to ", normpath(OUT))
