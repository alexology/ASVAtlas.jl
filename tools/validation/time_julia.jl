# Validation, step 7 (docs/VALIDATION.md): time of the network algorithms
# alone (HaplotypeNetworks.jl), on the alignments written by export_julia.jl,
# after a first run for compilation (best of 3).
#   julia tools/validation/time_julia.jl WORKDIR
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io = devnull)
import HaplotypeNetworks as HN

dir = joinpath(ARGS[1], "networks")
function readfasta(path)
    ids, seqs = String[], String[]
    for l in eachline(path)
        startswith(l, '>') ? (push!(ids, l[2:end]); push!(seqs, "")) : (seqs[end] *= strip(l))
    end
    return ids, seqs
end
files = sort([f[1:4] for f in readdir(dir) if endswith(f, ".members.tsv")])
species = Dict(split(l, '\t')[2] => split(l, '\t')[1] for l in Iterators.drop(eachline(joinpath(dir, "summary_julia.tsv")), 1))
open(joinpath(dir, "times_julia.tsv"), "w") do io
    println(io, "species\thaplotypes\tmst_s\tmsn_s\tmjn_s\ttcs_s")
    for f in files
        ids, seqs = readfasta(joinpath(dir, "$f.fasta"))
        w = [parse(Float64, split(l, '\t')[2]) for l in Iterators.drop(eachline(joinpath(dir, "$f.members.tsv")), 1)]
        hs = HN.haplotypes(ids, seqs; weights = w)
        t = Float64[]
        for m in (HN.mst, HN.msn, HN.mjn, HN.tcs)
            m(hs)                                  # compilation
            push!(t, minimum(@elapsed(m(hs)) for _ in 1:3))
        end
        println(io, species[f], '\t', HN.nhaplotypes(hs), '\t', join(round.(t; digits = 4), '\t'))
    end
end
print(read(joinpath(dir, "times_julia.tsv"), String))
