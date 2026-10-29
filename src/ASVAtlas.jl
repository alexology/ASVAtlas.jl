"""
    ASVAtlas

Search amplicon sequence variants (ASVs) against a reference set of ASVs whose
taxonomy and sampling locations are already known, and show the hits on a map.

[`load_reference`](@ref) reads the reference sequences (FASTA) and a metadata
table with their taxonomy and coordinates; [`search_asvs`](@ref) compares the
query ASVs with the references at a chosen percent identity (100, 97, ...)
and reports, for each hit, the taxonomy and the occurrences of the reference.
[`taxon_network`](@ref) builds the haplotype network of the reference ASVs of
a taxon (ProfileAlignments.jl and HaplotypeNetworks.jl), with the places where
each haplotype was found. [`serve`](@ref) starts a web application (Genie.jl)
with the search and the networks, on a world map. The search itself is done
by NucleotideSearch.jl.
"""
module ASVAtlas

using FASTX
import Genie
import HaplotypeNetworks as HN
using HTTP
using JSON3
using Logging
using Logging: AbstractLogger, ConsoleLogger, NullLogger, global_logger, with_logger
using NucleotideSearch
import ProfileAlignments as PA
using Sockets

export ReferenceDataset, Occurrence, ReferenceHit, ASVResult, IDENTITY_LEVELS, KNOWN_RANKS,
       load_reference, parse_queries, read_queries, search_asvs, best_taxa, taxon,
       summary_table, hit_table, occurrence_table, write_table, serve, shutdown,
       TaxonNetwork, NETWORK_METHODS, taxon_network, taxon_references, taxon_list,
       network_sites, reference_sequence, usage_report

include("table.jl")
include("reference.jl")
include("search.jl")
include("output.jl")
include("networks.jl")
include("usage.jl")
include("server.jl")
include("cli.jl")

end
