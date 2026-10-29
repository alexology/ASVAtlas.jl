# Starts the ASVAtlas web application; see `julia bin/asvatlas.jl --help`.
#     julia -t auto bin/asvatlas.jl --fasta reference.fasta --metadata metadata.tsv

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using ASVAtlas

exit(ASVAtlas.main(ARGS))
