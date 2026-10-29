# Summary of the usage log of ASVAtlas (`--usage-log`, ASVATLAS_USAGE_LOG):
# searches, query ASVs and networks by day or month, similarity levels,
# species and methods of the networks, search times.
#     julia bin/asvatlas-usage.jl usage.tsv
# In the container:
#     docker compose exec asvatlas julia /opt/ASVAtlas.jl/bin/asvatlas-usage.jl /opt/stats/usage.tsv

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using ASVAtlas

if length(ARGS) != 1
    println(stderr, "usage: julia bin/asvatlas-usage.jl USAGE_LOG")
    exit(1)
end
print(usage_report(ARGS[1]))
