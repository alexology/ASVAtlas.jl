# Starts the web application on the synthetic example data set, from a
# terminal:
#     julia -t auto examples/run_app.jl [port]
# or from a REPL (e.g. VS Code): include("examples/run_app.jl").
# Then open http://127.0.0.1:8080 and paste examples/data/queries.fasta.
# Stop the server with Ctrl-C.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using ASVAtlas

data = joinpath(@__DIR__, "data")
ref = load_reference(joinpath(data, "reference.fasta");
                     metadata = joinpath(data, "reference_metadata.tsv"),
                     name = "Synthetic example")
port = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8080
# The button "Try an example" of the page fills in these queries.
serve(ref; port, example = read(joinpath(data, "queries.fasta"), String))
