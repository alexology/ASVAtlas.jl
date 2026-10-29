# Container image of ASVAtlas. Build from the folder that contains
# ASVAtlas.jl and its unregistered dependencies NucleotideSearch.jl,
# HaplotypeNetworks.jl and ProfileAlignments.jl (see `[sources]` in
# ASVAtlas.jl/Project.toml):
#   docker build -f ASVAtlas.jl/deploy/Dockerfile -t asvatlas .
# Run with the reference data in /path/to/data:
#   docker run -d --name asvatlas -p 127.0.0.1:8080:8080 -v /path/to/data:/data:ro \
#     -e ASVATLAS_FASTA=/data/freshwater_invertebrates.fasta \
#     -e ASVATLAS_METADATA=/data/freshwater_invertebrates_metadata.tsv asvatlas
# See ASVAtlas.jl/deploy/docker-compose.yml for the same with Docker Compose.

FROM julia:1.12.2

RUN useradd --system --create-home asvatlas && mkdir /opt/depot && chown asvatlas /opt/depot
ENV JULIA_DEPOT_PATH=/opt/depot \
    JULIA_NUM_THREADS=auto \
    ASVATLAS_HOST=0.0.0.0 \
    ASVATLAS_PORT=8080

COPY --chown=asvatlas NucleotideSearch.jl /opt/NucleotideSearch.jl
COPY --chown=asvatlas HaplotypeNetworks.jl /opt/HaplotypeNetworks.jl
COPY --chown=asvatlas ProfileAlignments.jl /opt/ProfileAlignments.jl
COPY --chown=asvatlas ASVAtlas.jl /opt/ASVAtlas.jl
USER asvatlas
# The Manifest is resolved here (local Manifests are excluded from the build).
RUN julia --project=/opt/ASVAtlas.jl -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

EXPOSE 8080
CMD ["julia", "/opt/ASVAtlas.jl/bin/asvatlas.jl"]
