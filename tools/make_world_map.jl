# Writes src/assets/world.geojson, the offline world map of the web
# application, from the Natural Earth 1:50m admin-0 countries (public domain):
#     https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson
# Only the country name is kept and coordinates are rounded to 3 decimals
# (about 100 m), well below the detail of the 1:50m scale.
# Usage: julia --project=. tools/make_world_map.jl <ne_50m_admin_0_countries.geojson>

using JSON3

length(ARGS) == 1 || error("usage: make_world_map.jl <ne_50m_admin_0_countries.geojson>")
src = JSON3.read(read(ARGS[1], String))

# Rounds the coordinates and drops consecutive repeated points of each ring.
round_ring(ring) = foldl(ring; init = Vector{Vector{Float64}}()) do acc, p
    q = [round(Float64(p[1]); digits = 3), round(Float64(p[2]); digits = 3)]
    (isempty(acc) || last(acc) != q) && push!(acc, q)
    acc
end

function geometry(g)
    if g.type == "Polygon"
        rings = [round_ring(r) for r in g.coordinates]
    elseif g.type == "MultiPolygon"
        rings = [[round_ring(r) for r in poly] for poly in g.coordinates]
    else
        error("unexpected geometry $(g.type)")
    end
    return Dict("type" => g.type, "coordinates" => rings)
end

features = [Dict("type" => "Feature",
                 "properties" => Dict("name" => String(f.properties.NAME)),
                 "geometry" => geometry(f.geometry))
            for f in src.features if f.geometry !== nothing]
out = joinpath(@__DIR__, "..", "src", "assets", "world.geojson")
open(out, "w") do io
    JSON3.write(io, Dict("type" => "FeatureCollection",
                         "attribution" => "Natural Earth 1:50m admin-0 countries (public domain)",
                         "features" => features))
end
println("Wrote ", length(features), " countries, ", filesize(out), " bytes to ", normpath(out))
