# Usage statistics of the web application: one line per search or network in
# a tab-separated file (`serve(...; usage_log)`), with counts only: no address
# of the user, no sequence, no coordinate. `usage_report` sums the file up.

const _USAGE_COLUMNS = ["time", "event", "queries", "found", "identity", "seconds", "species",
                        "method", "sequences"]
const _usage_path = Ref("")
const _usage_lock = ReentrantLock()
const _usage_warned = Ref(false)

_usage_text(x) = replace(string(x), r"[\t\r\n]" => " ")

# Appends one event (nothing without a usage log). A failure to write is
# reported once and never stops the request.
function _log_usage(event::AbstractString; queries = 0, found = 0, identity = "", seconds = 0.0,
                    species = "", method = "", sequences = 0)
    path = _usage_path[]
    isempty(path) && return nothing
    row = [Libc.strftime("%Y-%m-%dT%H:%M", time()), event, queries, found, identity,
           round(seconds; digits = 2), species, method, sequences]
    try
        Base.@lock _usage_lock begin
            header = !isfile(path) || filesize(path) == 0
            open(path, "a") do io
                header && println(io, join(_USAGE_COLUMNS, '\t'))
                println(io, join(_usage_text.(row), '\t'))
            end
        end
    catch err
        _usage_warned[] || @warn "Usage log not written: $(sprint(showerror, err))"
        _usage_warned[] = true
    end
    return nothing
end

"""
    usage_report(path) -> String

Summary of the usage log written by [`serve`](@ref) with `usage_log = path`:
period, searches and query ASVs, networks, by month (or by day within two
months), similarity levels chosen, species and methods of the networks, and
search times. The log holds counts only (no address, sequence or coordinate
of the users).
"""
function usage_report(path::AbstractString)
    lines = filter(!isempty, readlines(path))
    (isempty(lines) || split(lines[1], '\t') != _USAGE_COLUMNS) &&
        throw(ArgumentError("$path is not an ASVAtlas usage log"))
    rows = [Dict(zip(_USAGE_COLUMNS, split(l, '\t'))) for l in lines[2:end]]
    rows = filter(r -> length(r) == length(_USAGE_COLUMNS), rows)
    io = IOBuffer()
    if isempty(rows)
        println(io, "No search or network yet.")
        return String(take!(io))
    end
    int(x) = something(tryparse(Int, x), 0)
    searches = filter(r -> r["event"] == "search", rows)
    networks = filter(r -> r["event"] == "network", rows)
    first_day, last_day = first(rows)["time"][1:10], last(rows)["time"][1:10]
    println(io, "ASVAtlas usage, $first_day to $last_day")
    println(io)
    println(io, "Searches:           ", length(searches))
    println(io, "Query ASVs:         ", sum(r -> int(r["queries"]), searches; init = 0))
    println(io, "ASVs with matches:  ", sum(r -> int(r["found"]), searches; init = 0))
    println(io, "Networks:           ", length(networks))
    # By day within two months, else by month.
    days = sort(unique(r["time"][1:10] for r in rows))
    by_day = length(days) <= 62
    key(r) = by_day ? r["time"][1:10] : r["time"][1:7]
    println(io)
    println(io, rpad(by_day ? "Day" : "Month", 12), lpad("searches", 10), lpad("ASVs", 10), lpad("networks", 10))
    for k in sort(unique(key.(rows)))
        s = filter(r -> key(r) == k, searches)
        println(io, rpad(k, 12), lpad(length(s), 10), lpad(sum(r -> int(r["queries"]), s; init = 0), 10),
                lpad(count(r -> key(r) == k, networks), 10))
    end
    table(title, values; n = 10) = begin
        isempty(values) && return
        counts = Dict{String,Int}()
        for v in values
            counts[v] = get(counts, v, 0) + 1
        end
        println(io)
        println(io, title)
        for (v, c) in first(sort!(collect(counts); by = x -> (-x[2], x[1])), n)
            println(io, "  ", rpad(v, 40), lpad(c, 6))
        end
    end
    table("Similarity of the searches (%)", [r["identity"] for r in searches])
    table("Species of the networks", [r["species"] for r in networks])
    table("Network methods", [uppercase(r["method"]) for r in networks])
    times = sort([something(tryparse(Float64, r["seconds"]), 0.0) for r in searches])
    if !isempty(times)
        println(io)
        println(io, "Search time: median $(times[(length(times) + 1) ÷ 2]) s, longest $(last(times)) s")
    end
    return String(take!(io))
end
