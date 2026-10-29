# Web application built on Genie.jl: the page and its static files (Leaflet
# and the world map are bundled, so no internet access is needed), the dataset
# description and the search API.
#
# Genie routes the requests; the connection layer is kept here so that it can
# refuse oversized bodies before reading them (Genie reads whole bodies into
# memory), check the Host header, bound the number of connections and give
# every response, Genie's own error pages included, the same safe headers.

const _ASSETS = joinpath(@__DIR__, "assets")
const _ASSET_TYPES = Dict(".html" => "text/html; charset=utf-8",
                          ".js" => "text/javascript; charset=utf-8",
                          ".css" => "text/css; charset=utf-8",
                          ".geojson" => "application/geo+json",
                          ".png" => "image/png")
# The only files served, read when the package is compiled.
const _ASSET_FILES = ["index.html", "app.js", "network.js", "about.js", "style.css", "world.geojson",
                      "leaflet/leaflet.js", "leaflet/leaflet.css",
                      "leaflet/images/layers.png", "leaflet/images/layers-2x.png",
                      "dnambientale.png"]
for f in _ASSET_FILES
    include_dependency(joinpath(_ASSETS, f))
end
const _STATIC = Dict((f == "index.html" ? "/" : "/" * f) =>
                     (read(joinpath(_ASSETS, f)), _ASSET_TYPES[splitext(f)[2]])
                     for f in _ASSET_FILES)

# Everything comes from this server except the optional OpenStreetMap tiles.
# The tile servers refuse requests without a Referer, so the origin (never the
# path) is sent: "no-referrer" makes every tile "Access blocked".
const _CSP = "default-src 'self'; img-src 'self' data: https://tile.openstreetmap.org; " *
             "connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'"
const _SECURITY_HEADERS = ["X-Content-Type-Options" => "nosniff",
                           "Content-Security-Policy" => _CSP,
                           "Referrer-Policy" => "strict-origin",
                           "X-Frame-Options" => "DENY",
                           "Cross-Origin-Opener-Policy" => "same-origin",
                           "Cross-Origin-Resource-Policy" => "same-origin"]
const _STATUS_TEXT = Dict(400 => "bad request", 401 => "unauthorized", 403 => "forbidden",
                          404 => "not found", 405 => "method not allowed",
                          413 => "request too large", 415 => "unsupported media type",
                          421 => "unknown host", 500 => "internal error")

# Default lowest identity of the web API (`serve(...; min_identity)`): low
# values would let anyone list the reference by searching unrelated
# sequences, and make every reference a candidate with very wide bands.
const _MIN_WEB_IDENTITY = 90
# At most this many hit rows (queries × hits per query) in one response.
const _MAX_HIT_ROWS = 1_000_000

_response(status, body, type; cache = "no-store") =
    HTTP.Response(status, ["Content-Type" => type, "Cache-Control" => cache], body)

_json(status, x) = _response(status, JSON3.write(x), "application/json; charset=utf-8")
_error(status, msg) = _json(status, (error = msg,))

# Same security headers on every response. Error pages not produced by
# ASVAtlas (Genie's own 404, 401 or 500 echo the requested path and take the
# content type of the request) are replaced by a fixed JSON message.
function _harden!(res::HTTP.Response)
    if res.status >= 400 && !startswith(HTTP.header(res, "Content-Type"), "application/json")
        res.body = Vector{UInt8}(JSON3.write((error = get(_STATUS_TEXT, res.status, "error"),)))
        HTTP.setheader(res, "Content-Type" => "application/json; charset=utf-8")
        HTTP.setheader(res, "Cache-Control" => "no-store")
    end
    for h in _SECURITY_HEADERS
        HTTP.setheader(res, h)
    end
    HTTP.setheader(res, "Content-Length" => string(length(res.body)))
    return res
end

_occurrence_json(o::Occurrence) =
    (lat = o.latitude, lon = o.longitude, fields = [(name = k, value = v) for (k, v) in o.fields])

# Only summary figures: the page never receives the list of taxa or any
# reference sequence.
function _dataset_json(ref::ReferenceDataset; max_queries, max_hits, max_length,
                       max_network_sequences = 300, min_identity = _MIN_WEB_IDENTITY,
                       example = "")
    (name = ref.name, marker = ref.marker, primers = ref.primers,
     about = isempty(ref.about) ? nothing : JSON3.read(ref.about), sources = _source_counts(ref),
     references = length(ref), ranks = ref.ranks,
     georeferenced = count(!isempty, ref.occurrences),
     occurrences = sum(length, ref.occurrences; init = 0),
     identity_levels = [l for l in IDENTITY_LEVELS if l >= min_identity], min_identity = min_identity,
     max_queries = max_queries, max_hits = max_hits, max_length = max_length,
     networks = "species" in ref.ranks,
     network_methods = String.(collect(NETWORK_METHODS)),
     max_network_sequences = max_network_sequences, example = example)
end

# Records, reference ASVs and places of each value of the occurrence column
# `dataset` (the studies of the reference), for the tab "Reference data".
function _source_counts(ref::ReferenceDataset)
    records = Dict{String,Int}()
    asvs = Dict{String,Set{Int}}()
    places = Dict{String,Set{Tuple{Float64,Float64}}}()
    for (i, occs) in enumerate(ref.occurrences), o in occs
        k = findfirst(p -> lowercase(first(p)) == "dataset", o.fields)
        k === nothing && continue
        d = last(o.fields[k])
        records[d] = get(records, d, 0) + 1
        push!(get!(asvs, d, Set{Int}()), i)
        push!(get!(places, d, Set{Tuple{Float64,Float64}}()), (o.latitude, o.longitude))
    end
    return Dict(d => (records = n, asvs = length(asvs[d]), places = length(places[d])) for (d, n) in records)
end

_parse_json(payload) = try
    JSON3.read(payload)
catch
    throw(ArgumentError("the request body is not valid JSON"))
end

function _string(body, key, default = nothing)
    v = get(body, key, default)
    v isa AbstractString || throw(ArgumentError("$key must be a string"))
    return String(v)
end

# Networks are built only for a species found by the user's ASVs: the query
# ASVs listed in `include` (identifiers of the FASTA text `sequences`) are
# searched again at `identity`, and those that hit `taxon` are placed in its
# network. Without a hit in that species the request is refused, so the page
# cannot be used to browse the reference.
function _network_request(ref, payload::AbstractString, lock; max_network_sequences, threads,
                          max_length = typemax(Int), min_identity = _MIN_WEB_IDENTITY)
    body = _parse_json(payload)
    body isa JSON3.Object || throw(ArgumentError("the request body must be a JSON object"))
    "species" in ref.ranks ||
        throw(ArgumentError("haplotype networks need species names in the reference"))
    taxon = _string(body, :taxon, "")
    method = _string(body, :method, "mjn")
    method in String.(NETWORK_METHODS) ||
        throw(ArgumentError("method must be one of $(join(NETWORK_METHODS, ", "))"))
    epsilon = Int(_number(body, :epsilon, 0, 0, 10))
    identity = _number(body, :identity, min_identity, min_identity, 100)
    include = get(body, :include, String[])
    (include isa AbstractVector && all(x -> x isa AbstractString, include)) ||
        throw(ArgumentError("include must be a list of query identifiers"))
    isempty(include) &&
        throw(ArgumentError("a network needs at least one of your ASVs (include)"))
    wanted = Set(String.(include))
    queries = filter(q -> first(q) in wanted, parse_queries(_string(body, :sequences, "")))
    length(queries) == length(wanted) ||
        throw(ArgumentError("some identifiers in include are not in sequences"))
    length(queries) <= max_network_sequences ||
        throw(ArgumentError("too many ASVs for a network; the limit is $max_network_sequences"))
    for (id, s) in queries
        length(s) <= max_length ||
            throw(ArgumentError("query $id is $(length(s)) bases long; the limit is $max_length"))
    end
    k = findfirst(==("species"), ref.ranks)
    t = @elapsed tn = Base.@lock lock begin
        hits = search_asvs(ref, queries; identity, threads)
        matching = [q for (q, r) in zip(queries, hits)
                    if any(h -> ref.lineages[h.reference][k] == taxon, r.hits)]
        isempty(matching) &&
            throw(ArgumentError("none of your ASVs matches $taxon at ≥ $identity% identity"))
        taxon_network(ref, "species", taxon; queries = matching, method = Symbol(method), epsilon,
                      max_sequences = max_network_sequences, threads)
    end
    @info "Network: $(uppercase(method)) of $taxon, $(length(tn.ids)) sequences, " *
          "$(HN.nhaplotypes(tn.network.hapset)) haplotypes, $(round(t; digits = 2)) s"
    return _network_json(ref, tn)
end

# Search results with the references they hit, each reference once.
function _results_json(ref::ReferenceDataset, results, params)
    hit_refs = sort!(unique(h.reference for r in results for h in r.hits))
    references = Dict(ref.ids[i] => (length = ref.lengths[i], taxon = taxon(ref, i),
                                      lineage = ref.lineages[i],
                                      occurrences = _occurrence_json.(ref.occurrences[i]))
                      for i in hit_refs)
    queries = [(id = r.id, length = length(r.sequence), exhaustive = r.exhaustive,
                best_taxa = best_taxa(ref, r),
                hits = [(reference = ref.ids[h.reference], identity = h.identity,
                         coverage = h.coverage, strand = string(h.strand),
                         query_start = h.query_start, query_end = h.query_end,
                         reference_start = h.target_start, reference_end = h.target_end,
                         cigar = h.cigar) for h in r.hits])
               for r in results]
    return (params = params, ranks = ref.ranks, queries = queries, references = references)
end

function _number(body, key, default, lo, hi)
    v = get(body, key, default)
    v isa Real && !(v isa Bool) && isfinite(v) && lo <= v <= hi ||
        throw(ArgumentError("$key must be a number between $lo and $hi"))
    return v
end

function _search_request(ref, payload::AbstractString, lock; max_queries, max_hits, threads,
                         max_length = typemax(Int), min_identity = _MIN_WEB_IDENTITY)
    body = _parse_json(payload)
    body isa JSON3.Object || throw(ArgumentError("the request body must be a JSON object"))
    text = get(body, :sequences, "")
    text isa AbstractString || throw(ArgumentError("sequences must be a string"))
    queries = parse_queries(text)
    length(queries) <= max_queries ||
        throw(ArgumentError("too many queries ($(length(queries))); the limit is $max_queries"))
    for (id, s) in queries
        length(s) <= max_length ||
            throw(ArgumentError("query $id is $(length(s)) bases long; the limit is $max_length"))
    end
    identity = _number(body, :identity, 100, min_identity, 100)
    query_cov = _number(body, :query_cov, 0.9, 0.01, 1)
    maxhits = Int(round(_number(body, :max_hits, min(50, max_hits), 1, max_hits)))
    length(queries) * maxhits <= _MAX_HIT_ROWS ||
        throw(ArgumentError("too many possible hits ($(length(queries)) queries × $maxhits); " *
                            "lower the hits per ASV or split the queries"))
    strand = get(body, :strand, "both")
    strand in ("both", "plus") || throw(ArgumentError("strand must be \"both\" or \"plus\""))
    t = @elapsed results = Base.@lock lock search_asvs(ref, queries; identity, query_cov,
                                                       maxhits, strand = Symbol(strand), threads)
    @info "Search: $(length(queries)) quer$(length(queries) == 1 ? "y" : "ies") at ≥ $identity% " *
          "identity, $(count(r -> !isempty(r.hits), results)) with hits, $(round(t; digits = 2)) s"
    params = (identity = identity, query_cov = query_cov, max_hits = maxhits, strand = strand)
    return _results_json(ref, results, params)
end

const _BASES = ('A', 'C', 'G', 'T')

# Runs one search and its JSON conversion before serving, so that the first
# request does not wait for compilation.
function _warmup(ref::ReferenceDataset, threads)
    db = ref.db
    codes = db.data[db.starts[1]:min(db.starts[2] - 1, db.starts[1] + 299)]
    seq = String([c < 0x04 ? _BASES[c + 1] : 'N' for c in codes])
    t = @elapsed begin
        payload = JSON3.write((sequences = ">warmup\n" * seq, identity = 100))
        with_logger(NullLogger()) do
            JSON3.write(_search_request(ref, payload, ReentrantLock(); max_queries = 1,
                                        max_hits = 10, threads))
        end
        _harden!(_error(404, "warm-up"))
        _warmup_network(ref, threads)
    end
    @info "Warm-up search done in $(round(t; digits = 1)) s"
    return nothing
end

# Compiles the network code on the smallest species with at least two
# reference ASVs, with its first reference as the user's ASV.
function _warmup_network(ref::ReferenceDataset, threads)
    "species" in ref.ranks || return nothing
    taxa = filter(t -> 2 <= last(t) <= 50, taxon_list(ref)["species"])
    isempty(taxa) && return nothing
    name = first(argmin(last, taxa))
    i = first(taxon_references(ref, "species", name))
    with_logger(NullLogger()) do
        payload = JSON3.write((taxon = name, include = ["warmup"],
                               sequences = ">warmup\n" * reference_sequence(ref, i)))
        JSON3.write(_network_request(ref, payload, ReentrantLock(); max_network_sequences = 51,
                                     threads))
    end
    return nothing
end

# Bad input gives 400 with its message; other errors are logged and hidden.
function _guarded(f)
    try
        return f()
    catch err
        err isa ArgumentError && return _error(400, err.msg)
        @error "Request failed" exception = (err, catch_backtrace())
        return _error(500, "internal error, see the server log")
    end
end

_is_json(req) = startswith(lowercase(strip(HTTP.header(req, "Content-Type"))), "application/json")

function _register_routes(ref::ReferenceDataset; max_queries, max_hits, max_length, threads,
                          max_network_sequences = 300, min_identity = _MIN_WEB_IDENTITY,
                          example = "")
    lock = ReentrantLock()
    dataset = _dataset_json(ref; max_queries, max_hits, max_length, max_network_sequences,
                            min_identity, example)
    not_allowed() = _error(405, "method not allowed")
    for (path, (body, type)) in _STATIC
        Genie.Router.route(path; named = Symbol(:asvatlas_static, path)) do
            _response(200, body, type; cache = "no-cache")
        end
        Genie.Router.route(not_allowed, path; method = Genie.Router.POST,
                           named = Symbol(:asvatlas_static_post, path))
    end
    Genie.Router.route("/api/dataset"; named = :asvatlas_dataset) do
        _json(200, dataset)
    end
    Genie.Router.route("/health"; named = :asvatlas_health) do
        _json(200, (status = "ok", references = length(ref)))
    end
    Genie.Router.route("/api/search"; method = Genie.Router.POST, named = :asvatlas_search) do
        _guarded() do
            # JSON only: a cross-site page cannot send it without a CORS
            # preflight, which this server never allows.
            _is_json(Genie.Requests.request()) ||
                return _error(415, "the request body must be JSON (Content-Type: application/json)")
            payload = Genie.Requests.rawpayload()
            _json(200, _search_request(ref, payload, lock; max_queries, max_hits, threads,
                                       max_length, min_identity))
        end
    end
    Genie.Router.route("/api/network"; method = Genie.Router.POST, named = :asvatlas_network) do
        _guarded() do
            _is_json(Genie.Requests.request()) ||
                return _error(415, "the request body must be JSON (Content-Type: application/json)")
            payload = Genie.Requests.rawpayload()
            _json(200, _network_request(ref, payload, lock; max_network_sequences, threads,
                                        max_length, min_identity))
        end
    end
    Genie.Router.route(not_allowed, "/api/network"; named = :asvatlas_network_get)
    Genie.Router.route(not_allowed, "/api/search"; named = :asvatlas_search_get)
    Genie.Router.route(not_allowed, "/api/dataset"; method = Genie.Router.POST,
                       named = :asvatlas_dataset_post)
    return nothing
end

struct _RequestError <: Exception
    status::Int
    msg::String
end

# Host name of a Host header value, without port ("[::1]:8080" -> "::1").
function _host_name(h::AbstractString)
    h = lowercase(strip(h))
    startswith(h, '[') && return String(first(split(h[2:end], ']')))
    return String(first(split(h, ':')))
end

function _check_request(req::HTTP.Request, allowed_hosts, max_body)
    if !isempty(allowed_hosts)
        h = HTTP.header(req, "Host")
        (lowercase(strip(h)) in allowed_hosts || _host_name(h) in allowed_hosts) ||
            throw(_RequestError(421, "unknown host"))
    end
    cl = HTTP.header(req, "Content-Length")
    if !isempty(cl)
        n = tryparse(Int, strip(cl))
        (n === nothing || n < 0) && throw(_RequestError(400, "invalid Content-Length"))
        n > max_body && throw(_RequestError(413, "request larger than $max_body bytes"))
    end
    return nothing
end

# Reads the body, stopping as soon as it exceeds `max_body` (chunked bodies
# have no Content-Length to check beforehand).
function _read_body(stream::HTTP.Stream, max_body::Int)
    buf = IOBuffer()
    while !eof(stream)
        write(buf, readavailable(stream))
        position(buf) > max_body && throw(_RequestError(413, "request larger than $max_body bytes"))
    end
    return take!(buf)
end

function _stream_handler(allowed_hosts::Set{String}, max_body::Int)
    return function (stream::HTTP.Stream)
        req = stream.message
        close_after = false
        res = try
            _check_request(req, allowed_hosts, max_body)
            req.body = _read_body(stream, max_body)
            HTTP.closeread(stream)
            Genie.Server.setup_http_listener(req, HTTP.Response(); stream)
        catch err
            if err isa _RequestError
                close_after = err.status in (400, 413)   # unread body left on the connection
                _error(err.status, err.msg)
            else
                @error "Request failed" exception = (err, catch_backtrace())
                _error(500, "internal error, see the server log")
            end
        end
        _harden!(res)
        close_after && HTTP.setheader(res, "Connection" => "close")
        req.response = res
        res.request = req
        HTTP.startwrite(stream)
        write(stream, res.body)
        # A refused body is not read: close the connection rather than let
        # HTTP.jl read (and discard) the rest of it, and mark the request as
        # read so that HTTP.jl does not log the unread body as an error.
        if close_after
            close(stream.stream)
            stream.ntoread = 0
        end
        return
    end
end

"""
    serve(ref::ReferenceDataset; host = "127.0.0.1", port = 8080, async = false,
          max_queries = 10_000, max_hits = 500, max_length = 20_000,
          max_body = 16 * 2^20, max_connections = 128, allowed_hosts = String[],
          max_network_sequences = 300, min_identity = 90, threads = Threads.nthreads(),
          warmup = true)

Start the web application (a Genie.jl application) on `http://host:port`.
The page answers the question "where in the world is my ASV?": it accepts
query ASVs (pasted FASTA or a file) and a similarity (100% by default, i.e.
identical records; 99, 98, 97, ... down to `min_identity`), and maps the
occurrences of the matching reference ASVs, with their taxonomy. A second tab
builds the haplotype network ([`taxon_network`](@ref), at most
`max_network_sequences` sequences) of a species matched by the user's ASVs,
with those ASVs placed in it, and maps the haplotypes found at each place.

The reference itself is not published: no reference sequence and no list of
taxa is ever sent; matches need a similarity of at least `min_identity`
percent (90 by default), so the reference cannot be browsed with unrelated
sequences; networks are built only for species matched by the ASVs sent.

With `async = false` the call blocks until interrupted (Ctrl-C); with
`async = true` it returns the server, to be stopped with
[`shutdown`](@ref). Genie has one set of routes per Julia session, so one
reference dataset is served at a time. Searches run one at a time on
`threads` tasks. With `warmup = true` a search is run before serving, so
that the first user does not wait for compilation. Each search is logged
with its size and duration. `GET /health` answers `{"status": "ok", ...}`.

Limits against misuse: `max_queries` sequences per search, each of at most
`max_length` bases, `max_hits` hits per query, request bodies of at most
`max_body` bytes (refused before being read) and `max_connections` open
connections. `allowed_hosts` (e.g. `["asvatlas.example.org"]`) rejects
requests whose `Host` header names another host, which protects against DNS
rebinding; empty accepts any. The page and its files are served from a fixed
list, with a restrictive Content Security Policy.

`example` is FASTA text offered by the button "Try an example" of the page
(no button when empty). It is sent to every visitor: use sequences that may
be public, not reference sequences that should stay private.

The default `host` accepts only local connections; use `"0.0.0.0"` to serve
the network. There is no authentication: on a public server put it behind a
reverse proxy with HTTPS (see the README and `deploy/`). The page uses
relative URLs, so the proxy can publish it under a sub-path.

The world map (Natural Earth 1:50m countries, public domain) and Leaflet are
bundled, so the page works offline; OpenStreetMap tiles can be chosen in the
map as a more detailed, online alternative.

The API is `GET /api/dataset`, `POST /api/search` with `Content-Type:
application/json` and a JSON object `{"sequences": "<FASTA>", "identity": 100,
"query_cov": 0.9, "max_hits": 50, "strand": "both"}`, and `POST
/api/network` with `{"taxon": "<species>", "sequences": "<FASTA>", "include":
["ASV1"], "identity": 97, "method": "mjn", "epsilon": 0}` (the ASVs in
`include` that match the species at `identity` are placed in its network).
"""
function serve(ref::ReferenceDataset; host::AbstractString = "127.0.0.1", port::Integer = 8080,
               async::Bool = false, max_queries::Integer = 10_000, max_hits::Integer = 500,
               max_length::Integer = 20_000, max_body::Integer = 16 * 2^20,
               max_connections::Integer = 128, allowed_hosts = String[],
               max_network_sequences::Integer = 300, min_identity::Real = _MIN_WEB_IDENTITY,
               threads::Integer = Threads.nthreads(), warmup::Bool = true,
               example::AbstractString = "")
    0 < min_identity <= 100 || throw(ArgumentError("min_identity must be in (0, 100]"))
    sizeof(example) <= 2^16 || throw(ArgumentError("the example must be at most 64 KiB"))
    ip = string(Sockets.getaddrinfo(host))
    _check_port(ip, Int(port))
    warmup && _warmup(ref, threads)
    Genie.config.app_env = "prod"     # no Revise call on every request
    # Genie's static file handler is not used: the files are routes, from a
    # fixed list (see `_STATIC`), and the document root is an empty folder.
    Genie.config.server_document_root = mktempdir()
    Genie.config.server_handle_static_files = false
    Genie.config.server_signature = ""
    Genie.config.log_requests = false
    Genie.config.websockets_server = false
    _register_routes(ref; max_queries, max_hits, max_length, threads, max_network_sequences,
                     min_identity, example = String(example))
    hosts = Set{String}(lowercase(strip(String(h))) for h in allowed_hosts)
    listener = HTTP.listen!(_stream_handler(hosts, Int(max_body)), parse(Sockets.IPAddr, ip),
                            Int(port); max_connections = Int(max_connections))
    server = Genie.Server.ServersCollection(webserver = listener)
    @info "ASVAtlas: $(ref) at http://$(host == "0.0.0.0" ? "localhost" : host):$port"
    # Everything in the occurrences is public: say what that is.
    columns = unique(k for occs in ref.occurrences for o in occs for (k, _) in o.fields)
    @info "Occurrence columns published with the hits: " *
          (isempty(columns) ? "none" : join(columns, ", "))
    async && return server
    _wait(listener)
    return server
end

# Blocks until Ctrl-C (also in scripts, where Julia otherwise exits at once on
# SIGINT without closing the server) or until the server is closed.
function _wait(listener)
    interactive = isinteractive()
    interactive || Base.exit_on_sigint(false)
    try
        while isopen(listener)
            sleep(0.5)
        end
    catch err
        err isa InterruptException || rethrow()
    finally
        close(listener)
        interactive || Base.exit_on_sigint(true)
    end
    return nothing
end

function _check_port(ip::AbstractString, port::Int)
    sock = try
        Sockets.listen(parse(Sockets.IPAddr, ip), port)
    catch err
        err isa Base.IOError || rethrow()
        error("cannot listen on $ip:$port ($(err.msg)): is another server (or an earlier " *
              "`serve`) using the port? Stop it or choose another `port`.")
    end
    close(sock)
    return nothing
end

"""
    shutdown(server)

Stop a server started with `serve(...; async = true)`.
"""
shutdown(server::Genie.Server.ServersCollection) = (Genie.Server.down(server); nothing)
