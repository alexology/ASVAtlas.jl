# Limits of the open web service (no accounts): per-client request rates and
# requests in progress, and a bound on the searches and networks waiting for
# the search lock, checked before the request body is read so that waiting
# requests cannot pile up bodies in memory. A client is its IP address (the
# last X-Forwarded-For entry behind a trusted reverse proxy).

# The only requests with a body; everything else is refused before reading one.
const _HEAVY = Dict(("POST", "/api/search") => :search, ("POST", "/api/network") => :network)

mutable struct _Guard
    lock::ReentrantLock
    events::Dict{Tuple{String,Symbol},Vector{Float64}}
    inflight::Dict{String,Int}
    heavy::Int                       # searches and networks running or waiting
    trust_forwarded::Bool
    window::Float64
    rates::Dict{Symbol,Int}          # per client and window (0: no limit)
    max_inflight_per_client::Int     # 0: no limit
    max_heavy::Int                   # running + waiting (0: no limit)
    body_timeout::Float64
end

_Guard(; trust_forwarded = false, window = 600, rate_searches = 60, rate_networks = 30,
       max_inflight_per_client = 2, max_queued = 8, body_timeout = 60) =
    _Guard(ReentrantLock(), Dict{Tuple{String,Symbol},Vector{Float64}}(), Dict{String,Int}(), 0,
           trust_forwarded, Float64(window), Dict(:search => rate_searches, :network => rate_networks),
           max_inflight_per_client, max_queued > 0 ? max_queued + 1 : 0, Float64(body_timeout))

function _client_address(g::_Guard, stream::HTTP.Stream)
    if g.trust_forwarded
        forwarded = HTTP.header(stream.message, "X-Forwarded-For")
        if !isempty(forwarded)
            last_hop = strip(last(split(forwarded, ',')))
            isempty(last_hop) || return String(first(last_hop, 64))
        end
    end
    return string(stream.stream.peerip)
end

"""Admit a search or network of `client`, or throw 429/503; release with `_release!`."""
function _admit!(g::_Guard, client::String, kind::Symbol; now::Float64 = time())
    lock(g.lock) do
        if length(g.events) > 10_000                 # forget idle clients now and then
            for (k, ts) in collect(g.events)
                filter!(t -> now - t < g.window, ts)
                isempty(ts) && delete!(g.events, k)
            end
        end
        times = get!(Vector{Float64}, g.events, (client, kind))
        filter!(t -> now - t < g.window, times)
        rate = get(g.rates, kind, 0)
        if rate > 0 && length(times) >= rate
            wait = max(1, ceil(Int, g.window - (now - first(times))))
            throw(_RequestError(429, "too many requests from this address; retry in $(cld(wait, 60)) minutes", wait))
        end
        g.max_inflight_per_client > 0 && get(g.inflight, client, 0) >= g.max_inflight_per_client &&
            throw(_RequestError(429, "requests from this address are already running; wait for them", 5))
        g.max_heavy > 0 && g.heavy >= g.max_heavy &&
            throw(_RequestError(503, "the server is busy; retry in a moment", 10))
        push!(times, now)
        g.inflight[client] = get(g.inflight, client, 0) + 1
        g.heavy += 1
    end
    return client
end

function _release!(g::_Guard, client::String)
    lock(g.lock) do
        n = get(g.inflight, client, 0) - 1
        n <= 0 ? delete!(g.inflight, client) : (g.inflight[client] = n)
        g.heavy = max(0, g.heavy - 1)
    end
    return nothing
end

# After `seconds` the connection is closed, which ends a pending read: slow or
# stalled uploads cannot hold the server's connections.
function _with_deadline(f, stream::HTTP.Stream, seconds::Real)
    expired = Threads.Atomic{Bool}(false)
    timer = Timer(seconds) do _
        expired[] = true
        try
            close(stream.stream)
        catch
        end
    end
    try
        return f()
    catch err
        expired[] && throw(_RequestError(408, "the request body took too long to arrive"))
        rethrow()
    finally
        close(timer)
    end
end

const _LOOPBACK = ("127.0.0.1", "localhost", "::1")

# Accepted Host names: the allow-list; without one, loopback names only on a
# loopback address (no DNS rebinding from a web page), otherwise any (warning).
function _accepted_hosts(allowed_hosts, host)
    isempty(allowed_hosts) || return Set{String}(lowercase(strip(String(h))) for h in allowed_hosts)
    host in _LOOPBACK && return Set{String}(_LOOPBACK)
    @warn "No allowed_hosts: requests with any Host header are accepted; set it to the public domain."
    return Set{String}()
end
