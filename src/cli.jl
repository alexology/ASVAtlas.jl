# Command-line entry point (bin/asvatlas.jl), configurable with arguments or
# environment variables, for running the application as a service.

const _CLI_OPTIONS = (
    # option, environment variable, help
    ("fasta", "ASVATLAS_FASTA", "reference ASVs (FASTA), required"),
    ("metadata", "ASVATLAS_METADATA", "taxonomy and occurrences (TSV, or CSV for .csv)"),
    ("fields", "ASVATLAS_FIELDS",
     "comma-separated occurrence columns to publish, or none (default: all)"),
    ("coordinate-digits", "ASVATLAS_COORDINATE_DIGITS",
     "round published coordinates to this many decimals (default: as read)"),
    ("name", "ASVATLAS_NAME", "dataset name shown in the page (default: FASTA file name)"),
    ("marker", "ASVATLAS_MARKER", "marker of the reference ASVs shown in the page (e.g. COI-5P)"),
    ("primers", "ASVATLAS_PRIMERS", "primer pair of the reference ASVs shown in the page (e.g. BF3-BR2)"),
    ("about", "ASVATLAS_ABOUT", "JSON description of the reference (tab \"Reference data\"; public)"),
    ("host", "ASVATLAS_HOST", "address to listen on (default 127.0.0.1; 0.0.0.0 for all)"),
    ("port", "ASVATLAS_PORT", "port (default 8080)"),
    ("max-queries", "ASVATLAS_MAX_QUERIES", "most queries per search (default 10000)"),
    ("max-hits", "ASVATLAS_MAX_HITS", "most hits per query (default 500)"),
    ("max-length", "ASVATLAS_MAX_LENGTH", "longest query in bases (default 20000)"),
    ("min-identity", "ASVATLAS_MIN_IDENTITY", "lowest similarity (%) accepted by the page (default 90)"),
    ("max-network-sequences", "ASVATLAS_MAX_NETWORK_SEQUENCES",
     "most sequences in a haplotype network (default 300)"),
    ("allowed-hosts", "ASVATLAS_ALLOWED_HOSTS",
     "comma-separated host names accepted in requests (default: any)"),
    ("example", "ASVATLAS_EXAMPLE",
     "FASTA file offered by \"Try an example\" (public; default: no example)"),
)

const _CLI_HELP = """
    Usage: julia -t auto bin/asvatlas.jl --fasta FILE [options]

    Starts the ASVAtlas web application. Each option can also be set with the
    environment variable shown; command-line values take precedence.

    """ * join(("  --$(rpad(o, 17)) $(rpad(e, 26)) $h" for (o, e, h) in _CLI_OPTIONS), "\n") * """

      --help                                         show this message

    The process runs until it is stopped (Ctrl-C, or SIGTERM from systemd/Docker).
    """

function _parse_cli(args::AbstractVector{<:AbstractString}, env = ENV)
    opts = Dict{String,String}()
    for (o, e, _) in _CLI_OPTIONS
        v = strip(get(env, e, ""))
        isempty(v) || (opts[o] = v)
    end
    known = Set(first.(_CLI_OPTIONS))
    i = 1
    while i <= length(args)
        a = args[i]
        a in ("-h", "--help") && return nothing
        startswith(a, "--") || throw(ArgumentError("unexpected argument \"$a\""))
        key, value = occursin('=', a) ? split(a[3:end], '='; limit = 2) :
                     (a[3:end], i < length(args) ? args[i += 1] : "")
        key in known || throw(ArgumentError("unknown option --$key"))
        isempty(value) && throw(ArgumentError("option --$key needs a value"))
        opts[key] = value
        i += 1
    end
    haskey(opts, "fasta") ||
        throw(ArgumentError("the reference FASTA is required (--fasta or ASVATLAS_FASTA)"))
    function int(key, default)
        haskey(opts, key) || return default
        v = tryparse(Int, opts[key])
        v === nothing && throw(ArgumentError("--$key must be an integer, got \"$(opts[key])\""))
        return v
    end
    fields = if !haskey(opts, "fields")
        nothing
    elseif lowercase(strip(opts["fields"])) == "none"
        String[]
    else
        String[strip(f) for f in split(opts["fields"], ',') if !isempty(strip(f))]
    end
    return (fasta = opts["fasta"], metadata = get(opts, "metadata", nothing), fields,
            coordinate_digits = int("coordinate-digits", nothing),
            name = get(opts, "name", basename(opts["fasta"])),
            marker = get(opts, "marker", ""), primers = get(opts, "primers", ""),
            about = get(opts, "about", nothing),
            host = get(opts, "host", "127.0.0.1"), port = int("port", 8080),
            max_queries = int("max-queries", 10_000), max_hits = int("max-hits", 500),
            max_length = int("max-length", 20_000),
            max_network_sequences = int("max-network-sequences", 300),
            min_identity = int("min-identity", 90),
            example = get(opts, "example", nothing),
            allowed_hosts = String[strip(h) for h in split(get(opts, "allowed-hosts", ""), ',')
                                   if !isempty(strip(h))])
end

# Writes each log record at once: stderr redirected to a file is buffered
# until exit, which would leave the log of a running server empty.
struct _FlushingLogger <: AbstractLogger
    logger::ConsoleLogger
end
Logging.min_enabled_level(l::_FlushingLogger) = Logging.min_enabled_level(l.logger)
Logging.shouldlog(l::_FlushingLogger, args...) = Logging.shouldlog(l.logger, args...)
Logging.catch_exceptions(l::_FlushingLogger) = Logging.catch_exceptions(l.logger)
function Logging.handle_message(l::_FlushingLogger, args...; kwargs...)
    Logging.handle_message(l.logger, args...; kwargs...)
    flush(stderr)
end

"""
    main(args = ARGS) -> Int

Command-line entry point used by `bin/asvatlas.jl`: loads the reference data
and runs [`serve`](@ref) until the process is stopped. Returns the exit code
(1 for invalid options or data).
"""
function main(args = ARGS)
    o = try
        _parse_cli(args)
    catch err
        err isa ArgumentError || rethrow()
        println(stderr, "asvatlas: ", err.msg, "\n\n", _CLI_HELP)
        return 1
    end
    if o === nothing
        print(_CLI_HELP)
        return 0
    end
    global_logger(_FlushingLogger(ConsoleLogger(stderr)))
    ref = try
        load_reference(o.fasta; metadata = o.metadata, name = o.name, marker = o.marker,
                       primers = o.primers, about = o.about, fields = o.fields,
                       coordinate_digits = o.coordinate_digits)
    catch err
        err isa Union{ArgumentError,SystemError} || rethrow()
        println(stderr, "asvatlas: cannot load the reference data: ", sprint(showerror, err))
        return 1
    end
    example = try
        o.example === nothing ? "" : read(o.example, String)
    catch err
        err isa SystemError || rethrow()
        println(stderr, "asvatlas: cannot read the example: ", sprint(showerror, err))
        return 1
    end
    serve(ref; example, host = o.host, port = o.port, max_queries = o.max_queries, max_hits = o.max_hits,
          max_length = o.max_length, allowed_hosts = o.allowed_hosts,
          max_network_sequences = o.max_network_sequences, min_identity = o.min_identity)
    return 0
end
