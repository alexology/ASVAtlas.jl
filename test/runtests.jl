using ASVAtlas
import HaplotypeNetworks as HN
using HTTP
using JSON3
using Random
using Sockets
using Test

const EXAMPLE = joinpath(@__DIR__, "..", "examples", "data")

rng = Xoshiro(7)
random_seq(n) = String([rand(rng, "ACGT") for _ in 1:n])
function mutate(s, n)
    c = collect(s)
    for i in randperm(rng, length(c))[1:n]
        c[i] = rand(rng, filter(!=(c[i]), collect("ACGT")))
    end
    String(c)
end
revcomp(s) = String([Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')[b] for b in reverse(s)])

function write_file(dir, name, text)
    path = joinpath(dir, name)
    write(path, text)
    path
end

@testset "ASVAtlas" begin

@testset "Delimited files" begin
    mktempdir() do dir
        p = write_file(dir, "a.csv", "﻿id, name ,note\r\nA,\"x, y\",\"say \"\"hi\"\"\"\r\n\r\nB,z,\"two\nlines\"\n")
        header, rows = ASVAtlas.read_delimited(p)
        @test header == ["id", "name", "note"]
        @test rows == [["A", "x, y", "say \"hi\""], ["B", "z", "two\nlines"]]
        p = write_file(dir, "b.tsv", "id\tv\nA\t1\t2\n")
        @test_throws ArgumentError ASVAtlas.read_delimited(p)
        p = write_file(dir, "c.csv", "id,v\nA,\"open\n")
        @test_throws ArgumentError ASVAtlas.read_delimited(p)
        p = write_file(dir, "d.tsv", "id\tid\nA\tB\n")
        @test_throws ArgumentError ASVAtlas.read_delimited(p)
        p = write_file(dir, "e.txt", "id;v\nA;1\n")
        @test ASVAtlas.read_delimited(p; delim = ';')[2] == [["A", "1"]]
    end
end

@testset "SINTAX headers" begin
    id, tax = ASVAtlas._parse_header("R1;tax=k:Animalia,g:Exemplus,s:Exemplus_primus;",
                                     "R1;tax=k:Animalia,g:Exemplus,s:Exemplus_primus; extra")
    @test id == "R1"
    @test tax == Dict("kingdom" => "Animalia", "genus" => "Exemplus", "species" => "Exemplus primus")
    @test ASVAtlas._parse_header("R2", "R2 some description") == ("R2", Dict{String,String}())
end

@testset "Reference datasets" begin
    ref = load_reference(joinpath(EXAMPLE, "reference.fasta");
                         metadata = joinpath(EXAMPLE, "reference_metadata.tsv"))
    @test length(ref) == 14
    @test ref.ranks == ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    @test sum(length, ref.occurrences) == 17
    @test count(isempty, ref.occurrences) == 5
    @test taxon(ref, 1) == "Fictiferus alpinus"
    o = first(ref.occurrences[1])
    @test o.latitude ≈ 46.34897 && o.longitude ≈ 10.70302
    @test o.fields == ["site" => "Site A01", "habitat" => "lake", "date" => "2025-08-26"]
    @test occursin("14 references", sprint(show, ref))

    mktempdir() do dir
        s1, s2, s3 = random_seq(300), random_seq(300), random_seq(300)
        fa = write_file(dir, "r.fasta", """
            >R1;tax=k:Animalia,f:Exemplidae,g:Exemplus,s:Exemplus_primus;
            $(lowercase(s1[1:150]))
            $(s1[151:end])
            >R2;tax=k:Animalia,g:Exemplus;
            $(replace(s2, 'T' => 'U'))
            >R3
            $s3
            """)
        # Metadata: Darwin Core coordinates, quoted CSV, a reference repeated,
        # one without coordinates, one row of an unknown identifier.
        meta = write_file(dir, "m.csv", """
            ASV_ID,decimalLatitude,decimalLongitude,Species,locality
            R1,45.1,9.2,Exemplus secundus,"Lake, north shore"
            R1,-10.5,-60.25,,
            R2,,,Exemplus tertius,somewhere
            R9,1,1,,
            """)
        ref = @test_logs (:warn, r"1 identifier") load_reference(fa; metadata = meta, name = "test")
        @test ref.name == "test"
        @test ref.ids == ["R1", "R2", "R3"]
        @test ref.lengths == [300, 300, 300]
        @test ref.ranks == ["kingdom", "family", "genus", "species"]
        # Metadata take precedence over the header, rank by rank.
        @test ref.lineages[1] == ["Animalia", "Exemplidae", "Exemplus", "Exemplus secundus"]
        @test ref.lineages[2] == ["Animalia", "", "Exemplus", "Exemplus tertius"]
        @test ref.lineages[3] == ["", "", "", ""]
        @test taxon(ref, 3) == ""
        @test [(o.latitude, o.longitude) for o in ref.occurrences[1]] == [(45.1, 9.2), (-10.5, -60.25)]
        @test ref.occurrences[1][1].fields == ["locality" => "Lake, north shore"]
        @test isempty(ref.occurrences[1][2].fields)
        @test isempty(ref.occurrences[2])

        # Published occurrence columns and coordinate precision
        meta3 = write_file(dir, "m3.tsv", """
            id\tlat\tlon\tsite\tcollector\tnotes
            R1\t45.12345\t9.98765\tS1\tJane Doe\tprivate
            R2\t-10.55\t-60.04\tS2\t\t
            """)
        ref3 = load_reference(fa; metadata = meta3, fields = ["SITE"], coordinate_digits = 1)
        @test ref3.occurrences[1][1].fields == ["site" => "S1"]
        @test ref3.occurrences[2][1].fields == ["site" => "S2"]
        @test (ref3.occurrences[1][1].latitude, ref3.occurrences[1][1].longitude) == (45.1, 10.0)
        @test (ref3.occurrences[2][1].latitude, ref3.occurrences[2][1].longitude) == (-10.6, -60.0)
        ref3 = load_reference(fa; metadata = meta3, fields = String[])
        @test all(isempty(o.fields) for occs in ref3.occurrences for o in occs)
        @test ref3.occurrences[1][1].latitude == 45.12345
        ref3 = load_reference(fa; metadata = meta3)
        @test first.(ref3.occurrences[1][1].fields) == ["site", "collector", "notes"]
        @test_throws ArgumentError load_reference(fa; metadata = meta3, fields = ["email"])
        @test_throws ArgumentError load_reference(fa; metadata = meta3, fields = ["lat"])
        @test_throws ArgumentError load_reference(fa; metadata = meta3, fields = "site")
        @test_throws ArgumentError load_reference(fa; metadata = meta3, coordinate_digits = -1)
        @test_throws ArgumentError load_reference(fa; fields = ["site"])

        # Errors
        bad(text) = load_reference(fa; metadata = write_file(dir, "bad.tsv", text))
        @test_throws ArgumentError bad("id\tspecies\nR1\tA b\nR1\tC d\n")
        @test_throws ArgumentError bad("id\tlat\tlon\nR1\t95\t0\n")
        @test_throws ArgumentError bad("id\tlat\tlon\nR1\t1\t181\n")
        @test_throws ArgumentError bad("id\tlat\tlon\nR1\tnorth\t0\n")
        @test_throws ArgumentError bad("id\tlat\tlon\nR1\t1\t\n")
        @test_throws ArgumentError bad("id\tlat\nR1\t1\n")
        @test_throws ArgumentError bad("name\tlat\tlon\nR1\t1\t1\n")
        @test_throws ArgumentError bad("id\tgenus\tGenus\nR1\tA\tA\n")
        @test_throws ArgumentError bad("\tlat\tlon\n\t1\t1\n")
        ref2 = load_reference(fa; metadata = write_file(dir, "ok.tsv", "name\ty\tx\nR3\t1\t2\n"),
                              id_column = "name", latitude_column = "y", longitude_column = "x")
        @test ref2.occurrences[3][1].latitude == 1 && ref2.occurrences[3][1].longitude == 2
        @test ref2.ranks == ["kingdom", "family", "genus", "species"]
        @test_throws ArgumentError load_reference(fa; metadata = meta, id_column = "nope")
        @test_throws ArgumentError load_reference(write_file(dir, "dup.fa", ">A\nACGT\n>A\nACGG\n"))
        @test_throws ArgumentError load_reference(write_file(dir, "inv.fa", ">A\nACGTZ\n"))
        @test_throws ArgumentError load_reference(write_file(dir, "empty.fa", ""))
    end
end

@testset "Queries" begin
    q = parse_queries(">a desc\nacg-u\nTT\n\n>b\nGGGG\r\n")
    @test q == ["a" => "ACGTTT", "b" => "GGGG"]
    @test parse_queries("ACGT\n\n  tgca  \n") == ["query_1" => "ACGT", "query_2" => "TGCA"]
    @test_throws ArgumentError parse_queries("")
    @test_throws ArgumentError parse_queries(">a\nACGT\n>a\nACGT\n")
    @test_throws ArgumentError parse_queries(">a\nACGTX\n")
    @test_throws ArgumentError parse_queries("ACGT\n>a\nACGT\n")
    mktempdir() do dir
        @test read_queries(write_file(dir, "q.fa", ">x\nAC\n")) == ["x" => "AC"]
    end
end

@testset "Search" begin
    base = random_seq(450)
    other = random_seq(450)
    seqs = [base, mutate(base, 4), mutate(base, 20), mutate(base, 30), other, other]
    taxa = ["Exemplus primus", "Exemplus primus", "Exemplus secundus", "Exemplus", "Alius unus",
            "Alius duo"]
    mktempdir() do dir
        header(t) = occursin(' ', t) ? "g:$(split(t)[1]),s:$(replace(t, ' ' => '_'))" : "g:$t"
        fa = write_file(dir, "r.fa", join((">R$i;tax=$(header(t));\n$s\n"
                                           for (i, (s, t)) in enumerate(zip(seqs, taxa)))))
        meta = write_file(dir, "m.tsv", "id\tlat\tlon\tsite\nR1\t45\t9\tS1\nR1\t46\t10\tS2\nR3\t40\t15\tS3\n")
        ref = load_reference(fa; metadata = meta)
        @test ref.ranks == ["genus", "species"]
        @test ref.lineages[4] == ["Exemplus", ""]
        queries = ["q_exact" => base, "q_rc" => revcomp(seqs[3]), "q_none" => random_seq(450),
                   "q_tie" => other, "q_part" => base[51:400]]
        r100 = search_asvs(ref, queries; identity = 100, threads = 1)
        @test [x.id for x in r100] == first.(queries)
        @test [ref.ids[h.reference] for h in r100[1].hits] == ["R1"]
        @test only(r100[1].hits).identity == 100 && only(r100[1].hits).coverage == 1
        @test only(r100[2].hits).strand == '-'
        @test ref.ids[only(r100[2].hits).reference] == "R3"
        @test isempty(r100[3].hits)
        @test Set(ref.ids[h.reference] for h in r100[4].hits) == Set(["R5", "R6"])
        @test Set(best_taxa(ref, r100[4])) == Set(["Alius unus", "Alius duo"])
        h = only(r100[5].hits)
        @test (h.query_start, h.query_end, h.target_start, h.target_end) == (1, 350, 51, 400)
        @test h.cigar == "350M"
        @test best_taxa(ref, r100[3]) == String[]

        r97 = search_asvs(ref, queries; identity = 97)
        @test [ref.ids[h.reference] for h in r97[1].hits] == ["R1", "R2"]
        r90 = search_asvs(ref, queries; identity = 90, threads = 4)
        @test [ref.ids[h.reference] for h in r90[1].hits] == ["R1", "R2", "R3", "R4"]
        @test issorted([h.identity for h in r90[1].hits]; rev = true)
        @test r90[1].hits[2].identity ≈ 100 * 446 / 450
        for (lo, hi) in ((r97, r100), (r90, r97))
            for (a, b) in zip(lo, hi)
                @test Set(h.reference for h in b.hits) ⊆ Set(h.reference for h in a.hits)
            end
        end
        for (t, r) in ((100, r100), (97, r97), (90, r90)), x in r, h in x.hits
            @test h.identity >= t - 1e-9
        end
        @test length(search_asvs(ref, queries; identity = 90, maxhits = 2)[1].hits) == 2
        @test !isempty(search_asvs(ref, ["q" => seqs[3]]; identity = 100, strand = :plus)[1].hits)
        @test isempty(search_asvs(ref, ["q" => revcomp(seqs[3])]; identity = 100, strand = :plus)[1].hits)
        @test search_asvs(ref, [base])[1].id == "query_1"
        @test isempty(search_asvs(ref, Pair{String,String}[]))
        @test_throws ArgumentError search_asvs(ref, queries; identity = 0)
        @test_throws ArgumentError search_asvs(ref, queries; identity = 101)
        @test_throws ArgumentError search_asvs(ref, queries; query_cov = 0)
        @test_throws ArgumentError search_asvs(ref, queries; maxhits = 0)
        @test_throws ArgumentError search_asvs(ref, ["a" => base, "a" => base])

        @testset "Tables" begin
            s = summary_table(ref, r90)
            @test s.query == first.(queries)
            @test s.n_hits == length.(getfield.(r90, :hits))
            @test s.best_identity[1] == 100 && ismissing(s.best_identity[3])
            @test s.best_references[1] == "R1"
            @test split(s.best_references[4], ';') ⊆ ["R5", "R6"]
            @test s.n_georeferenced[1] == 2 && s.n_occurrences[1] == 3
            h = hit_table(ref, r90)
            @test length(h.query) == sum(s.n_hits)
            @test h.reference[1:4] == ["R1", "R2", "R3", "R4"] && h.hit[1:4] == 1:4
            @test h.species[1] == "Exemplus primus" && h.genus[4] == "Exemplus" && h.species[4] == ""
            @test h.taxon[4] == "Exemplus"
            o = occurrence_table(ref, r90)
            @test o.query[1:3] == ["q_exact", "q_exact", "q_exact"]
            @test o.fields[1] == "site=S1"
            @test length(o.query) == sum(s.n_occurrences)
            path = joinpath(dir, "hits.tsv")
            write_table(path, h)
            header, rows = ASVAtlas.read_delimited(path)
            @test header == string.(collect(keys(h)))
            @test length(rows) == length(h.query)
            @test rows[1][4] == "100.0"
            @test occursin("best_identity", sprint(write_table, s))
            @test sprint(write_table, (a = ["x\ty"], b = [missing])) == "a\tb\nx y\t\n"
        end
    end
end

@testset "Web application" begin
    ref = load_reference(joinpath(EXAMPLE, "reference.fasta");
                         metadata = joinpath(EXAMPLE, "reference_metadata.tsv"), name = "Example")
    port, sock = listenany(ip"127.0.0.1", 20000)
    close(sock)
    server = @test_logs (:info, r"ASVAtlas") (:info, r"Warm-up") match_mode = :any serve(
        ref; port, async = true, max_queries = 10, max_hits = 20, max_body = 100_000)
    url = "http://127.0.0.1:$port"
    get(p) = HTTP.get(url * p; status_exception = false)
    post(body) = HTTP.post(url * "/api/search", ["Content-Type" => "application/json"], body;
                           status_exception = false)
    try
        r = get("/")
        @test r.status == 200
        @test occursin("text/html", HTTP.header(r, "Content-Type"))
        @test occursin("\"leaflet/leaflet.js\"", String(r.body))
        @test !occursin("unpkg", String(r.body))
        # Relative URLs only, so that a reverse proxy can use a sub-path.
        @test !occursin(r"(href|src)=\"/", String(r.body))
        @test !occursin("fetch(\"/", String(get("/app.js").body))
        h = JSON3.read(get("/health").body)
        @test h.status == "ok" && h.references == 14
        @test !isempty(HTTP.header(r, "Content-Security-Policy"))
        for p in ("/app.js", "/style.css", "/leaflet/leaflet.js", "/leaflet/leaflet.css",
                  "/leaflet/images/layers.png")
            @test get(p).status == 200
        end
        # Optional sampling locations: a table read in the page, never sent.
        @test occursin("id=\"coord-file\"", String(get("/").body))
        @test occursin("function readCoordTable", String(get("/app.js").body))
        # Tab "Reference data".
        @test occursin("id=\"tab-about\"", String(get("/").body)) && get("/about.js").status == 200
        # Powered by dnambientale.it, with the group's icon.
        @test occursin("href=\"https://dnambientale.it/en/index.html\"", String(get("/").body))
        @test get("/dnambientale.png").body == read(joinpath(@__DIR__, "..", "src", "assets", "dnambientale.png"))
        @test HTTP.header(get("/dnambientale.png"), "Content-Type") == "image/png"
        @test get("/leaflet/leaflet.js").body ==
              read(joinpath(@__DIR__, "..", "src", "assets", "leaflet", "leaflet.js"))
        world = JSON3.read(get("/world.geojson").body)
        @test world.type == "FeatureCollection" && length(world.features) > 200
        @test any(f -> f.properties.name == "Italy", world.features)
        @test get("/nothing").status == 404
        for p in ("/../Project.toml", "/%2e%2e/Project.toml", "/..%2fProject.toml")
            @test get(p).status != 200
        end
        @test HTTP.post(url * "/"; status_exception = false).status == 405

        d = JSON3.read(get("/api/dataset").body)
        @test d.name == "Example" && d.references == 14 && d.occurrences == 17
        @test d.georeferenced == 9 && d.max_queries == 10 && d.max_hits == 20
        @test collect(d.identity_levels) == [100, 99, 98, 97, 95, 90] && d.min_identity == 90

        fasta = read(joinpath(EXAMPLE, "queries.fasta"), String)
        r = post(JSON3.write((sequences = fasta, identity = 97, query_cov = 0.9, max_hits = 10)))
        @test r.status == 200
        res = JSON3.read(r.body)
        @test res.params.identity == 97 && res.params.strand == "both"
        @test collect(res.ranks) == ref.ranks
        expected = search_asvs(ref, read_queries(joinpath(EXAMPLE, "queries.fasta"));
                               identity = 97, maxhits = 10)
        @test [q.id for q in res.queries] == [x.id for x in expected]
        @test [length(q.hits) for q in res.queries] == [length(x.hits) for x in expected]
        q1 = res.queries[1]
        @test q1.hits[1].reference == "REF_001" && q1.hits[1].identity == 100
        @test collect(q1.best_taxa) == ["Fictiferus alpinus"]
        r1 = res.references[:REF_001]
        @test r1.taxon == "Fictiferus alpinus" && length(r1.occurrences) == 2
        @test r1.occurrences[1].fields[1].name == "site"
        hit_ids = Set(h.reference for q in res.queries for h in q.hits)
        @test Set(String.(keys(res.references))) == hit_ids
        @test res.queries[4].hits[1].strand == "-"

        # Each search is logged (checked on the request function: the server
        # logs from its own tasks).
        @test_logs (:info, r"Search: 6 queries at ≥ 100% identity, 2 with hits") ASVAtlas._search_request(
            ref, JSON3.write((sequences = fasta, identity = 100)), ReentrantLock();
            max_queries = 10, max_hits = 20, threads = 1)
        r = post(JSON3.write((sequences = fasta, identity = 100)))
        @test sum(length(q.hits) for q in JSON3.read(r.body).queries) == 2
        r = post(JSON3.write((sequences = fasta, strand = "plus", identity = 100)))
        @test isempty(JSON3.read(r.body).queries[4].hits)

        for body in ("not json", "[1]", JSON3.write((sequences = 5,)), JSON3.write((sequences = "",)),
                     JSON3.write((sequences = fasta, identity = 120)),
                     JSON3.write((sequences = fasta, identity = "97")),
                     JSON3.write((sequences = fasta, max_hits = 21)),
                     JSON3.write((sequences = fasta, strand = "minus")),
                     JSON3.write((sequences = ">a\nACGTJ\n",)),
                     JSON3.write((sequences = join((">q$i\nACGT" for i in 1:11), '\n'),)))
            r = post(body)
            @test r.status == 400
            @test haskey(JSON3.read(r.body), :error)
        end
        @test post(JSON3.write((sequences = "A"^200_000,))).status == 413
        @test HTTP.get(url * "/api/search"; status_exception = false).status == 405
        # The port is taken: a clear error instead of a silent return.
        err = try
            serve(ref; port, async = true)
            nothing
        catch e
            e
        end
        @test err isa ErrorException && occursin("cannot listen", err.msg)
        @test get("/api/dataset").status == 200
    finally
        shutdown(server)
    end
end

@testset "Web security" begin
    ref = load_reference(joinpath(EXAMPLE, "reference.fasta");
                         metadata = joinpath(EXAMPLE, "reference_metadata.tsv"), name = "Example")
    port, sock = listenany(ip"127.0.0.1", 21000)
    close(sock)
    server = serve(ref; port, async = true, warmup = false, max_queries = 10, max_hits = 20,
                   max_length = 1000, max_body = 100_000,
                   allowed_hosts = ["asvatlas.example.org", "127.0.0.1"])
    url = "http://127.0.0.1:$port"
    get(p, h = Pair{String,String}[]) = HTTP.get(url * p, h; status_exception = false)
    post(body; type = "application/json", h = Pair{String,String}[]) =
        HTTP.post(url * "/api/search", ["Content-Type" => type, h...], body; status_exception = false)
    # Status line of a request written by hand on a socket.
    function raw_status(text)
        s = connect(ip"127.0.0.1", port)
        try
            write(s, text)
            return readline(s)
        finally
            close(s)
        end
    end
    fasta = read(joinpath(EXAMPLE, "queries.fasta"), String)
    try
        # Security headers on pages, files, API answers and errors.
        for r in (get("/"), get("/app.js"), get("/world.geojson"), get("/api/dataset"), get("/nothing"))
            @test occursin("default-src 'self'", HTTP.header(r, "Content-Security-Policy"))
            @test HTTP.header(r, "X-Content-Type-Options") == "nosniff"
            @test HTTP.header(r, "X-Frame-Options") == "DENY"
            @test HTTP.header(r, "Referrer-Policy") == "strict-origin"
        end
        @test HTTP.header(get("/world.geojson"), "Content-Type") == "application/geo+json"
        # Genie's error pages echo the path with the content type of the request:
        # replaced by a fixed JSON message.
        r = get("/%3Cscript%3Ealert(1)%3C/script%3E", ["Accept" => "text/html"])
        @test r.status == 404
        @test startswith(HTTP.header(r, "Content-Type"), "application/json")
        body = String(copy(r.body))
        @test !occursin("script", body) && JSON3.read(body).error == "not found"
        # Only the listed files are served.
        for p in ("/../Project.toml", "/..%2f..%2fProject.toml", "/..%5c..%5cProject.toml",
                  "/server.jl", "/assets/app.js", "/leaflet/LICENSE", "/leaflet")
            @test get(p).status == 404
        end
        # Host header (DNS rebinding).
        @test get("/health", ["Host" => "evil.example"]).status == 421
        @test get("/health", ["Host" => "asvatlas.example.org:443"]).status == 200
        @test get("/health", ["Host" => "ASVATLAS.example.org"]).status == 200
        # JSON only (no cross-site form posts).
        @test post(JSON3.write((sequences = fasta,)); type = "text/plain").status == 415
        @test post(JSON3.write((sequences = fasta,)); type = "application/x-www-form-urlencoded").status == 415
        @test post(JSON3.write((sequences = fasta,)); type = "application/json; charset=utf-8").status == 200
        # Limits on the work of one request.
        @test post(JSON3.write((sequences = ">long\n" * "ACGT"^251,))).status == 400
        @test post(JSON3.write((sequences = fasta, identity = 40))).status == 400
        @test post(JSON3.write((sequences = fasta, identity = true))).status == 400
        @test post(JSON3.write((sequences = fasta, max_hits = true))).status == 400
        @test_throws ArgumentError ASVAtlas._search_request(
            ref, JSON3.write((sequences = fasta, max_hits = 200_000)), ReentrantLock();
            max_queries = 10, max_hits = 1_000_000, threads = 1)
        # Bodies over the limit are refused before being read: by Content-Length...
        @test raw_status("POST /api/search HTTP/1.1\r\nHost: 127.0.0.1\r\n" *
                         "Content-Type: application/json\r\nContent-Length: 999999999\r\n\r\n{")[1:12] ==
              "HTTP/1.1 413"
        @test raw_status("POST /api/search HTTP/1.1\r\nHost: 127.0.0.1\r\n" *
                         "Content-Type: application/json\r\nContent-Length: -5\r\n\r\n") ==
              "HTTP/1.1 400 Bad Request"
        # ...or while reading a chunked body.
        chunk = "C350\r\n" * "A"^50_000 * "\r\n"
        @test raw_status("POST /api/search HTTP/1.1\r\nHost: 127.0.0.1\r\n" *
                         "Content-Type: application/json\r\nTransfer-Encoding: chunked\r\n\r\n" *
                         chunk^3)[1:12] == "HTTP/1.1 413"
        # The server still answers afterwards.
        @test get("/health").status == 200
        @test post(JSON3.write((sequences = fasta,))).status == 200
    finally
        shutdown(server)
    end
end

@testset "Haplotype networks" begin
    ref = load_reference(joinpath(EXAMPLE, "reference.fasta");
                         metadata = joinpath(EXAMPLE, "reference_metadata.tsv"))
    taxa = taxon_list(ref)
    @test taxa["genus"] == ["Exemplus" => 5, "Fictiferus" => 7, "Pseudofictus" => 2]
    @test taxon_list(ref; min_references = 3)["species"] ==
          ["Exemplus primus" => 3, "Fictiferus alpinus" => 3, "Fictiferus lacustris" => 3]
    @test taxon_references(ref, "Species", "Fictiferus alpinus") == [1, 2, 3]
    @test_throws ArgumentError taxon_references(ref, "subgenus", "X")
    @test reference_sequence(ref, 4) ==
          ASVAtlas._read_reference_fasta(joinpath(EXAMPLE, "reference.fasta"))[2][4]

    for m in NETWORK_METHODS
        tn = taxon_network(ref, "genus", "Fictiferus"; method = m, threads = 1)
        @test tn.method == m && length(tn.ids) == 7
        @test HN.nhaplotypes(tn.network.hapset) == 7
        @test size(tn.layout) == (HN.nnodes(tn.network), 2)
        @test all(isfinite, tn.layout)
        @test occursin("Fictiferus", sprint(show, tn))
        # Haplotypes are numbered by decreasing number of occurrences.
        f = HN.haplotype_frequencies(tn.network.hapset)
        @test issorted(f; rev = true)
        @test f[1] == maximum(max(1, length(ref.occurrences[i])) for i in taxon_references(ref, "genus", "Fictiferus"))
    end
    # MST: one edge fewer than haplotypes, steps = Hamming distances of the alignment.
    tn = taxon_network(ref, "genus", "Fictiferus"; method = :mst)
    @test HN.nedges(tn.network) == 6
    # Query ASVs: identical to a reference -> same haplotype, frequency unchanged.
    q = read_queries(joinpath(EXAMPLE, "queries.fasta"))
    tn = taxon_network(ref, "genus", "Fictiferus"; queries = q[1:2])
    hs = tn.network.hapset
    j = findfirst(==("ASV_identical"), tn.ids)
    k = findfirst(==("REF_001"), tn.ids)
    @test hs.assignment[j] == hs.assignment[k]
    @test tn.reference[j] == 0 && tn.reference[k] == 1
    @test HN.nhaplotypes(hs) == 8    # ASV_99 differs from every reference
    @test HN.haplotype_frequencies(hs)[hs.assignment[findfirst(==("ASV_99"), tn.ids)]] == 0
    # Sequences of different lengths are aligned: the trimmed query joins its reference.
    tn = taxon_network(ref, "species", "Exemplus primus"; queries = ["cut" => reference_sequence(ref, 11)[11:410]])
    hs = tn.network.hapset
    @test hs.assignment[findfirst(==("cut"), tn.ids)] == hs.assignment[findfirst(==("REF_011"), tn.ids)]
    @test length(unique(length.(tn.alignment))) == 1
    # A single sequence.
    tn = taxon_network(ref, "species", "Fictiferus sp.")
    @test HN.nnodes(tn.network) == 1 && HN.nedges(tn.network) == 0
    # Sites: every occurrence of the references counted once, by haplotype.
    tn = taxon_network(ref, "genus", "Fictiferus")
    sites = network_sites(ref, tn)
    refs = taxon_references(ref, "genus", "Fictiferus")
    @test sum(sum(values(s.counts)) for s in sites) == sum(length(ref.occurrences[i]) for i in refs)
    @test sum(length(s.occurrences) for s in sites) == sum(length(ref.occurrences[i]) for i in refs)
    # Errors.
    @test_throws ArgumentError taxon_network(ref, "genus", "Nonexistens")
    @test_throws ArgumentError taxon_network(ref, "genus", "Fictiferus"; method = :nj)
    @test_throws ArgumentError taxon_network(ref, "genus", "Fictiferus"; max_sequences = 5)
    @test_throws ArgumentError taxon_network(ref, "genus", "Fictiferus"; epsilon = -1)

    # Web API: networks only for species matched by the ASVs sent.
    port, sock = listenany(ip"127.0.0.1", 22000)
    close(sock)
    fasta = read(joinpath(EXAMPLE, "queries.fasta"), String)
    server = serve(ref; port, async = true, warmup = true, max_network_sequences = 5, example = fasta)
    url = "http://127.0.0.1:$port"
    post(body; type = "application/json") =
        HTTP.post(url * "/api/network", ["Content-Type" => type], body; status_exception = false)
    try
        d = JSON3.read(HTTP.get(url * "/api/dataset").body)
        # The example offered by the page, chosen by whoever runs the server.
        @test d.example == fasta
        @test collect(d.network_methods) == ["mjn", "tcs", "msn", "mst"] && d.max_network_sequences == 5
        @test d.networks == true && d.min_identity == 90
        @test d.marker == "" && d.primers == ""   # not given for the example
        # The list of taxa is not published.
        @test HTTP.get(url * "/api/taxa"; status_exception = false).status == 404
        r = post(JSON3.write((taxon = "Fictiferus alpinus", sequences = fasta, include = ["ASV_identical"],
                              identity = 100, method = "msn")))
        @test r.status == 200
        n = JSON3.read(r.body)
        @test n.method == "msn" && n.sequences == 4 && length(n.haplotypes) == 3
        @test any(h -> collect(h.queries) == ["ASV_identical"] && collect(h.references) == ["REF_001"], n.haplotypes)
        @test n.haplotypes[1].taxa[1].name == "Fictiferus alpinus" && n.haplotypes[1].taxa[1].rank == "species"
        @test all(e -> e.steps >= 1, n.edges)
        @test sum(c.n for s in n.sites for c in s.counts) ==
              sum(length(ref.occurrences[i]) for i in taxon_references(ref, "species", "Fictiferus alpinus"))
        @test n.sites[1].occurrences[1].fields[1].name == "site"
        # ASV_99 matches its species only from 99% on.
        lac = (taxon = "Fictiferus lacustris", sequences = fasta, include = ["ASV_99"])
        @test post(JSON3.write((; lac..., identity = 100))).status == 400
        r = post(JSON3.write((; lac..., identity = 99)))
        @test r.status == 200 && JSON3.read(r.body).sequences == 4
        # Of the ASVs listed, only those matching the species enter its network.
        r = post(JSON3.write((taxon = "Fictiferus alpinus", sequences = fasta, identity = 99,
                              include = ["ASV_identical", "ASV_99"])))
        n = JSON3.read(r.body)
        @test r.status == 200 && n.sequences == 4 && !any(h -> "ASV_99" in h.queries, n.haplotypes)
        # A query read on the other strand is turned around.
        r = post(JSON3.write((taxon = "Pseudofictus minor", sequences = fasta, include = ["ASV_reverse"],
                              identity = 100)))
        @test r.status == 200 && collect(JSON3.read(r.body).reversed) == ["ASV_reverse"]
        # Errors and limits.
        alp = (taxon = "Fictiferus alpinus", sequences = fasta, include = ["ASV_identical"])
        for body in ((taxon = "Fictiferus alpinus",),                          # no ASV of the user
                     (taxon = "Fictiferus alpinus", sequences = fasta, include = String[]),
                     (taxon = "Exemplus secundus", sequences = fasta, include = ["ASV_identical"]),
                     (taxon = "Nonexistens", sequences = fasta, include = ["ASV_identical"]),
                     (; alp..., identity = 80),                                # below min_identity
                     (; alp..., method = "nj"), (; alp..., epsilon = 11), (; alp..., epsilon = true),
                     (; alp..., taxon = 5),
                     (; alp..., include = ["missing"]), (; alp..., include = "ASV_identical"),
                     (; alp..., include = ["ASV_identical", "ASV_99", "ASV_97", "ASV_reverse",
                                           "ASV_genus", "ASV_unknown"]))   # 6 > 5 ASVs
            r = post(JSON3.write(body))
            @test r.status == 400
            @test haskey(JSON3.read(r.body), :error)
        end
        @test post(JSON3.write(alp); type = "text/plain").status == 415
        @test HTTP.get(url * "/api/network"; status_exception = false).status == 405
        @test HTTP.get(url * "/network.js").status == 200
    finally
        shutdown(server)
    end
end

@testset "Network orientation and numbering (AT1, AT3)" begin
    rng2 = Xoshiro(26)
    a = String(rand(rng2, collect("ACGT"), 120))
    b = a[1:59] * (a[60] == 'A' ? "C" : "A") * a[61:end]         # one SNP
    c = a[1:20] * (a[21] == 'A' ? "G" : "A") * a[22:89] * (a[90] == 'C' ? "T" : "C") * a[91:end]   # two SNPs
    u = String(rand(rng2, collect("ACGT"), 120))                    # unrelated
    mktempdir() do dir
        # AT1: two references one SNP apart; the reverse complement of the
        # first one as query, and a third reference stored on the other strand.
        fa = write_file(dir, "r.fa", ">R1;tax=s:Exemplus_primus;\n$a\n>R2;tax=s:Exemplus_primus;\n$b\n" *
                                     ">R3;tax=s:Exemplus_primus;\n$(revcomp(b))\n")
        ref = load_reference(fa)
        hit = only(search_asvs(ref, ["q" => revcomp(a)]; identity = 100)[1].hits)
        @test hit.strand == '-' && ref.ids[hit.reference] == "R1"
        tn = taxon_network(ref, "species", "Exemplus primus"; queries = ["q" => revcomp(a)])
        hs = tn.network.hapset
        @test length(hs.sites) == 120                     # nothing masked by a wrong strand
        @test HN.nhaplotypes(hs) == 2
        j(id) = findfirst(==(id), tn.ids)
        @test hs.assignment[j("q")] == hs.assignment[j("R1")]
        @test hs.assignment[j("R3")] == hs.assignment[j("R2")]
        @test tn.reversed[j("q")] && tn.reversed[j("R3")] && !tn.reversed[j("R1")] && !tn.reversed[j("R2")]
        @test HN.nedges(tn.network) == 1 && only(tn.network.steps) == 1
        @test Set(ASVAtlas._network_json(ref, tn).reversed) == Set(["q", "R3"])
        # Distances between haplotypes, used by the grouping slider of the page.
        @test ASVAtlas._network_json(ref, tn).distances == [[0, 1], [1, 0]]
        exref = load_reference(joinpath(EXAMPLE, "reference.fasta");
                               metadata = joinpath(EXAMPLE, "reference_metadata.tsv"))
        tn = taxon_network(exref, "genus", "Fictiferus"; method = :mst)
        D = reduce(hcat, ASVAtlas._network_json(exref, tn).distances)
        @test D == D' && all(iszero, [D[i, i] for i in axes(D, 1)])
        @test [D[tn.network.src[e], tn.network.dst[e]] for e in 1:HN.nedges(tn.network)] == tn.network.steps
        # A sequence unrelated to the others on both strands is left out and reported.
        tn = @test_logs (:warn, r"Left out .*: x") taxon_network(ref, "species", "Exemplus primus";
                                                                 queries = ["x" => "A"^60])
        @test tn.excluded == ["x"] && !("x" in tn.ids) && length(tn.network.hapset.sites) == 120
        @test ASVAtlas._network_json(ref, tn).excluded == ["x"]
        # Even when the unrelated sequence comes first (most occurrences): the
        # orientation starts from the sequence linked to most others.
        fa = write_file(dir, "u.fa", ">U1\n$u\n>R1\n$a\n>R2\n$(revcomp(b))\n")
        meta = "id\tspecies\tlatitude\tlongitude\n" * join(("U1\tExemplus primus\t$k\t$k\n" for k in 1:5)) *
               "R1\tExemplus primus\t\t\nR2\tExemplus primus\t\t\n"
        refu = load_reference(fa; metadata = write_file(dir, "u.tsv", meta))
        tn = @test_logs (:warn, r"U1") taxon_network(refu, "species", "Exemplus primus")
        @test tn.excluded == ["U1"] && Set(tn.ids) == Set(["R1", "R2"])
        @test HN.nhaplotypes(tn.network.hapset) == 2 && tn.reversed[findfirst(==("R2"), tn.ids)] !=
                                                          tn.reversed[findfirst(==("R1"), tn.ids)]

        # AT3: one reference with 5 occurrences, two identical ones with 4 each:
        # their haplotype (8) must be H1, before the 5.
        fa = write_file(dir, "n.fa", ">A1\n$c\n>B1\n$a\n>B2\n$a\n")
        rows = ["id\tspecies\tlatitude\tlongitude\tsite"]
        for (id, n, lat) in (("A1", 5, 10.0), ("B1", 4, 20.0), ("B2", 4, 30.0)), k in 1:n
            push!(rows, "$id\tExemplus primus\t$(lat + k)\t$k\t$id-$k")
        end
        ref = load_reference(fa; metadata = write_file(dir, "n.tsv", join(rows, '\n') * "\n"))
        tn = taxon_network(ref, "species", "Exemplus primus")
        hs = tn.network.hapset
        @test HN.haplotype_frequencies(hs) == [8.0, 5.0]
        @test Set(tn.ids[hs.members[1]]) == Set(["B1", "B2"]) && tn.ids[only(hs.members[2])] == "A1"
        @test issorted(HN.haplotype_frequencies(taxon_network(ref, "species", "Exemplus primus";
                                                              method = :tcs).network.hapset); rev = true)
        # Sites, JSON and network follow the same numbering.
        counts = Dict{Int,Int}()
        for s in network_sites(ref, tn), (h, n) in s.counts
            counts[h] = get(counts, h, 0) + n
        end
        @test counts == Dict(1 => 8, 2 => 5)
        js = ASVAtlas._network_json(ref, tn)
        @test [h.name for h in js.haplotypes] == ["H1", "H2"] && [h.occurrences for h in js.haplotypes] == [8, 5]
        @test sort(js.haplotypes[1].references) == ["B1", "B2"]
        @test tn.network.names[1:2] == ["H1", "H2"]
    end
end

@testset "Command line" begin
    noenv = Dict{String,String}()
    o = ASVAtlas._parse_cli(["--fasta", "r.fa", "--metadata=m.tsv", "--port", "9000",
                             "--host", "0.0.0.0", "--max-queries", "5"], noenv)
    @test o == (fasta = "r.fa", metadata = "m.tsv", fields = nothing, coordinate_digits = nothing,
                name = "r.fa", marker = "", primers = "", about = nothing, host = "0.0.0.0", port = 9000,
                max_queries = 5, max_hits = 500,
                max_length = 20_000, max_network_sequences = 300, min_identity = 90,
                example = nothing, allowed_hosts = String[])
    @test ASVAtlas._parse_cli(["--fasta", "r.fa", "--min-identity", "97"], noenv).min_identity == 97
    @test ASVAtlas._parse_cli(["--fasta", "r.fa", "--example", "q.fa"], noenv).example == "q.fa"
    o = ASVAtlas._parse_cli(["--fasta", "r.fa", "--marker", "COI-5P"], Dict("ASVATLAS_PRIMERS" => "BF3-BR2"))
    @test o.marker == "COI-5P" && o.primers == "BF3-BR2"
    r = load_reference(joinpath(EXAMPLE, "reference.fasta"); marker = " COI-5P ", primers = "BF3-BR2")
    @test r.marker == "COI-5P" && r.primers == "BF3-BR2" && occursin("COI-5P, BF3-BR2", sprint(show, r))
    @test JSON3.read(JSON3.write(ASVAtlas._dataset_json(r; max_queries = 1, max_hits = 1, max_length = 1))).marker == "COI-5P"
    @test ASVAtlas._parse_cli(["--fasta", "r.fa", "--about", "a.json"], noenv).about == "a.json"
    # Description of the reference (tab "Reference data") and counts per study.
    mktempdir() do dir
        meta = write_file(dir, "m.tsv", "id	latitude	longitude	dataset
" *
            "REF_001	46.3	10.7	study_a
REF_001	46.0	11.3	study_a
" *
            "REF_002	46.0	11.3	study_a
REF_003	45.0	9.0	study_b
")
        about = write_file(dir, "a.json", """{"description": "Test reference.",
            "sources": [{"dataset": "study_a", "coordinates": "approximate", "url": "https://example.org/a"}]}""")
        r = load_reference(joinpath(EXAMPLE, "reference.fasta"); metadata = meta, about)
        d = JSON3.read(JSON3.write(ASVAtlas._dataset_json(r; max_queries = 1, max_hits = 1, max_length = 1)))
        @test d.about.description == "Test reference." && d.about.sources[1].coordinates == "approximate"
        @test d.sources.study_a.records == 3 && d.sources.study_a.asvs == 2 && d.sources.study_a.places == 2
        @test d.sources.study_b.records == 1
        @test JSON3.read(JSON3.write(ASVAtlas._dataset_json(load_reference(joinpath(EXAMPLE, "reference.fasta"));
                                                            max_queries = 1, max_hits = 1, max_length = 1))).about === nothing
        for (name, text) in (("bad.json", "{not json"), ("list.json", "[1, 2]"),
                             ("kind.json", """{"sources": [{"coordinates": "guessed"}]}"""),
                             ("url.json", """{"sources": [{"url": "javascript:alert(1)"}]}"""))
            @test_throws ArgumentError load_reference(joinpath(EXAMPLE, "reference.fasta");
                                                      about = write_file(dir, name, text))
        end
    end
    @test ASVAtlas._parse_cli(["--fasta", "r.fa"], noenv).example === nothing
    @test_throws ArgumentError serve(load_reference(joinpath(EXAMPLE, "reference.fasta"));
                                     example = "A"^(2^16 + 1))
    @test ASVAtlas._parse_cli(["--fasta", "r.fa", "--max-network-sequences", "50"],
                              noenv).max_network_sequences == 50
    o = ASVAtlas._parse_cli(["--fasta", "r.fa", "--fields", "site, date,",
                             "--coordinate-digits", "2"], noenv)
    @test o.fields == ["site", "date"] && o.coordinate_digits == 2
    @test ASVAtlas._parse_cli(["--fasta", "r.fa", "--fields", "None"], noenv).fields == String[]
    @test ASVAtlas._parse_cli(String[], Dict("ASVATLAS_FASTA" => "r.fa",
                                             "ASVATLAS_FIELDS" => "site")).fields == ["site"]
    o = ASVAtlas._parse_cli(["--fasta", "r.fa", "--max-length", "5000",
                             "--allowed-hosts", "a.example.org, b.example.org,"], noenv)
    @test o.max_length == 5000 && o.allowed_hosts == ["a.example.org", "b.example.org"]
    env = Dict("ASVATLAS_FASTA" => "/data/r.fa", "ASVATLAS_PORT" => "8081",
               "ASVATLAS_NAME" => "Lakes", "ASVATLAS_MAX_HITS" => " 20 ")
    o = ASVAtlas._parse_cli(String[], env)
    @test o.fasta == "/data/r.fa" && o.port == 8081 && o.name == "Lakes" && o.max_hits == 20
    @test o.metadata === nothing && o.host == "127.0.0.1"
    @test ASVAtlas._parse_cli(["--port", "9"], env).port == 9      # arguments win
    @test ASVAtlas._parse_cli(["--help"], noenv) === nothing
    for args in (String[], ["r.fa"], ["--fasta"], ["--fasta", "r.fa", "--port", "x"],
                 ["--fasta", "r.fa", "--colour", "red"],
                 ["--fasta", "r.fa", "--coordinate-digits", "x"])
        @test_throws ArgumentError ASVAtlas._parse_cli(args, noenv)
    end
    @test redirect_stdout(devnull) do
        ASVAtlas.main(["--help"])
    end == 0
    @test redirect_stderr(devnull) do
        ASVAtlas.main(["--port", "x"])
    end == 1
    @test redirect_stderr(devnull) do
        ASVAtlas.main(["--fasta", joinpath(EXAMPLE, "missing.fasta")])
    end == 1
end

end
