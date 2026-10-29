# Delimited text files with a header row. Fields may be quoted with `"` (a
# doubled quote stands for one quote), so that quoted fields can hold the
# delimiter and line breaks.

_guess_delim(path::AbstractString) = lowercase(splitext(path)[2]) == ".csv" ? ',' : '\t'

"""
    read_delimited(path; delim = nothing) -> (header, rows)

Read a delimited file with a header row. The delimiter is `,` for `.csv`
files and a tab otherwise, unless given. Header names are stripped of
surrounding spaces; blank lines are skipped. Every row must have as many
fields as the header.
"""
function read_delimited(path::AbstractString; delim::Union{Nothing,Char} = nothing)
    d = delim === nothing ? _guess_delim(path) : delim
    text = read(path, String)
    startswith(text, '﻿') && (text = text[nextind(text, 1):end])
    rows = _parse_delimited(text, d)
    isempty(rows) && throw(ArgumentError("$path is empty"))
    header = [String(strip(h)) for h in popfirst!(rows)]
    allunique(header) || throw(ArgumentError("$path has repeated column names"))
    for (i, r) in enumerate(rows)
        length(r) == length(header) ||
            throw(ArgumentError("$path: data row $i has $(length(r)) fields, the header $(length(header))"))
    end
    return header, rows
end

function _parse_delimited(text::AbstractString, d::Char)
    rows = Vector{Vector{String}}()
    row = String[]
    field = IOBuffer()
    quoted = false
    i = firstindex(text)
    n = lastindex(text)
    while i <= n
        c = text[i]
        if quoted
            if c == '"'
                j = nextind(text, i)
                if j <= n && text[j] == '"'
                    write(field, '"')
                    i = j
                else
                    quoted = false
                end
            else
                write(field, c)
            end
        elseif c == '"' && position(field) == 0
            quoted = true
        elseif c == d
            push!(row, String(take!(field)))
        elseif c == '\n' || c == '\r'
            push!(row, String(take!(field)))
            push!(rows, row)
            row = String[]
            j = nextind(text, i)
            c == '\r' && j <= n && text[j] == '\n' && (i = j)
        else
            write(field, c)
        end
        i = nextind(text, i)
    end
    quoted && throw(ArgumentError("unterminated quoted field"))
    if position(field) > 0 || !isempty(row)
        push!(row, String(take!(field)))
        push!(rows, row)
    end
    filter!(r -> !(length(r) == 1 && isempty(strip(r[1]))), rows)
    return rows
end

# Field for a tab-separated output file: tabs and line breaks become spaces.
_tsv_field(x) = replace(string(x), r"[\t\r\n]" => " ")
_tsv_field(::Missing) = ""
_tsv_field(x::AbstractFloat) = isnan(x) ? "" : string(round(x; digits = 6))
