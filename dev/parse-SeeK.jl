using DelimitedFiles
using Downloads
using Tar
using CodecZlib

# ---------------------------------------------------------------------------------------- #
# MISC CONVERSION UTILITIES

if !(@isdefined GREEKS)
    const GREEKS = Dict("GAMMA" => 'Γ', "LAMBDA" => 'Λ', "SIGMA" => 'Σ', "DELTA"=>'Δ')
end
if !(@isdefined SUBSCRIPTS)
    const SUBSCRIPTS = Dict('0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄',
                            '5'=>'₅', '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉')
end

function greekify(s::String)
    for (key, substitute) in GREEKS
        if occursin(key, s)
            return replace(s, key=>substitute)
        end
    end
    return s
end

function subscriptify(s::String)
    m = match(r"_[0-9]", s)
    m === nothing && return s
    sub    = m.match
    newsub = SUBSCRIPTS[sub[end]]
    return replace(s, sub => newsub)   
end

# ---------------------------------------------------------------------------------------- #
# DOWNLOAD THE SEEK-PATH REPO; FIND DIRECTORY WITH BAND_PATH_DATA
user = "giovannipizzi"
repo = "seekpath"
url  = "https://github.com/$(user)/$(repo)/tarball/master" # tarball url

# create a temporary file path to download in and then manipulate the data
mktemp(mktempdir()) do path, io
    # download a (shallow) tarball version of SeeK's source code
    Downloads.download(url, io)
    # extract compressed tarball to to tmpdir/ using Tar.jl and CodecZlib.jl
    tmpdir = open(path) do targz_io
        tar = GzipDecompressorStream(targz_io)
        Tar.extract(tar)
    end

    # find the extracted content, which is always of the form $user-$repo-$hash
    contents = readdir(tmpdir)
    i = findfirst(str -> contains(str, user), contents)::Int
    localdir = contents[i]
    basedir  = joinpath(tmpdir, localdir) # path to "root" of repo directory
    
    # directory containing all the extended Bravais type data we want
    datadir  = joinpath(basedir, "seekpath", "hpkot", "band_path_data")

    # license file
    license_file = joinpath(basedir, "LICENSE.txt")

    # ------------------------------------------------------------------------------------ #
    # EXTRACT ALL EXTENDED BRAVAIS TYPES
    bravaistypes = Symbol.(filter!(s->endswith(s, r"[1-3]"), readdir(datadir)))

    pathsd  = Dict{Symbol, Vector{Vector{Symbol}}}()
    pointsd = Dict{Symbol, Vector{Pair{Symbol, Expr}}}()
    paramsd = Dict{Symbol, Vector{Pair{Symbol, Expr}}}()

    for bt in bravaistypes
        btstr = string(bt)
        # -------------------------------------------------------------------------------- #
        # READ PATHS/SEGMENTS
        raw_paths = readdlm(joinpath(datadir, btstr, "path.txt"), String)
        # adjust formatting
        raw_paths .= greekify.(raw_paths)
        raw_paths .= subscriptify.(raw_paths)
        raw_paths = Symbol.(raw_paths)

        # split into connected paths
        pathsd[bt] = [[(first(raw_paths))]]
        last_stop = first(raw_paths)
        i = 1
        for (start, stop) in eachrow(raw_paths)
            if start == last_stop
                push!(pathsd[bt][i], stop) # extend path
            else
                push!(pathsd[bt], [start, stop]) # start new path
                i += 1 # continue on new path
            end
            last_stop = stop
        end
        
        # -------------------------------------------------------------------------------- #
        # READ POINTS
        raw_points = readdlm(joinpath(datadir, btstr, "points.txt"), String)
        labels = Symbol.(subscriptify.(greekify.(raw_points[:,1])))
        values_str  = "[".*join.(eachrow(raw_points[:,2:4]), ',').*']' # as `String`s
        values_repr = Meta.parse.(values_str) # as `Expr`s
        pointsd[bt] = [label => value for (label, value) in zip(labels, values_repr)]
        
        # -------------------------------------------------------------------------------- #
        # READ VARIABLES/PARAMETERS, IF ANY
        param_filepath = joinpath(datadir, btstr, "k_vector_parameters.txt")
        (read(param_filepath, String) == "") && continue # for empty files, nothing to add
        raw_params = readdlm(param_filepath, String)
        paramsd[bt] = [Symbol(label) => Meta.parse(replace(value, "beta"=>'β')) 
                                                    for (label, value) in eachrow(raw_params)]
    end

    # ------------------------------------------------------------------------------------ #
    # FILTER OUT POINTS FROM `pointsd` THAT ARE NOT REFERENCED IN `pathsd`
    # NB: For a few extended Bravais types, SeeK includes points that are not referenced
    #     in the paths; I haven't found a reason/use-case for these extra points, so it
    #     seems easier just to get rid of them entirely (then we can e.g. guarantee that
    #     every point point will be referenced)
    for bt in keys(pointsd)
        points = pointsd[bt]
        paths  = pathsd[bt]
        remove_idxs = Int[]
        for (idx, (klab, kvexpr)) in enumerate(points)
            if klab ∉ Iterators.flatten(paths)
                push!(remove_idxs, idx)
            end
        end
        deleteat!(points, remove_idxs)
    end

    # ------------------------------------------------------------------------------------ #
    # WRITE DICTIONARIES TO A .JL FILE; USE JULIA'S PARSER AS FILE FORMAT
    data_dir = joinpath(@__DIR__, "..", "assets")
    isdir(data_dir) || mkpath(data_dir)
    open(joinpath(data_dir, "data-SeeK.jl"), "w") do io
        println(io,
        "# ", "-"^76, "\n",
        """
        # This file was autogenerated by /$(last(splitpath(@__DIR__)))/$(basename(@__FILE__)).
        # It contains data included in the SeeK-path package (github.com/$user/$repo).
        # The SeeK-path package is licensed under the MIT license, a copy of which is
        # included below.
        """,
        "# ", "-"^76, "\n"
        )

        for name in (:pathsd, :pointsd, :paramsd)
            # obtain the variable via its name (can't use `getfield(Main, name)` since this 
            # is a local scope)
            d = Base.@locals()[name]

            print(io, "const ", name, "_3d = Dict(")
            for (key, val) in d
                println(io)
                print(io, "   :", key, " => ")
                if val isa Vector{Pair{Symbol, Expr}}
                    print(io, '[')
                    for (i, (sym, ex)) in enumerate(val)
                        print(io, ':', sym, " => :(", ex, ")")
                        i != length(val) && print(io, ", ")
                    end
                    print(io, "],")
                else
                    print(io, val, ',')
                end
            end
            println(io, "\n)\n")
        end

        # copy over licence information from SeeK-path
        license_lines = readlines(license_file)
        println(io, "# ------------------ COPY OF ", uppercase(user*"/"*repo),
                    " LICENSE ------------------")
        foreach(line -> println(io, "# ", line), license_lines)
        print(io, "# ", "-"^76)
    end
end
# ---------------------------------------------------------------------------------------- #