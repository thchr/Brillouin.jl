var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = Brillouin","category":"page"},{"location":"api/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Brillouin.WignerSeitz.Cell\nBrillouin.KPaths.KPath\nBrillouin.KPaths.KPathInterpolant","category":"page"},{"location":"api/#Brillouin.WignerSeitz.Cell","page":"API","title":"Brillouin.WignerSeitz.Cell","text":"struct Cell{D} <: AbstractArray{Array{StaticArrays.SArray{Tuple{D}, Float64, 1, D}, 1}, 1}\n\nverts::Array{StaticArrays.SVector{D, Float64}, 1} where D\nfaces::Vector{Vector{Int64}}\nbasis::StaticArrays.SArray{Tuple{D}, StaticArrays.SVector{D, Float64}, 1, D} where D\nsetting::Base.RefValue{Brillouin.BasisEnum}\n\n\n\n\n\n","category":"type"},{"location":"api/#Brillouin.KPaths.KPath","page":"API","title":"Brillouin.KPaths.KPath","text":"struct KPath{D} <: Brillouin.KPaths.AbstractPath{Pair{Symbol, StaticArrays.SArray{Tuple{D}, Float64, 1, D}}}\n\npoints::Dict{Symbol, StaticArrays.SVector{D, Float64}} where D\npaths::Vector{Vector{Symbol}}\nbasis::Bravais.ReciprocalBasis\nsetting::Base.RefValue{Brillouin.BasisEnum}\n\n\n\n\n\n","category":"type"},{"location":"api/#Brillouin.KPaths.KPathInterpolant","page":"API","title":"Brillouin.KPaths.KPathInterpolant","text":"struct KPathInterpolant{D} <: Brillouin.KPaths.AbstractPath{StaticArrays.SArray{Tuple{D}, Float64, 1, D}}\n\nkpaths::Array{Array{StaticArrays.SVector{D, Float64}, 1}, 1} where D\nlabels::Vector{Dict{Int64, Symbol}}\nbasis::Bravais.ReciprocalBasis\nsetting::Base.RefValue{Brillouin.BasisEnum}\n\n\n\n\n\n","category":"type"},{"location":"api/#Exported-methods","page":"API","title":"Exported methods","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [Brillouin, Brillouin.WignerSeitz, Brillouin.KPaths]\nPrivate = false\nOrder   = [:function]","category":"page"},{"location":"api/#Brillouin.basis-Tuple{Any}","page":"API","title":"Brillouin.basis","text":"basis(x::Union{KPath, KPathInterpolant, Cell})\n\nReturn the (reciprocal or direct) lattice basis associated with x, in Cartesian coordinates.\n\nMethods in Brillouin will by default return points in the lattice basis, i.e., points are referred to basis(x). This corresponds to the setting(x) == LATTICE. Coordinates may instead be referred to a Cartesian basis, corresponding to setting(x) == CARTESIAN by using cartesianize. The result of  basis(x), however, is invariant to this and always refers to the lattice basis in Cartesian coordinates.\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.cartesianize","page":"API","title":"Brillouin.cartesianize","text":"cartesianize\n\nTransform an object with coordinates in a (explicitly or implicitly specified) lattice basis to an object with coordinates in a Cartesian basis.\n\n\n\n\n\n","category":"function"},{"location":"api/#Brillouin.cartesianize!","page":"API","title":"Brillouin.cartesianize!","text":"cartesianize!\n\nIn-place transform an object with coordinates in a (explicitly or implicitly specified) lattice basis to an object with coordinates in a Cartesian basis.\n\n\n\n\n\n","category":"function"},{"location":"api/#Brillouin.latticize","page":"API","title":"Brillouin.latticize","text":"latticize\n\nTransform an object with coordinates in a Cartesian basis to an object with coordinates in a (explicitly or implicitly specified) lattice basis.\n\n\n\n\n\n","category":"function"},{"location":"api/#Brillouin.latticize!","page":"API","title":"Brillouin.latticize!","text":"latticize!\n\nIn-place transform object with coordinates in a Cartesian basis to an object with coordinates in a (explicitly or implicitly specified) lattice basis.\n\n\n\n\n\n","category":"function"},{"location":"api/#Brillouin.setting-Tuple{Any}","page":"API","title":"Brillouin.setting","text":"setting(x::Union{KPath, KPathInterpolant, Cell})\n\nReturn the basis setting of coordinates in x. The returned value is a member of the BasisEnum enum with member values LATTICE (i.e. coordinates in the basis of the lattice vectors) or CARTESIAN (i.e. coordinates in the Cartesian basis). By default, methods in Brillouin will return coordinates in the LATTICE setting.\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.WignerSeitz.reduce_to_wignerseitz-Union{Tuple{D}, Tuple{StaticArrays.StaticVector{D, <:Real}, AbstractVector{<:StaticArrays.SVector{D, <:Real}}}} where D","page":"API","title":"Brillouin.WignerSeitz.reduce_to_wignerseitz","text":"reduce_to_wignerseitz(v::StaticVector, Vs::BasisLike)  -->  v′\n\nReturn the periodic image v′ of the point v in the basis Vs.\n\nv is assumed to be provided in the lattice basis (i.e., relative to Vs) and v′ is returned in similar fashion.\n\nThe returned point v′ lies in the Wigner-Seitz cell (or its boundary) defined by Vs, has the least possible norm among all equivalent images of v, and differs from v at most by integer lattice translations such that mod(v, 1) ≈ mod(v′, 1).\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.WignerSeitz.wignerseitz-Union{Tuple{AbstractVector{<:StaticArrays.SVector{D, <:Real}}}, Tuple{D}} where D","page":"API","title":"Brillouin.WignerSeitz.wignerseitz","text":"wignerseitz(basis::AbstractVector{<:SVector{D}}; merge::Bool = true, Nmax = 3)\nwignerseitz(basis::AbstractVector{<:AbstractVector}; merge::Bool = true, Nmax = 3)\n                                                            --> Cell{D}\n\nGiven a provided basis, return a Cell{D} containing the vertices and associated (outward oriented) faces of the Wigner-Seitz cell defined by basis in D dimensions. The returned vertices are given in the the basis of basis (see cartesianize! for conversion).\n\nKeyword arguments\n\nmerge (default, true): if true, co-planar faces are merged to form polygonal planar faces (e.g., triangles, quadrilaterals, and ngons generally). If false, raw \"unprocessed\" triangles (D=3) and segments (D=2) are returned instead. merge has no impact for D=1.\nNmax (default, 3): includes -Nmax:Nmax points in the initial lattice used to generate the underlying Voronoi tesselation. It is unwise to set this to anything lower than 3 without explicitly testing convergence; and probably unnecessary to increase it beyond 3 as well.\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.KPaths.interpolate-Tuple{AbstractVector{<:AbstractVector{<:Real}}, Integer}","page":"API","title":"Brillouin.KPaths.interpolate","text":"interpolate(kvs::AbstractVector{<:AbstractVector{<:Real}}, N::Integer)\n                                                --> Vector{Vector{<:Real}}\n\nReturn an interpolated k-path between discrete k-points in kvs, with approximately N interpolation points in total (typically fewer).\n\nNote that, in general, it is not possible to do this so that all interpolated k-points are equidistant; samples are however exactly equidistant across each linear segment defined by points in kvs and approximately equidistant across all segments.\n\nSee also interpolate(::KPath, ::Integer) and splice.\n\nwarning: Future deprecation\nThis method signature is likely to be deprecated in future versions of Brillouin.jl.\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.KPaths.interpolate-Union{Tuple{D}, Tuple{KPath{D}, Integer}} where D","page":"API","title":"Brillouin.KPaths.interpolate","text":"interpolate(kp::KPath, N::Integer)\ninterpolate(kp::KPath; N::Integer, density::Real) --> KPathInterpolant\n\nReturn an interpolant of kp with N points distributed approximately equidistantly across the full k-path (equidistance is measured in a Cartesian metric).\n\nNote that the interpolant may contain slightly fewer or more points than N (typically fewer) in order to improve equidistance. N can also be provided as a keyword argument.\n\nAs an alternative to specifying the desired total number of interpolate points via N, a desired density per unit (reciprocal) length can be specified via the keyword argument density.\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.KPaths.irrfbz_path-Union{Tuple{D}, Tuple{Integer, Any}, Tuple{Integer, Any, Val{D}}} where D","page":"API","title":"Brillouin.KPaths.irrfbz_path","text":"irrfbz_path(sgnum::Integer, Rs, [::Union{Val(D), Integer},]=Val(3))  -->  ::KPath{D}\n\nReturns a k-path (::KPath) in the (primitive) irreducible Brillouin zone for a space group with number sgnum, (conventional) direct lattice vectors Rs, and dimension D. The path includes all distinct high-symmetry lines and points as well as relevant parts of the Brillouin zone boundary.\n\nThe dimension D (1, 2, or 3) is specified as the third input argument, preferably as a static Val{D} type parameter (or, type-unstably, as an <:Integer). Defaults to Val(3).\n\nRs refers to the direct basis of the conventional unit cell, i.e., not the primitive  direct basis vectors. The setting of Rs must agree with the conventional setting choices in the International Tables of Crystallography, Volume A (the \"ITA conventional setting\"). If Rs is a subtype of a StaticVector or NTuple, the dimension can be inferred from its (static) size; in this case, this dimension will take precedence (i.e. override, if different) over any dimension specified in the third input argument.\n\nNotes\n\nThe returned k-points are given in the basis of the primitive reciprocal basis in the CDML setting. To obtain the associated transformation matrices between the ITA conventional setting and the CDML primitive setting, see primitivebasismatrix of Bravais.jl](https://thchr.github.io/Crystalline.jl/stable/bravais/) (or, equivalently, the relations defined Table 2 of [1]). To transform to a Cartesian basis, see cartesianize!.\nTo interpolate a KPath, see interpolate(::KPath, ::Integer) and splice(::KPath, ::Integer).\nAll paths currently assume time-reversal symmetry (or, equivalently, inversion symmetry). If neither are present, include the \"inverted\" -k paths manually.\n\nData and referencing\n\n3D paths are sourced from the SeeK-path publication: please cite the original work [2].\n\nReferences\n\n[1] Aroyo et al., Acta Cryst. A70, 126 (2014). [2] Hinuma, Pizzi, Kumagai, Oba, & Tanaka, Band structure diagram paths based on     crystallography, Comp. Mat. Sci. 128, 140 (2017).\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.KPaths.paths-Tuple{KPath}","page":"API","title":"Brillouin.KPaths.paths","text":"paths(kp::KPath) -> Vector{Vector{Symbol}}\n\nReturn a vector of vectors, with each vector describing a connected path between between k-points referenced in kp (see also points(::KPath)).\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.KPaths.points-Tuple{KPath}","page":"API","title":"Brillouin.KPaths.points","text":"points(kp::KPath{D}) -> Dict{Symbol, SVector{D,Float64}}\n\nReturn a dictionary of the k-points (values) and associated k-labels (keys) referenced in kp.\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.KPaths.splice-Tuple{AbstractVector{<:AbstractVector{<:Real}}, Integer}","page":"API","title":"Brillouin.KPaths.splice","text":"splice(kvs::AbstractVector{<:AbstractVector{<:Real}}, N::Integer)\n                                                --> Vector{Vector{<:Real}}\n\nReturn an interpolated k-path between the discrete k-points in kvs, with N interpolation points inserted in each segment defined by pairs of adjacent k-points.\n\nSee also splice(::KPath, ::Integer) and interpolate.\n\nwarning: Future deprecation\nThis method signature is likely to be deprecated in future versions of Brillouin.jl.\n\n\n\n\n\n","category":"method"},{"location":"api/#Brillouin.KPaths.splice-Union{Tuple{D}, Tuple{KPath{D}, Integer}} where D","page":"API","title":"Brillouin.KPaths.splice","text":"splice(kp::KPath, N::Integer) --> KPathInterpolant\n\nReturn an interpolant of kp with N points inserted into each k-path segment of kp.\n\n\n\n\n\n","category":"method"},{"location":"kpaths/#**k**-space-paths","page":"k-space paths","title":"k-space paths","text":"","category":"section"},{"location":"kpaths/#Default-paths","page":"k-space paths","title":"Default paths","text":"","category":"section"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"To generate a k-path for, say, the space group of diamond (space group 227; a cubic face-centered Bravais lattice), we can call irrfbz_path, which will return a minimal path in the irreducible Brillouin zone:","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"using Brillouin\nsgnum = 227\nRs = [[1,0,0], [0,1,0], [0,0,1]] # conventional direct basis\nkp = irrfbz_path(sgnum, Rs)","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"The path data is sourced from the HPKOT paper (or, equivalently, the SeeK-path Python package). ","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"The coordinates of the path are given with respect to the primitive reciprocal basis (here, [[-2π,2π,2π], [2π,-2π,2π], [2π,2π,-2π]]). To convert to a Cartesian basis, we can use cartesianize or cartesianize! (in-place):","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"cartesianize(kp)","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"We can visualize the k-path using PlotlyJS.jl (conversion to a Cartesian basis for plotting is automatic):","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"using PlotlyJS\nPᵏ = plot(kp)\nMain.HTMLPlot(Pᵏ) # hide","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"Usually, it'll be more helpful to understand the path's geometry in the context of the associated Brillouin zone. To visualize this, we can plot the combination of a Cell (created via wignerseitz) and a KPath:","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"pGs = basis(kp)      # primitive reciprocal basis associated with k-path\nc = wignerseitz(pGs) # associated Brillouin zone\nPᶜ⁺ᵏ = plot(c, kp)\nMain.HTMLPlot(Pᶜ⁺ᵏ) # hide","category":"page"},{"location":"kpaths/#Interpolation","page":"k-space paths","title":"Interpolation","text":"","category":"section"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"Interpolation of a KPath structure can be achieved using either interpolate(::KPath, ::Integer) or splice(::KPath, ::Integer), returning a KPathInterpolant. As an example, interpolate(kp, N) returns an interpolation with a target of N interpolation points, distributed as equidistantly as possible (with the distance metric evaluated in Cartesian space):","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"kpi = interpolate(kp, 100)","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"The returned KPathInterpolant implements the AbstractVector interface, with iterants returning SVector{D, Float64} elements. To get a conventional \"flat\" vector, we can simply call collect(kpi).","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"Internally, KPathInterpolant includes additional structure and information: namely, the high-symmetry points and associated labels along the path as well as a partitioning into connected vs. disconnected path segments.","category":"page"},{"location":"kpaths/#Band-structure","page":"k-space paths","title":"Band structure","text":"","category":"section"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"The additional structure of KPathInterpolation enables convenient and clear visualizations of band structure diagrams in combination with PlotlyJS.","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"To illustrate this, suppose we are considering a tight-binding problem for an s-orbital situated at the 1a Wyckoff position. Such a problem has a single band with dispersion [1] (assuming a cubic side length a = 1):","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"epsilon(mathbfk)\n=\n4gammaBigl(\n    cos tfrac12k_x cos tfrac12k_y +\n    cos tfrac12k_y cos tfrac12k_z +\n    cos tfrac12k_z cos tfrac12k_x\n    Bigr)","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"with k_xyz denoting coordinates in a Cartesian basis (which are related to the coordinates k_123 in a primitive reciprocal basis by k_x = 2π(-k_1+k_2+k_3), k_x = 2π(k_1-k_2+k_3), and k_z = 2π(k_1+k_2-k_3)).","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"We can calculate the energy band along our k-path using the interpolation object kpi. To do so, we define a function that implements epsilon(mathbfk) and broadcast it over the elements of kpi:","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"function ϵ(k; γ::Real=1.0)\n    kx = 2π*(-k[1]+k[2]+k[3])\n    ky = 2π*(+k[1]-k[2]+k[3])\n    kz = 2π*(+k[1]+k[2]-k[3])\n    return 4γ * (cos(kx/2)*cos(ky/2) + cos(ky/2)*cos(kz/2) + cos(kz/2)*cos(kx/2))\nend\nband = ϵ.(kpi)","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"Finally, we can visualize the associated band using a Brillouin-overloaded PlotlyJS plot-call:","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"P = plot(kpi, [band])\nMain.HTMLPlot(P, 525) # hide","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"If we have multiple bands, say epsilon_1(mathbfk) = epsilon(mathbfk) and epsilon_2(mathbfk) = 20 - tfrac12epsilon(mathbfk), we can easily plot that by collecting the bands in a single vector (or concatenating into a matrix):","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"band1 = ϵ.(kpi)\nband2 = 20 .- (1/2).*band1\nP¹² = plot(kpi, [band1, band2])\nMain.HTMLPlot(P¹², 525) # hide","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"plot(::KPathInterpolant, ::Any, ::Layout)","category":"page"},{"location":"kpaths/#PlotlyJS.plot-Tuple{KPathInterpolant, Any, Layout}","page":"k-space paths","title":"PlotlyJS.plot","text":"plot(kpi::KPathInterpolant, bands, [layout]; kwargs...)\n\nPlot a dispersion diagram for provided bands and k-path interpolant kpi.\n\nbands must be an iterable of iterables of <:Reals (e.g., a Vector{Vector{Float64}}), with the first iteration running over distinct energy bands, and the second running over distinct k-points in kpi. Note that the length of each iterant of bands must equal length(kpi).\n\nAlternatively, bands can be an AbstractMatrix{<:Real}, with columns interpreted as distinct energy bands and rows as distinct k-points.\n\nA layout can be supplied to overwrite default layout choices (set by Brillouin.DEFAULT_PLOTLY_LAYOUT_DISPERSION). Alternatively, some simple settings can be set directly via keyword arguments (see below).\n\nKeyword arguments kwargs\n\nylims: y-axis limits (default: quasi-tight around bands's range)\nylabel: y-axis label (default: \"Energy\")\ntitle: plot title (default: nothing); can be a String or an attr dictionary of PlotlyJS properties\nband_highlights: dictionary of non-default styling for specified band ranges (default: nothing, indicating all-default styling).\nExample: To color bands 2 and 3 black, set band_highlights = Dict(2:3 => attr(color=:black, width=3)). Unlisted band ranges are shown in default band styling.\nannotations: dictionary of hover-text annotations for labeled high-symmetry points in kpi (default: nothing, indicating no annotations). Suitable for labeling of irreps.\nExample: Assume bands 1 and 2 touch at X, but not at Γ. To label this, we set: annotations = Dict(:X => [1:2 => \"touching!\"], :Γ => [1 => \"isolated\", 2 => \"isolated\"]). If a band-range is provided, a single annotation is placed at the mean of the energies at these band-ranges. Alternatively, if the first element of each pair is a non-Integer Real number, it is interpreted as referring to the frequency of the annotation.\n\n\n\n\n\n","category":"method"},{"location":"kpaths/#Treating-unit-cells-in-non-standard-settings","page":"k-space paths","title":"Treating unit cells in non-standard settings","text":"","category":"section"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"irrfbz_path(sgnum, Rs) requires Rs to be provided in a standard setting. Often, the setting of Rs may not be standard and it can be a hassle to convert existing calculations to such a setting. To avoid this, we can provide the unit cell information to irrfbz via Spglib's Cell format (this functionality depends on separately loading the Spglib.jl package). If the provided unit cell is a supercell of a smaller primitive cell, irrfbz_path returns the standard k-path of the smaller primitive cell in the basis of the supercell reciprocal lattice vectors.","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"For example, to construct a k-path for a non-standard trigonal lattice in space group 166 (R-3m):","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"using Spglib\na = 1.0\nc = 8.0\npRs_standard    = [[a*√3/2,  a/2, c/3],\n                   [-a*√3/2, a/2, c/3],\n                   [0, -a, c/3]]\npRs_nonstandard = [[a*√3/2, -a/2, c/3],\n                   [0, a, c/3],\n                   [-a*√3/2, -a/2, c/3]]\n\ncell_standard = Spglib.Cell(pRs_standard, [[0, 0, 0]], [0])\ncell_nonstandard = Spglib.Cell(pRs_nonstandard, [[0, 0, 0]], [0])\n\nkp_standard = irrfbz_path(cell_standard)\nkp_nonstandard = irrfbz_path(cell_nonstandard)","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"Note that the space group symmetry is inferred by Spglib from the atomic positions and provided basis.","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"We can check that the generated k-paths for the standard and non-standard lattices are equivalent by plotting and comparing their k-paths and associated Wigner-Seitz cells:","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"using Bravais # for `reciprocalbasis`\nPˢ  = plot(wignerseitz(reciprocalbasis(pRs_standard)), kp_standard, Layout(title=\"standard cell\"))\nMain.HTMLPlot(Pˢ) # hide","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"Pⁿˢ = plot(wignerseitz(reciprocalbasis(pRs_nonstandard)), kp_nonstandard, Layout(title=\"non-standard cell\"))\nMain.HTMLPlot(Pⁿˢ) # hide","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"","category":"page"},{"location":"kpaths/","page":"k-space paths","title":"k-space paths","text":"[1] See e.g. http://www.physics.rutgers.edu/~eandrei/chengdu/reading/tight-binding.pdf","category":"page"},{"location":"wignerseitz/#Wigner–Seitz-cells","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"","category":"section"},{"location":"wignerseitz/#Unit-cells","page":"Wigner–Seitz cells","title":"Unit cells","text":"","category":"section"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"The Wigner–Seitz cell associated with an arbitrary lattice basis can be generated via wignerseitz. For example, to generate the unit cell of a (primitive) lattice with Bravais type cF, we might write:","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"using Brillouin\n\nRs = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]\ncᴿ = wignerseitz(Rs)","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"Note that the coordinates of the Wigner-Seitz cell vertices are referred to the basis Rs; to convert to Cartesian space, see cartesianize and cartesianize! (in-place).","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"We can plot the generated cells using e.g. PlotlyJS.jl via plot(cᴿ) (or, alternatively, via a 3D-capable backend of Makie.jl):","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"using PlotlyJS\nPᴿ = plot(cᴿ)\nMain.HTMLPlot(Pᴿ) # hide","category":"page"},{"location":"wignerseitz/#Brillouin-zones","page":"Wigner–Seitz cells","title":"Brillouin zones","text":"","category":"section"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"To generate Brillouin zones, we simply give the corresponding reciprocal lattice Gs:","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"Gs = 2π.*[[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]] # reciprocal basis of `Rs`\ncᴳ = wignerseitz(Gs)\nPᴳ = plot(cᴳ)\nMain.HTMLPlot(Pᴳ) # hide","category":"page"},{"location":"wignerseitz/#Two-dimensions","page":"Wigner–Seitz cells","title":"Two dimensions","text":"","category":"section"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"wignerseitz and plot(::Cell) works in two dimensions as well. As an example, we can illustrate the Wigner–Seitz unit cell of graphene (which has a hexagonal \"hp\" Bravais type):","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"using Brillouin, PlotlyJS\n\nRs = [[1.0, 0.0], [-0.5, √3/2]]\ncᴿ = wignerseitz(Rs)\nPᴿ = plot(cᴿ)\nMain.HTMLPlot(Pᴿ) # hide","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"and its associated Brillouin zone:","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"Gs = 2π.*[[1.0, 1/√3], [0.0, 2/√3]]\ncᴳ = wignerseitz(Gs)\nPᴳ = plot(cᴳ)\nMain.HTMLPlot(Pᴳ) # hide","category":"page"},{"location":"wignerseitz/#Reducing-points-to-the-interior-of-Wigner-Seitz-cell","page":"Wigner–Seitz cells","title":"Reducing points to the interior of Wigner-Seitz cell","text":"","category":"section"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"Points can be reduced to their equivalent image in the Wigner-Seitz cell via  reduce_to_wignerseitz. The resulting point has the least possible norm of all equivalent images.","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"As an example, consider the reduction of a point v initially outside the 2D Brillouin zone defined earlier:","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"v = [0.8, 0.2]\nv′ = reduce_to_wignerseitz(v, Gs)","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"We can visualize the reduction by adding the points to the previous plot of the Brillouin zone:","category":"page"},{"location":"wignerseitz/","page":"Wigner–Seitz cells","title":"Wigner–Seitz cells","text":"vᶜ = cartesianize(v, Gs)\nv′ᶜ = cartesianize(v′, Gs)\n\naddtraces!(Pᴳ, scatter(x=vᶜ[1:1],  y=vᶜ[2:2],  hovertext=\"v (unreduced)\"))\naddtraces!(Pᴳ, scatter(x=v′ᶜ[1:1], y=v′ᶜ[2:2], hovertext=\"v′ (reduced)\"))\nMain.HTMLPlot(Pᴳ) # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Brillouin","category":"page"},{"location":"#Brillouin.jl","page":"Home","title":"Brillouin.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Brillouin.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"wignerseitz.md\",\n         \"kpaths.md\",\n         \"api.md\"]","category":"page"}]
}