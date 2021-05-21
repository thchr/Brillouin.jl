using Brillouin
using Test
using StaticArrays

# defines `test_show` to test `show` methods with good error diffs, using DeepDiffs.jl
include("testutils_show.jl")

# ---------------------------------------------------------------------------------------- #
include("kpaths.jl")
include("wignerseitz.jl")