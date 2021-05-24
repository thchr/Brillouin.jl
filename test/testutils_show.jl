# test print with nicely printed diff on failures (from 
# https://github.com/invenia/PkgTemplates.jl/blob/master/test/runtests.jl): access point
# is `test_show`.

using DeepDiffs: deepdiff
function print_diff(a, b)
    old = Base.have_color
    @eval Base have_color = true
    try
        println(deepdiff(a, b))
    finally
        @eval Base have_color = $old
    end
end

function test_text(expected::AbstractString, observed::AbstractString)
    if expected == observed
        @test true
    else
        print_diff(expected, observed)
        @test expected == observed
    end
end

function test_show(v, observed::AbstractString; mime=MIME"text/plain"())
    test_text(repr(mime, v), observed)
end