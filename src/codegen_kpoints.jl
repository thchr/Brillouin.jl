using StaticArrays: SVector
using LinearAlgebra: dot, cross, norm

# ---------------------------------------------------------------------------------------- #
# EXPR UTILITIES

# extract all the `Symbol`s referenced in an `Expr`
function get_expr_symbols!(syms::Vector{Symbol}, ex::Expr)
    for arg in ex.args
        if arg isa Symbol
            arg ∈ syms || push!(syms, arg)
        elseif arg isa Expr
            get_expr_symbols!(syms, arg)
        end
    end
    return syms
end
get_expr_symbols(ex::Expr) = get_expr_symbols!(Symbol[], ex)

# ---------------------------------------------------------------------------------------- #
# LOAD `pathsd`, `pointsd`, & `paramsd` DICTIONARIES
include(joinpath(@__DIR__, "..", "assets", "data", "data-SeeK.jl"))

# ---------------------------------------------------------------------------------------- #
# CODEGEN TO CREATE FUNCTIONS FOR EXTENDED BRAVAIS TYPES

# Generate method definitions `$(bt)_points` (`bt` denotes the extended Bravais type) for
# each dataset using metaprogramming. If the datasets require arguments, the arguments are
# listed in the order :a, :b, :c, :β and includes only the featured set of arguments.
for (bt, points) in pointsd
    fn = Symbol(string(bt)*"_points")
    params = get(paramsd, bt, nothing)
    paths  = pathsd[bt]

    # --- build up an expression for the dictionary of k-point labels and vectors ---
    ex = :(Dict())
    for (klab, point_expr) in points
        # due to Julia issue #39240, it is more efficient if we explicitly put `Float64`
        # into the Expr, due to its heterogeneous elements
        point_expr′ = Expr(:ref, :Float64, point_expr.args...) # [x,y,...] to Float64[x,y,...]
        push!(ex.args, Expr(:call, :(=>), QuoteNode(klab), point_expr′))
    end

    if params === nothing
        # --- only constant points; define "constant" function of `Rs` ---
        @eval $fn(Rs::AbstractVector{<:SVector{3, <:Real}}) = $ex

    else
        # --- free arguments; define constant dependent on `Rs` ---

        # write out the all the parameter definitions as expressions
        params_ex = quote end
        for (paramlab, param_expr) in params
            push!(params_ex.args, Expr(:(=), :($(paramlab)), param_expr))
        end
        # find out what the arguments are; can only be {:a, :b, :c} and {:sinβ, :cosβ}
        syms = get_expr_symbols(params_ex) # get all referenced symbol in `params_ex`
        # create list of typed arguments, e.g. (:(a::Real), :(b::Real), ...)
        setup_ex = quote end
        for (i, sym) in enumerate((:a, :b, :c))
            if sym ∈ syms
                push!(setup_ex.args, Expr(:(=), sym, :(norm(Rs[$i]))))
            end
            #sym ∈ syms && push!(args, Expr( :(::), sym, :Real))
        end
        # define "cosβ = cos(β)" and "sinβ = sin(β)" if we can find :cosβ or :sinβ in the
        # argument list; in practice, they only ever occur together in SeeK, so we can just
        # check for one and add both.
        if :cosβ ∈ syms
            for (i, sym) in zip((1,3), (:a, :c))
                # if :a or :c are not in argslist, then we need to add them, since we need
                # them to compute :cosβ and :sinβ as dot and cross products; β is ∠(R₃, R₁)
                sym ∉ syms && push!(setup_ex.args, Expr(:(=), sym, :(norm(Rs[$i]))))
            end
            push!(setup_ex.args, Expr(:(=), :cosβ, :(dot(Rs[3], Rs[1])/(c*a))))
            push!(setup_ex.args, Expr(:(=), :sinβ, :(sqrt(1 - cosβ^2)))) # acos(sin(x)) = √(1-x²)
        end

        # --- make function of the featured subset of :a, :b, :c, and :β, in that order ---
        @eval function $fn(Rs::AbstractVector{<:SVector{3, <:Real}})
            $setup_ex
            $params_ex
            return $ex
        end
    end

    #=
    # see what methods we generated and check that they "work":
    fn′ = getfield(Main, fn) # get the actual function, not just its symbol
    mth = only(methods(fn′))
    println("\n", mth)    
    display(fn′([rand(SVector{3,Float64}) for _ in 1:3])) 
    =#
end