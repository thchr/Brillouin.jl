@precompile_setup begin
    using Bravais: directbasis

    Rs2D_v = [(@SVector [1.0,0]), (@SVector [0,1.0])]
    Rs2D = directbasis(11, Val(2))
    Rs3D_v = [(@SVector [1.0,0,0]), (@SVector [0,1.0,0]), (@SVector [0,0,1.0])]
    Rs3D = directbasis(227, Val(3))

    @precompile_all_calls begin
        function _precompile(sgnum, Rs, Dᵛ::Val{D}) where D
            kp = irrfbz_path(sgnum, Rs, Dᵛ)
            cartesianize(kp)
            interpolate(kp, 10)
            interpolate(kp; density=.1)

            # pGs = basis(kp)      # FIXME: We cannot precompile `wignerseitz` due to some
            # wignerseitz(pGs)     #        pointer shenanigans in DirectQhull.jl.
            # wignerseitz(Rs)      #        If uncommented, throws "ERROR: pointerref:
                                   #        invalid pointer"
            nothing
        end

        _precompile(11,  Rs2D_v, Val(2))
        _precompile(11,  Rs2D,   Val(2))
        _precompile(227, Rs3D_v, Val(3))
        _precompile(227, Rs3D,   Val(3))
    end
end