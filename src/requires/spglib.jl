# Functionalities that are available when Spglib is loaded

using Bravais: reciprocalbasis
using StaticArrays
using LinearAlgebra
import ..KPaths: irrfbz_path

"""
    irrfbz_path(cell::Spglib.Cell)  -->  ::KPath{D}

Returns a **k**-path for an arbitrary 3D unit cell `cell` provided as an `Spglib.Cell`.
The **k**-path is equivalent to the standard **k**-path given by `irrfbz_path` for the
associated standard primitive lattice, but adapted to the possibly non-standard setting of
`cell`.

If `cell` is a supercell of a smaller primitive cell, the standard **k**-path of the
associated primitive cell is returned (in the basis of supercell reciprocal lattice).
"""
function irrfbz_path(cell::Spglib.Cell)
    # extract a standardized primitive basis `pRs` assoc. w/ `cell` via Spglib
    dset = Spglib.get_dataset(cell)
    sgnum = dset.spacegroup_number
    pRs = SVector{3}(SVector{3,Float64}(col) for col in eachcol(dset.std_lattice))

    # print warning if the input cell is a supercell (without any distortion)
    if !isapprox(det(cell.lattice), det(dset.primitive_lattice)) # check volumes agree
        @warn "The provided cell is a supercell: the returned k-path is the standard k-path " *
              "of the associated primitive cell in the basis of the supercell reciprocal lattice." cell
    end

    # Calculate kpath for standard primitive lattice
    kp = irrfbz_path(sgnum, pRs)

    # Now, we convert from the standard primitive lattice to the original lattice.
    # The conversion formula is:
    #     `cell.lattice = rotation * dset.primitive_lattice * transformation`
    # - `rotation` rotates the lattice vectors in 3D space.
    # - `transformation` is an integer-valued matrix that represents the linear combination
    #   of the lattice basis vectors.
    # `transformation` is not used, but is commented here just for information. It can be
    # obtained via
    #   `transformation = inv(primitivebasismatrix(centering(sgnum, 3))) * dset.transformation_matrix'`
    rotation = dset.std_rotation_matrix

    # Rotate k-points in Cartesian space by `rotation`
    pRs_original = SVector{3}(SVector{3,Float64}(col) for col in eachcol(cell.lattice))
    pGs_original = reciprocalbasis(pRs_original)
    cartesianize!(kp)
    for (lab, kv) in points(kp)
        points(kp)[lab] = rotation * kv
    end
    return latticize(kp, pGs_original)
end