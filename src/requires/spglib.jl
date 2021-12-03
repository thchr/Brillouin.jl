# Functionalities that are available when Spglib is loaded

using Bravais: reciprocalbasis, DirectBasis
using LinearAlgebra
import ..KPaths: irrfbz_path

"""
    irrfbz_path(cell::Spglib.Cell)  -->  ::KPath{D}

Returns a **k**-path for an arbitrary 3D cell in the Spglib format. The **k**-path is equivalent to the
"standard" **k**-path given by `irrfbz_path` for the associated standard primitive lattice,
but adapted to the possibly non-standard setting of `cell`.

`cell` cannot represent a supercell of a smaller primitive cell.
"""
function irrfbz_path(cell::Spglib.Cell)
    # standardize cell
    dset = Spglib.get_dataset(cell)
    sgnum = dset.spacegroup_number
    std_lattice = DirectBasis{3}(collect(eachcol(dset.std_lattice)))

    # If the input cell is a supercell (without any distortion), then the irrfbz algorithm cannot work.
    if !isapprox(det(cell.lattice), det(dset.primitive_lattice)) # check volumes agree
        error(DomainError(cell, "`cell` is a supercell and untreatable by `irrfbz_path`. " *
            "If the cell represents a standard conventional lattice, use `irrfbz_path(sgnum::Integer, Rs)`"))
    end

    # Calculate kpath for standard primitive lattice
    kp_standard = irrfbz_path(sgnum, std_lattice)

    # Now, we convert from the standard primitive lattice to the original lattice.
    # The conversion formula is `cell.lattice = rotation * dset.primitive_lattice * transformation`.
    # `rotation` rotates the lattice vectors in 3D space.
    # `transformation` is an integer-valued matrix that represents the linear combination of the lattice basis vectors.
    # `transformation` is not used, but is commented here just for information.
    # transformation = inv(Bravais.primitivebasismatrix(Bravais.centering(sgnum, 3))) * dset.transformation_matrix'
    rotation = dset.std_rotation_matrix

    # Rotate k points in Cartesian space by `rotation`
    recip_basis = reciprocalbasis(collect.(eachcol(cell.lattice)))
    kp_cart = cartesianize(kp_standard)
    for (lab, kv) in points(kp_cart)
        points(kp_cart)[lab] = rotation * kv
    end
    kp = latticize(kp_cart, recip_basis)
    return kp
end