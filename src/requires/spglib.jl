# Functionalities that are available when Spglib is loaded

using Bravais: reciprocalbasis, DirectBasis
using LinearAlgebra
import ..KPaths: irrfbz_path

export irrfbz_path

"""
    irrfbz_path(cell::Spglib.Cell)  -->  ::KPath{D}

Returns a k-path for arbitrary 3D cell in the Spglib format. The k-path is equivalent to the
"standard" k-path given by `irrfbz_path` for the equivalent standard primitive lattice.
Does not give the correct k-path if the input cell is a supercell of a smaller primitive
cell (a warning is printed).
"""
function irrfbz_path(cell::Spglib.Cell)
    # standardize cell
    dset = Spglib.get_dataset(cell)
    sgnum = dset.spacegroup_number
    std_lattice = DirectBasis(collect(eachcol(dset.std_lattice)))

    # If the input cell is a supercell (without any distortion), then the irrfbz algorithm cannot work.
    # Check the volume of the input cell and primitive cell are equal.
    if round(Int, det(cell.lattice) / det(dset.primitive_lattice)) != 1
        @warn "input cell is a supercell. irrfbz Does not give a correct k path."
    end

    # Calculate kpath for standard primitive lattice
    kp_standard = Brillouin.irrfbz_path(sgnum, std_lattice)

    # Now, we convert from the standard primitive lattice to the original lattice.
    # The conversion formula is `cell.lattice = rotation * dset.primitive_lattice * transformation`.
    # `transformation` is not used, but is commented here just for information.
    # transformation = inv(Bravais.primitivebasismatrix(Bravais.centering(sgnum, 3))) * dset.transformation_matrix'
    rotation = dset.std_rotation_matrix

    # Rotate k points in Cartesian space by `rotation`
    recip_basis = reciprocalbasis(DirectBasis(collect(eachcol(Matrix(cell.lattice)))))
    kp_cart = cartesianize(kp_standard)
    for (lab, kv) in points(kp_cart)
        points(kp_cart)[lab] = rotation * kv
    end
    kp = latticize(kp_cart, recip_basis)
    kp
end