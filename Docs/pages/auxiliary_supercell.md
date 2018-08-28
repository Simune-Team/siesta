title: Auxiliary Supercell

For periodic systems SIESTA uses an auxiliary supercell as an indexing
tool to represent all the orbital interactions present in the
Hamiltonian (and other matrices). These interactions can extend to
orbitals beyond the unit cell, so extra bookeeping is needed
beyond the simple \(H_{\mu\nu}\) notation used previously.

Consider the simplified system depicted in the figure below. The unit
cell contains two atoms (open and black circles) with a single orbital
each (orbitals are not marked to avoid clutter). The orbitals of the
atoms in the unit cell at the center of the figure overlap with those
of other atoms in the crystal whenever the distance between atoms is
less than the sum of the orbital ranges. In the figure these
interactions extend to the first-neighboring cells only. In this
two-dimensional example, the atoms in the unit cell interact with
those in nine cells, and this geometrical construct is fit to
represent the complete set of non-zero Hamiltonian matrix elements
\(H_{\mu\nu}({\bf R}) \), where \(\mu\) and \(\nu\) are orbital
indexes (ranging over the unit cell only) and \({\bf R}\) is a
relative vector from the unit cell. This vector could just be a cell
lattice vector, but in practice SIESTA uses the relative vector
between the atoms involved.

![Orbital Interactions](|media|/Interactions.png "Orbital interactions")
{: style="text-align: center" ; width="50%" }

There are multiple image interactions for the same pair, each
associated with a different \({\bf R}\).  The set of nine cells is the
auxiliary supercell. In this case it is a 3x3 repetition of the unit
cell, and the "3" can be seen as "2x1+1": the central cell plus one
neighboring cell in either side.  We can also represent the
interactions with a rectangular matrix:

![RectangularMatrix](|media|/RectangularMatrix.png "Interaction
Matrix")


which has a row for each unit-cell orbital, and as many columns as
orbitals in the nine cells involved in the interactions. Any matrix
element $$H_{\mu\nu}({\bf R}) = <\phi_\mu({\bf R=0})|H|\phi_\nu({\bf
R})>$$ can be stored by itself in the appropriate slot.

When it it time to build the Hamiltonian matrix for a given k-point \({\bf k}\):

$$ H_{\mu\nu}({\bf k}) =
\sum_{\bf R} { H_{\mu\nu}({\bf R}) e^{i{\bf k}\cdot{\bf R}}} $$

every slot's contribution, with the appropriate phase, is folded back
and reduced into the left-most square matrix, as the arrows in the
figure indicate.

For periodic systems with large unit cells, k-point sampling is not
really necessary, and the phases involved are all 1 (formally only the
\(\Gamma\) point \({\bf k=0}\) is used). In this case the auxiliary supercell
is not strictly necessary: the matrix elements can be reduced on the
fly to the unit-cell square interaction matrix, and other operations
can be similarly folded automatically throughout. It is possible,
however, to request that an auxiliary supercell be used, since the
extra level of bookeeping can be useful for other purposes (e.g. for
COHP analysis of the orbital-pairs contributions to the energy).

In the program, the auxiliary supercell is handled in several key routines:

* The need for a supercell is assessed in [[siesta_init]], by checking whether k-points are
  going to be used anywhere in the program.
* The size of the supercell needed is stored in the [[siesta_geom:nsc(variable)]] array
* The supercell offsets for each supercell index are stored in the [[siesta_geom:isc_off(variable)]]
  array.
* Indexing arrays live in [[sparse_matrices]] and are initialized in [[state_init]]

When an auxiliary supercell is used, the "column" stored in the
`listh` array defined in the [sparsity](./sparse.html) page ranges
over the long dimension of the rectangular interaction matrix. The
"unit-cell" column to which the folding is done is recorded in the
array `indxuo`. Its contents, however, can be really computed on the
fly, since the column indexes in the rectangular matrix are just
juxtapositions of blocks of size `no_u`.

The following idiom can be used to go through the arrays:
```fortran
    do io = 1, no_u                                                                 
      do ind = listhptr(io) + 1, lishptr(io) + numh(io)                             
        is = (listh(ind) - 1) / no_u + 1                                            
        col = mod(listh(ind) - 1, no_u) + 1     ! this is col=indxuo(listh(ind))
        H(io, col) = H(io, col) + H_sparse(ind) * exp(-1j * sum(xijo(:, ind) * k(:)))   
      end do                                                                        
    end do                  
```

