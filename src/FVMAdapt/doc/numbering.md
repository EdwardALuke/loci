# FVMAdapt Numbering Conventions {#fvmadapt_numbering}

FVMAdapt uses local numbering for each refinement entity. These numbers are not
global mesh entity IDs. They label positions within one edge, face, or cell
refinement tree as that tree is built and replayed.

The coordinate names below are local reference directions. They are not physical
`x`, `y`, and `z` coordinates from the mesh, and they do not imply orthogonal
axes or equal edge lengths. They identify ordered node, edge, face, and child
positions used to interpret split codes, child ordering, and face/edge
mappings.

The diagrams below show canonical reference shapes, not a claim that the
physical mesh cells are squares or cubes. A skewed quadrilateral still has the
same local `x` and `y` directions; those directions follow the deformed physical
edges.

When the entity type matters, this page uses subscripts for the entity family
and superscripts for the local number. For example,
n<sub>q</sub><sup>(0)</sup> is QuadFace-local node zero, and
e<sub>H</sub><sup>(6)</sup> is HexCell-local edge six.

## QuadFace

A `QuadFace` has four corner-node positions and four boundary-edge positions.
Both are numbered from zero. These positions define the face-local reference view
used by `QuadFace` split codes.

<img src="figures/numbering/quadface_numbering.svg"
     alt="QuadFace node and edge numbering"
     width="420">

The face-local origin is n<sub>q</sub><sup>(0)</sup>. Face-local `x` points
from n<sub>q</sub><sup>(0)</sup> to n<sub>q</sub><sup>(1)</sup>; face-local
`y` points from n<sub>q</sub><sup>(0)</sup> to
n<sub>q</sub><sup>(3)</sup>.

This describes local ordering only. The physical edges from
n<sub>q</sub><sup>(0)</sup> to n<sub>q</sub><sup>(1)</sup> and from
n<sub>q</sub><sup>(0)</sup> to n<sub>q</sub><sup>(3)</sup> do not have to be
perpendicular.

The boundary-edge numbers are:

- e<sub>q</sub><sup>(0)</sup>: n<sub>q</sub><sup>(0)</sup> to n<sub>q</sub><sup>(1)</sup>
- e<sub>q</sub><sup>(1)</sup>: n<sub>q</sub><sup>(1)</sup> to n<sub>q</sub><sup>(2)</sup>
- e<sub>q</sub><sup>(2)</sup>: n<sub>q</sub><sup>(3)</sup> to n<sub>q</sub><sup>(2)</sup>
- e<sub>q</sub><sup>(3)</sup>: n<sub>q</sub><sup>(0)</sup> to n<sub>q</sub><sup>(3)</sup>


## HexCell

A `HexCell` is made of six quad faces. Its local view is a reference
hexahedron. Node n<sub>H</sub><sup>(0)</sup> is the cell-local origin with
`(xi, eta, zeta) = (0, 0, 0)`. A physical hexahedron can be skewed or
stretched; the local coordinates below name positions in the reference shape.

The face numbers are named by the constant local coordinate on that face:

<img src="figures/numbering/hexcell_numbering.svg"
     alt="HexCell local node and face numbering"
     width="620">

In the above diagram, numbers bounded by squares are HexCell face numbers
f<sub>H</sub><sup>(i)</sup>, and the unboxed numbers are HexCell node numbers
n<sub>H</sub><sup>(i)</sup>.

| HexCell face | Local coordinate |
| --- | --- |
| f<sub>H</sub><sup>(0)</sup> | `xi = 1` |
| f<sub>H</sub><sup>(1)</sup> | `xi = 0` |
| f<sub>H</sub><sup>(2)</sup> | `eta = 1` |
| f<sub>H</sub><sup>(3)</sup> | `eta = 0` |
| f<sub>H</sub><sup>(4)</sup> | `zeta = 1` |
| f<sub>H</sub><sup>(5)</sup> | `zeta = 0` |

The local node numbering is:

| HexCell node | `xi` | `eta` | `zeta` |
| --- | --- | --- | --- |
| n<sub>H</sub><sup>(0)</sup> | 0 | 0 | 0 |
| n<sub>H</sub><sup>(1)</sup> | 0 | 0 | 1 |
| n<sub>H</sub><sup>(2)</sup> | 0 | 1 | 0 |
| n<sub>H</sub><sup>(3)</sup> | 0 | 1 | 1 |
| n<sub>H</sub><sup>(4)</sup> | 1 | 0 | 0 |
| n<sub>H</sub><sup>(5)</sup> | 1 | 0 | 1 |
| n<sub>H</sub><sup>(6)</sup> | 1 | 1 | 0 |
| n<sub>H</sub><sup>(7)</sup> | 1 | 1 | 1 |

Equivalently, the local node number is `4*xi + 2*eta + zeta` for `xi`, `eta`,
and `zeta` values of zero or one.

The `HexCell` also has 12 local edges. These are HexCell-local edge ids, not
global mesh edge ids. The implementation numbers them by their endpoint nodes:

| HexCell edge | Endpoints |
| --- | --- |
| e<sub>H</sub><sup>(0)</sup> | n<sub>H</sub><sup>(0)</sup> to n<sub>H</sub><sup>(4)</sup> |
| e<sub>H</sub><sup>(1)</sup> | n<sub>H</sub><sup>(1)</sup> to n<sub>H</sub><sup>(5)</sup> |
| e<sub>H</sub><sup>(2)</sup> | n<sub>H</sub><sup>(2)</sup> to n<sub>H</sub><sup>(6)</sup> |
| e<sub>H</sub><sup>(3)</sup> | n<sub>H</sub><sup>(3)</sup> to n<sub>H</sub><sup>(7)</sup> |
| e<sub>H</sub><sup>(4)</sup> | n<sub>H</sub><sup>(0)</sup> to n<sub>H</sub><sup>(2)</sup> |
| e<sub>H</sub><sup>(5)</sup> | n<sub>H</sub><sup>(1)</sup> to n<sub>H</sub><sup>(3)</sup> |
| e<sub>H</sub><sup>(6)</sup> | n<sub>H</sub><sup>(4)</sup> to n<sub>H</sub><sup>(6)</sup> |
| e<sub>H</sub><sup>(7)</sup> | n<sub>H</sub><sup>(5)</sup> to n<sub>H</sub><sup>(7)</sup> |
| e<sub>H</sub><sup>(8)</sup> | n<sub>H</sub><sup>(0)</sup> to n<sub>H</sub><sup>(1)</sup> |
| e<sub>H</sub><sup>(9)</sup> | n<sub>H</sub><sup>(2)</sup> to n<sub>H</sub><sup>(3)</sup> |
| e<sub>H</sub><sup>(10)</sup> | n<sub>H</sub><sup>(4)</sup> to n<sub>H</sub><sup>(5)</sup> |
| e<sub>H</sub><sup>(11)</sup> | n<sub>H</sub><sup>(6)</sup> to n<sub>H</sub><sup>(7)</sup> |

## Face-Local Orientation

A mesh face also has its own `face2node` order. That order can start at a
different corner of the same physical face, and it can run around the face in
the opposite direction from the `HexCell` view. The cyclic order, and therefore
the normal associated with that order, is what fixes the face-local
orientation. Node `0` alone is not enough.

The examples below use f<sub>H</sub><sup>(0)</sup>, the `xi = 1` face. The left
side of each figure shows the face in the HexCell, with HexCell node numbers
and HexCell edge numbers on the highlighted face. The projected copy shows one
possible QuadFace-local view of that same physical face. On the projected copy,
the blue node and edge labels are QuadFace-local numbers, and the blue arrows
show QuadFace-local edge directions.

<img src="figures/numbering/quadface_orientation_q0_h4.svg"
     alt="HexCell xi equals 1 face projected with QuadFace n_q zero at n_H four"
     width="820">

<em>QuadFace-local node zero lands on HexCell node four, and the local edge
directions follow that orientation.</em>

<img src="figures/numbering/quadface_orientation_q0_h6.svg"
     alt="HexCell xi equals 1 face projected with QuadFace n_q zero at n_H six"
     width="820">

<em>The same face-local orientation has been rotated so that QuadFace-local node
zero lands on HexCell node six.</em>

<img src="figures/numbering/quadface_orientation_reversed.svg"
     alt="HexCell xi equals 1 face projected with reversed QuadFace local order"
     width="820">

<em>Here QuadFace-local node zero is back on HexCell node four, but the local
edge directions are reversed, so the face-local normal points in the opposite
direction.</em>

When a hex cell is built, FVMAdapt compares the mesh face's `face2node` order
with the HexCell face order and records the result as `hexOrientCode`. That
orientation code is used when face plans are applied to cells, so split codes,
child ids, and edge ids refer to the same physical subfaces even when the
face-local and cell-local views differ.
