#import "@preview/ctheorems:1.1.2": *
#import "@preview/physica:0.9.3": *
#show: thmrules

#set heading(numbering: "1.")
#set par(
  leading: 1em,
  justify: true,
  first-line-indent: 1em,
)
#set enum(numbering: "(1)")
#set math.equation(numbering: "(1)")

#set document(
  title: "Implementation note of moyo",
  author: "Kohei Shinohara (@lan496)",
)

#outline()

= Space group

== Symmetry operation search

== Space-group type identification

=== Geometric-crystal-class identification

=== Arithmetic-crystal-class identification

=== Space-group type identification

== Cell standardization

After we identify the space-group type of a given space group, we standardize the cell, which includes symmetrization of the positions and basis vectors and adjustment of a transformation matrix.

Suppose we have a transformation to a DB space group in a primitive basis.
$
  (bold(P), bold(p))^(-1) cal(G)_"prim" (bold(P), bold(p)) = cal(G)_"DB,prim"
$
Since $(bold(P), bold(p))^(-1) (bold(I), bold(t)) (bold(P), bold(p)) = (bold(I), bold(P)^(-1) bold(t)) in cal(T)_"DB,prim"$ for all $bold(t) in bb(Z)^3$, $bold(P)$ is a unimodular matrix.

Let $bold(A)_"prim"$ be basis vectors for the primitive input cell.
Here, we would like to standardize the transformed basis vectors $bold(A)_"prim" bold(P)$ not to be skewed.

=== Triclinic

In this case, $cal(G) = cal(T)$ ($P 1$) or $cal(G) = {1, overline(1)} times cal(T)$ ($P overline(1)$).
Since, $(bold(P), bold(p))^(-1) overline(1) (bold(P), bold(p)) = overline(1)$, there is no need to care about $(bold(P), bold(p))$ in fact.
Therefore, we substitute $(bold(P), bold(p))$ with $(bold(P)_"Niggli", bold(0))$, where $bold(A)_"prim" bold(P)_"Niggli"$ is a Niggli-reduced basis of $bold(A)_"prim"$.

=== Monoclinic

Let $bold(Q)_c in bb(Z)^(3 times 3)$ is a conventional centering matrix from a primitive to conventional cell.
By definition, we can write
$
  bold(A)_"prim" bold(P) bold(Q)_c =: bold(A)_"prim" mat( bold(m)_1, bold(m)_2, bold(n) ),
$
where $bold(A)_"prim" bold(n)$ is a unique axis of the monoclinic cell.
Here, we assume unique axis $bold(c)$ for notation simplicity.

We denote a Minkowski basis of $bold(A)_"prim" angle.l bold(m)_1, bold(m)_2 angle.r$ as ${ bold(A)_"prim" bold(m)'_1, bold(A)_"prim" bold(m)'_2 }$,
$
  bold(A)_"prim" mat( bold(m)'_1, bold(m)'_2 ) = bold(A)_"prim" mat( bold(m)_1, bold(m)_2 ) bold(S).
$

$
  bold(P)' := mat(
    bold(S), bold(0);
    bold(0), 1,
  ) in "SL"_3(bb(Z)) \
  bold(A)_"prim" bold(P) bold(Q)_c bold(P)' = bold(A)_"prim" mat( bold(m)'_1, bold(m)'_2, bold(n) )
$

// Next, we prove $bold(P)'$ does not change the space group as a set.

$
  bold(W) =
    mat(
      1, 0, 0;
      0, 1, 0;
      0, 0, 1
    ), // 1
    mat(
      -1, 0, 0;
      0, -1, 0;
      0, 0, 1;
    ), // 2
    mat(
      -1, 0, 0;
      0, -1, 0;
      0, 0, -1
    ), // -1
    mat(
      1, 0, 0;
      0, 1, 0;
      0, 0, -1
    )  // m = -2
$

$
  (bold(P), bold(p)) (bold(P)', bold(0)) = (bold(P)bold(P)', bold(p))
$

=== Other lattice systems

In these cases, $bold(A)_"prim" bold(P)$ is not skewed.
That is, all cell angles are fixed by the symmetry.


= Magnetic space group
