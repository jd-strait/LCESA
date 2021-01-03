# LCESA
Landmark-constrained elastic shape analysis tools


**Main programs**

`mygeod2.m` - landmark-constrained registration of 2 input planar curves, where only the curves (represented as two-dimensional vectors) are provided

`mygeod2pieces.m` - landmark-constrained registration of 2 input planar curves, where the curves are already provided as corresponding pieces (i.e., if landmarks were used previously to "split" up curves)

**Main sub programs**

`alignDP.m`, `alignpiecesDP.m` - uses dynamic programming to perform landmark-constrained registration

`alignGD.m`, `alignpiecesGD.m` - use gradient-descent to perform landmark-constrained registration

`cursor.m` - program to select landmarks interactively on curves

`GDReparLC.m` - gradient-descent to solve for optimal re-parameterization function with landmark constraints (used within `alignGD.m` and `alignpiecesGD.m`)

`KarcherMean.m`, `KarcherMeanPieces.m` - finds landmark-constrained elastic shape mean (provided curve is available in full form or in pieces based on landmarks)

`PCAshape.m` - performs tangent PCA for landmark-constrained shapes (requires Karcher mean to be computed already)


E-mail Justin Strait (justin.strait@uga.edu) with any questions about the code or for additional details.
