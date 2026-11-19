gap> START_TEST("Forms: matobj/quadformbymatrix.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[1,0,0,0],[0,3,0,0],[0,0,0,6],[0,0,6,0]]*Z(7)^0, GF(7));
<4x4 plist matrix over GF(7)>
gap> form := QuadraticFormByMatrix(mat,GF(7));
< quadratic form >
gap> Display(form);
Quadratic form
Gram Matrix:
<immutable 4x4 plist matrix over GF(7):
[[ Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7) ]
 [ 0*Z(7), Z(7), 0*Z(7), 0*Z(7) ]
 [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^5 ]
 [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ]
]>
gap> gf := GF(2^2);
GF(2^2)
gap> mat := ToMatObj(InvariantQuadraticForm( SO(-1, 4, 4) )!.matrix, gf);
<4x4 plist matrix over GF(2^2)>
gap> form := QuadraticFormByMatrix( mat, gf );
< quadratic form >
gap> Display(form);
Quadratic form
Gram Matrix:
<immutable 4x4 plist matrix over GF(2^2):
[[ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), Z(2^2)^2, Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^2)^2 ]
]>
gap> STOP_TEST("matobj/quadformbymatrix.tst", 10000 );
