gap> START_TEST("Forms: matobj/test_forms12.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> f := GF(5);
GF(5)
gap> mat := ToMatObj(IdentityMat(3,f), f);
<3x3 plist matrix over GF(5)>
gap> form := BilinearFormByMatrix(mat,f);
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
<immutable 3x3 plist matrix over GF(5):
[[ Z(5)^0, 0*Z(5), 0*Z(5) ]
 [ 0*Z(5), Z(5)^0, 0*Z(5) ]
 [ 0*Z(5), 0*Z(5), Z(5)^0 ]
]>
gap> mat[1,2] := One(f);
Z(5)^0
gap> Display(form);
Bilinear form
Gram Matrix:
<immutable 3x3 plist matrix over GF(5):
[[ Z(5)^0, 0*Z(5), 0*Z(5) ]
 [ 0*Z(5), Z(5)^0, 0*Z(5) ]
 [ 0*Z(5), 0*Z(5), Z(5)^0 ]
]>
gap> f := GF(5);
GF(5)
gap> mat := ToMatObj(IdentityMat(3,f), f);
<3x3 plist matrix over GF(5)>
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> Display(form);
Quadratic form
Gram Matrix:
<immutable 3x3 plist matrix over GF(5):
[[ Z(5)^0, 0*Z(5), 0*Z(5) ]
 [ 0*Z(5), Z(5)^0, 0*Z(5) ]
 [ 0*Z(5), 0*Z(5), Z(5)^0 ]
]>
gap> mat[1,2] := One(f);
Z(5)^0
gap> Display(form);
Quadratic form
Gram Matrix:
<immutable 3x3 plist matrix over GF(5):
[[ Z(5)^0, 0*Z(5), 0*Z(5) ]
 [ 0*Z(5), Z(5)^0, 0*Z(5) ]
 [ 0*Z(5), 0*Z(5), Z(5)^0 ]
]>
gap> f := GF(81);
GF(3^4)
gap> mat := ToMatObj(IdentityMat(4,f), f);
<4x4 plist matrix over GF(3^4)>
gap> form := HermitianFormByMatrix(mat,f);
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
<immutable 4x4 plist matrix over GF(3^4):
[[ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
]>
gap> mat[1,2] := One(f);
Z(3)^0
gap> Display(form);
Hermitian form
Gram Matrix:
<immutable 4x4 plist matrix over GF(3^4):
[[ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
]>
gap> STOP_TEST("matobj/test_forms12.tst", 10000 );
