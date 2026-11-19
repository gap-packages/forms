gap> START_TEST("Forms: matobj/hermitianformbymatrix.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> gf := GF(3^2);
GF(3^2)
gap> mat := ToMatObj(IdentityMat(4, gf), gf);
<4x4 plist matrix over GF(3^2)>
gap> form := HermitianFormByMatrix( mat, gf );
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
<immutable 4x4 plist matrix over GF(3^2):
[[ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
]>
gap> mat := ToMatObj([[Z(11)^0,0*Z(11),0*Z(11)],[0*Z(11),0*Z(11),Z(11)],
>     [0*Z(11),Z(11),0*Z(11)]], GF(121));
<3x3 plist matrix over GF(11^2)>
gap> form := HermitianFormByMatrix(mat,GF(121));
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
<immutable 3x3 plist matrix over GF(11^2):
[[ Z(11)^0, 0*Z(11), 0*Z(11) ]
 [ 0*Z(11), 0*Z(11), Z(11) ]
 [ 0*Z(11), Z(11), 0*Z(11) ]
]>
gap> STOP_TEST("matobj/hermitianformbymatrix.tst", 10000 );
