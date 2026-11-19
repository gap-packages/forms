gap> START_TEST("Forms: matobj/bilformbymatrix.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj(IdentityMat(4, GF(9)), GF(9));
<4x4 plist matrix over GF(3^2)>
gap> form := BilinearFormByMatrix(mat,GF(9));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
<immutable 4x4 plist matrix over GF(3^2):
[[ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
]>
gap> mat := ToMatObj([[0*Z(2),Z(16)^12,0*Z(2),Z(4)^2,Z(16)^13],
>    [Z(16)^12,0*Z(2),0*Z(2),Z(16)^11,Z(16)],
>    [0*Z(2),0*Z(2),0*Z(2),Z(4)^2,Z(16)^3],
>    [Z(4)^2,Z(16)^11,Z(4)^2,0*Z(2),Z(16)^3],
>    [Z(16)^13,Z(16),Z(16)^3,Z(16)^3,0*Z(2) ]], GF(16));
<5x5 plist matrix over GF(2^4)>
gap> form := BilinearFormByMatrix(mat,GF(16));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
<immutable 5x5 plist matrix over GF(2^4):
[[ 0*Z(2), Z(2^4)^12, 0*Z(2), Z(2^2)^2, Z(2^4)^13 ]
 [ Z(2^4)^12, 0*Z(2), 0*Z(2), Z(2^4)^11, Z(2^4) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^2)^2, Z(2^4)^3 ]
 [ Z(2^2)^2, Z(2^4)^11, Z(2^2)^2, 0*Z(2), Z(2^4)^3 ]
 [ Z(2^4)^13, Z(2^4), Z(2^4)^3, Z(2^4)^3, 0*Z(2) ]
]>
gap> mat := ToMatObj([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]*Z(7)^0, GF(7));
<4x4 plist matrix over GF(7)>
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> WittIndex(form);
1
gap> form := BilinearFormByMatrix(ChangedBaseDomain(mat, GF(49)),GF(49));
< bilinear form >
gap> WittIndex(form);
2
gap> STOP_TEST("matobj/bilformbymatrix.tst", 10000 );
