gap> START_TEST("Forms: matobj/bg_th_ex8.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[Z(16)^3,1,0,0],[0,Z(16)^5,0,0],
>              [0,0,Z(16)^3,1],[0,0,0,Z(16)^12]]*Z(16)^0, GF(16));
<4x4 plist matrix over GF(2^4)>
gap> qform := QuadraticFormByMatrix(mat,GF(16));
< quadratic form >
gap> Display( qform );
Quadratic form
Gram Matrix:
<immutable 4x4 plist matrix over GF(2^4):
[[ Z(2^4)^3, Z(2)^0, 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), Z(2^2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), Z(2^4)^3, Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^4)^12 ]
]>
gap> mat2 := ToMatObj([[Z(16)^7,1,0,0],[0,0,0,0],
>              [0,0,Z(16)^2,1],[0,0,0,Z(16)^9]]*Z(16)^0, GF(16));
<4x4 plist matrix over GF(2^4)>
gap> qform2 := QuadraticFormByMatrix(mat2, GF(16));
< quadratic form >
gap> Display( qform2 );
Quadratic form
Gram Matrix:
<immutable 4x4 plist matrix over GF(2^4):
[[ Z(2^4)^7, Z(2)^0, 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), Z(2^4)^2, Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^4)^9 ]
]>
gap> biform := AssociatedBilinearForm( qform2 );
< bilinear form >
gap> Display( biform );
Bilinear form
Gram Matrix:
<immutable 4x4 plist matrix over GF(2^4):
[[ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ]
 [ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]
]>
gap> STOP_TEST("matobj/bg_th_ex8.tst", 10000 );
