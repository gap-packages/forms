gap> START_TEST("Forms: matobj/isometriccanonicalform.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([ [ Z(8) , 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
> [ 0*Z(2), Z(2)^0, Z(2^3)^5, 0*Z(2), 0*Z(2) ], 
> [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
> [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
> [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ], GF(8));;
gap> form := QuadraticFormByMatrix(mat,GF(8));
< quadratic form >
gap> iso := IsometricCanonicalForm(form);
< parabolic quadratic form >
gap> Display(form);
Parabolic quadratic form
Gram Matrix:
<immutable 5x5 plist matrix over GF(2^3):
[[ Z(2^3), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), Z(2)^0, Z(2^3)^5, 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
]>
Witt Index: 2
gap> Display(iso);
Parabolic quadratic form
Gram Matrix:
<immutable 5x5 plist matrix over GF(2^3):
[[ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
]>
Witt Index: 2
gap> STOP_TEST("matobj/isometriccanonicalform.tst", 10000 );
