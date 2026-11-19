gap> START_TEST("Forms: matobj/test_forms3.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> f := GF(8);
GF(2^3)
gap> mat := ToMatObj([[Z(8),0*Z(2),0*Z(2),0*Z(2),0*Z(2)],[0*Z(2),Z(2)^0,Z(2^3)^5,0*Z(2),0*Z(2)],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),0*Z(2)],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),Z(2)^0],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),0*Z(2)]], f);
<5x5 plist matrix over GF(2^3)>
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> TypeOfForm(form);
0
gap> BaseChangeToCanonical(form);
<immutable 5x5 plist matrix over GF(2^3)>
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
gap> STOP_TEST("matobj/test_forms3.tst", 10000 );
