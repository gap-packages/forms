gap> START_TEST("Forms: matobj/test_forms4.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> f := GF(8);
GF(2^3)
gap> mat := ToMatObj([[Z(8),0,0,0],[0,0,Z(8)^4,0],[0,0,0,1],[0,0,0,0]]*Z(8)^0, f);
<4x4 plist matrix over GF(2^3)>
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> IsSingularForm(form);
true
gap> TypeOfForm(form);
0
gap> iso := IsometricCanonicalForm(form);
< parabolic quadratic form >
gap> Display(form);
Singular parabolic quadratic form
Gram Matrix:
<immutable 4x4 plist matrix over GF(2^3):
[[ Z(2^3), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), Z(2^3)^4, 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
]>
Witt Index: 1
gap> Display(iso);
Parabolic quadratic form
Gram Matrix:
<immutable 4x4 plist matrix over GF(2^3):
[[ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
]>
Witt Index: 1
gap> IsDegenerateForm(iso);
#I  Testing degeneracy of the *associated bilinear form*
true
gap> RadicalOfForm(iso);
<vector space over GF(2^3), with 7 generators>
gap> STOP_TEST("matobj/test_forms4.tst", 10000 );
