gap> START_TEST("Forms: test_forms4.tst");
gap> f := GF(8);
GF(2^3)
gap> mat := [[Z(8),0,0,0],[0,0,Z(8)^4,0],[0,0,0,1],[0,0,0,0]]*Z(8)^0;
[ [ Z(2^3), 0*Z(2), 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), Z(2^3)^4, 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
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
z = Z(8)
 z^1   .   .   .
   .   . z^4   .
   .   .   .   1
   .   .   .   .
Witt Index: 1
gap> Display(iso);
Parabolic quadratic form
Gram Matrix:
 1 . . .
 . . 1 .
 . . . .
 . . . .
Witt Index: 1
gap> IsDegenerateForm(iso);
#I  Testing degeneracy of the *associated bilinear form*
true
gap> RadicalOfForm(iso);
<vector space over GF(2^3), with 7 generators>
gap> STOP_TEST("test_forms4.tst", 10000 );
