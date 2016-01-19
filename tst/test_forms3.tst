gap> START_TEST("Forms: test_forms3.tst");
gap> f := GF(8);
GF(2^3)
gap> mat := [[Z(8),0*Z(2),0*Z(2),0*Z(2),0*Z(2)],[0*Z(2),Z(2)^0,Z(2^3)^5,0*Z(2),0*Z(2)],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),0*Z(2)],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),Z(2)^0],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),0*Z(2)]];
[ [ Z(2^3), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), Z(2)^0, Z(2^3)^5, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> TypeOfForm(form);
0
gap> BaseChangeToCanonical(form);
[ [ Z(2^3)^3, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), Z(2^3)^2, Z(2^3)^4, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ] ]
gap> iso := IsometricCanonicalForm(form);
< parabolic quadratic form >
gap> Display(form);
Parabolic quadratic form
Gram Matrix:
z = Z(8)
 z^1   .   .   .   .
   .   1 z^5   .   .
   .   .   .   .   .
   .   .   .   .   1
   .   .   .   .   .
Witt Index: 2
gap> Display(iso);
Parabolic quadratic form
Gram Matrix:
 1 . . . .
 . . 1 . .
 . . . . .
 . . . . 1
 . . . . .
Witt Index: 2
gap> STOP_TEST("test_forms3.tst", 10000 );
