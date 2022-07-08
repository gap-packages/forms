gap> START_TEST("Forms: test_forms5.tst");
gap> f := GF(9);
GF(3^2)
gap> gram := [[0,0,0,Z(9)^2],[0,0,1,0],[0,1,0,0],[-Z(9)^2,0,0,0]]*Z(9)^0;
[ [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3^2)^2 ], [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ], 
  [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], [ Z(3^2)^6, 0*Z(3), 0*Z(3), 0*Z(3) ] ]
gap> form := HermitianFormByMatrix(gram,f);
< hermitian form >
gap> TypeOfForm(form);
1/2
gap> BaseChangeToCanonical(form);
[ [ Z(3)^0, 0*Z(3), 0*Z(3), Z(3^2)^3 ], [ Z(3^2), 0*Z(3), 0*Z(3), Z(3^2)^2 ], 
  [ 0*Z(3), Z(3^2), Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3^2)^2, Z(3^2)^3, 0*Z(3) ] ]
gap> iso := IsometricCanonicalForm(form);
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
z = Z(9)
   .   .   . z^2
   .   .   1   .
   .   1   .   .
 z^6   .   .   .
Witt Index: 2
gap> Display(iso);
Hermitian form
Gram Matrix:
 1 . . .
 . 1 . .
 . . 1 .
 . . . 1
Witt Index: 2
gap> STOP_TEST("test_forms5.tst", 10000 );
