gap> START_TEST("Forms: quadformbymatrix.tst");
gap> mat := [[1,0,0,0],[0,3,0,0],[0,0,0,6],[0,0,6,0]]*Z(7)^0;
[ [ Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^3 ], [ 0*Z(7), 0*Z(7), Z(7)^3, 0*Z(7) ] ]
gap> form := QuadraticFormByMatrix(mat,GF(7));
< quadratic form >
gap> Display(form);
Quadratic form
Gram Matrix:
 1 . . .
 . 3 . .
 . . . 5
 . . . .
gap> gf := GF(2^2);
GF(2^2)
gap> mat := InvariantQuadraticForm( SO(-1, 4, 4) )!.matrix;
[ [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2^2)^2, Z(2)^0 ], [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^2)^2 ] ]
gap> form := QuadraticFormByMatrix( mat, gf );
< quadratic form >
gap> Display(form);
Quadratic form
Gram Matrix:
z = Z(4)
   .   1   .   .
   .   .   .   .
   .   . z^2   1
   .   .   . z^2
gap> STOP_TEST("quadformbymatrix.tst", 10000 );
