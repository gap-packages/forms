gap> START_TEST("Forms: isometriccanonicalform.tst");
gap> mat := [ [ Z(8) , 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
> [ 0*Z(2), Z(2)^0, Z(2^3)^5, 0*Z(2), 0*Z(2) ], 
> [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
> [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
> [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ];;
gap> form := QuadraticFormByMatrix(mat,GF(8));
< quadratic form >
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
gap> STOP_TEST("isometriccanonicalform.tst", 10000 );
