gap> START_TEST("Forms: bg_th_ex8.tst");
gap> #Background theory: example 8
gap> mat := [[Z(16)^3,1,0,0],[0,Z(16)^5,0,0],
>              [0,0,Z(16)^3,1],[0,0,0,Z(16)^12]]*Z(16)^0;
[ [ Z(2^4)^3, Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), Z(2^2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2^4)^3, Z(2)^0 ], [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^4)^12 ] 
 ]
gap> qform := QuadraticFormByMatrix(mat,GF(16));
< quadratic form >
gap> Display( qform );
Quadratic form
Gram Matrix:
z = Z(16)
  z^3    1    .    .
    .  z^5    .    .
    .    .  z^3    1
    .    .    . z^12
gap> mat2 := [[Z(16)^7,1,0,0],[0,0,0,0],
>              [0,0,Z(16)^2,1],[0,0,0,Z(16)^9]]*Z(16)^0;
[ [ Z(2^4)^7, Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2^4)^2, Z(2)^0 ], [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^4)^9 ] ]
gap> qform2 := QuadraticFormByMatrix(mat2, GF(16));
< quadratic form >
gap> Display( qform2 );
Quadratic form
Gram Matrix:
z = Z(16)
  z^7    1    .    .
    .    .    .    .
    .    .  z^2    1
    .    .    .  z^9
gap> biform := AssociatedBilinearForm( qform2 );
< bilinear form >
gap> Display( biform );
Bilinear form
Gram Matrix:
 . 1 . .
 1 . . .
 . . . 1
 . . 1 .
gap> STOP_TEST("bg_th_ex8.tst", 10000 );
