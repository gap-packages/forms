gap> START_TEST("Forms: bg_th_ex6.tst");
gap> V := GF(4)^3;                           
( GF(2^2)^3 )
gap> mat := [[Z(2^2)^2,Z(2^2),Z(2^2)^2],[Z(2^2)^2,Z(2)^0,Z(2)^0],
>         [0*Z(2),Z(2)^0,0*Z(2)]];
[ [ Z(2^2)^2, Z(2^2), Z(2^2)^2 ], [ Z(2^2)^2, Z(2)^0, Z(2)^0 ], 
  [ 0*Z(2), Z(2)^0, 0*Z(2) ] ]
gap> qform := QuadraticFormByMatrix(mat, GF(4));
< quadratic form >
gap> Display( qform );
Quadratic form
Gram Matrix:
z = Z(4)
 z^2   1 z^2
   .   1   .
   .   .   .
gap> PolynomialOfForm( qform );
Z(2^2)^2*x_1^2+x_1*x_2+Z(2^2)^2*x_1*x_3+x_2^2
gap> STOP_TEST("bg_th_ex6.tst", 10000 );
