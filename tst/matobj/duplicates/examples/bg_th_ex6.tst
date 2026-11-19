gap> START_TEST("Forms: matobj/bg_th_ex6.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> V := GF(4)^3;                           
( GF(2^2)^3 )
gap> mat := ToMatObj([[Z(2^2)^2,Z(2^2),Z(2^2)^2],[Z(2^2)^2,Z(2)^0,Z(2)^0],
>         [0*Z(2),Z(2)^0,0*Z(2)]], GF(4));
<3x3 plist matrix over GF(2^2)>
gap> qform := QuadraticFormByMatrix(mat, GF(4));
< quadratic form >
gap> Display( qform );
Quadratic form
Gram Matrix:
<immutable 3x3 plist matrix over GF(2^2):
[[ Z(2^2)^2, Z(2)^0, Z(2^2)^2 ]
 [ 0*Z(2), Z(2)^0, 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2) ]
]>
gap> PolynomialOfForm( qform );
Z(2^2)^2*x_1^2+x_1*x_2+Z(2^2)^2*x_1*x_3+x_2^2
gap> STOP_TEST("matobj/bg_th_ex6.tst", 10000 );
