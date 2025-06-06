gap> START_TEST("Forms: hermitianformbypoly.tst");
gap> r := PolynomialRing( GF(9), 4);
GF(3^2)[x_1,x_2,x_3,x_4]
gap> vars := IndeterminatesOfPolynomialRing( r );
[ x_1, x_2, x_3, x_4 ]
gap> poly := vars[1]*vars[2]^3+vars[1]^3*vars[2]+
>              vars[3]*vars[4]^3+vars[3]^3*vars[4];
x_1^3*x_2+x_1*x_2^3+x_3^3*x_4+x_3*x_4^3
gap> form := HermitianFormByPolynomial(poly,r);
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
 . 1 . .
 1 . . .
 . . . 1
 . . 1 .
Polynomial: x_1^3*x_2+x_1*x_2^3+x_3^3*x_4+x_3*x_4^3

gap> STOP_TEST("hermitianformbypoly.tst", 10000 );
