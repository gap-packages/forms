gap> START_TEST("Forms: quadformbypoly.tst");
gap> r := PolynomialRing( GF(8), 3);
GF(2^3)[x_1,x_2,x_3]
gap> poly := r.1^2 + r.2^2 + r.3^2;
x_1^2+x_2^2+x_3^2
gap> form := QuadraticFormByPolynomial(poly, r);
< quadratic form >
gap> RadicalOfForm(form);
<vector space over GF(2^3), with 63 generators>
gap> r := PolynomialRing(GF(9),4);
GF(3^2)[x_1,x_2,x_3,x_4]
gap> poly := Z(3)^2*r.1^2+r.2^2+r.3*r.4;
x_1^2+x_2^2+x_3*x_4
gap> qform := QuadraticFormByPolynomial(poly,r);
< quadratic form >
gap> Display(qform);
Quadratic form
Gram Matrix:
 1 . . .
 . 1 . .
 . . . 1
 . . . .
Polynomial: x_1^2+x_2^2+x_3*x_4

gap> STOP_TEST("quadformbypoly.tst", 10000 );
