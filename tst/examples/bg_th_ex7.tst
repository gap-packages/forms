gap> START_TEST("Forms: bg_th_ex7.tst");
gap> r := PolynomialRing(GF(8),4);
GF(2^3)[x_1,x_2,x_3,x_4]
gap> poly := r.1*r.2+r.3*r.4;
x_1*x_2+x_3*x_4
gap> qform := QuadraticFormByPolynomial(poly, r);
< quadratic form >
gap> Display(qform);
Quadratic form
Gram Matrix:
 . 1 . .
 . . . .
 . . . 1
 . . . .
Polynomial: x_1*x_2+x_3*x_4

gap> RadicalOfForm(qform);
<vector space of dimension 0 over GF(2^3)>
gap> STOP_TEST("bg_th_ex7.tst", 10000 );
