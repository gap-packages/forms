gap> START_TEST("Forms: bilformbyquadform.tst");
gap> r := PolynomialRing(GF(9),4);
GF(3^2)[x_1,x_2,x_3,x_4]
gap> poly := -r.1*r.2+Z(3^2)*r.3^2+r.4^2;
-x_1*x_2+Z(3^2)*x_3^2+x_4^2
gap> qform := QuadraticFormByPolynomial(poly,r);
< quadratic form >
gap> Display( qform );
Quadratic form
Gram Matrix:
z = Z(9)
   .   2   .   .
   .   .   .   .
   .   . z^1   .
   .   .   .   1
Polynomial: -x_1*x_2+Z(3^2)*x_3^2+x_4^2

gap> form := BilinearFormByQuadraticForm( qform );
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
z = Z(9)
   .   1   .   .
   1   .   .   .
   .   . z^1   .
   .   .   .   1
gap> Set(GF(9)^4, x -> [x,x]^form = x^qform);
[ true ]
gap> STOP_TEST("bilformbyquadform.tst", 10000 );
