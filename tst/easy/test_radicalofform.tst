gap> START_TEST("Forms: test_radicalofform.tst");
gap> r := PolynomialRing( GF(8), 3 );
GF(2^3)[x_1,x_2,x_3]
gap> poly := r.1^2 + r.2 * r.3;
x_1^2+x_2*x_3
gap> form := QuadraticFormByPolynomial( poly, r );
< quadratic form >
gap> r := RadicalOfForm( form );;
gap> Dimension(r);
0
gap> STOP_TEST("test_radicalofform.tst", 10000 );
