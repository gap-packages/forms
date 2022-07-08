gap> START_TEST("Forms: test_forms7.tst");
gap> r := PolynomialRing( GF(7), 6);
GF(7)[x_1,x_2,x_3,x_4,x_5,x_6]
gap> vars := IndeterminatesOfPolynomialRing( r );
[ x_1, x_2, x_3, x_4, x_5, x_6 ]
gap> pol := (Z(7)^4)*vars[1]^2-vars[2]*vars[3]-vars[4]*vars[6]+vars[5]^2;
Z(7)^4*x_1^2-x_2*x_3-x_4*x_6+x_5^2
gap> form := BilinearFormByPolynomial(pol, r);
< bilinear form >
gap> IsEllipticForm(form);
true
gap> Display(form);
Elliptic bilinear form
Gram Matrix:
 4 . . . . .
 . . 3 . . .
 . 3 . . . .
 . . . . . 3
 . . . . 1 .
 . . . 3 . .
Polynomial: Z(7)^4*x_1^2-x_2*x_3-x_4*x_6+x_5^2

Witt Index: 2
gap> STOP_TEST("test_forms7.tst", 10000 );
