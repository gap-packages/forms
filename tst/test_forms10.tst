gap> START_TEST("Forms: test_forms10.tst");
gap> r := PolynomialRing( GF(8), 6);
GF(2^3)[x_1,x_2,x_3,x_4,x_5,x_6]
gap> vars := IndeterminatesOfPolynomialRing( r );
[ x_1, x_2, x_3, x_4, x_5, x_6 ]
gap> pol := vars[1]*vars[6]+vars[2]*vars[5]+vars[3]*vars[4];
x_1*x_6+x_2*x_5+x_3*x_4
gap> form := QuadraticFormByPolynomial(pol,r,6);
< quadratic form >
gap> B := BaseChangeToCanonical(form);
< immutable compressed matrix 6x6 over GF(8) >
gap> mat := form!.matrix;;
gap> Display(mat);
 . . . . . 1
 . . . . 1 .
 . . . 1 . .
 . . . . . .
 . . . . . .
 . . . . . .
gap> mat2 := B*mat*TransposedMat(B);;
gap> Display(mat2);
 . 1 . . . .
 . . . . . .
 . . . 1 . .
 . . . . . .
 . . . . . .
 . . . . 1 .
gap> iso := IsometricCanonicalForm(form);
< hyperbolic quadratic form >
gap> Display(iso);
Hyperbolic quadratic form
Gram Matrix:
 . 1 . . . .
 . . . . . .
 . . . 1 . .
 . . . . . .
 . . . . . 1
 . . . . . .
Witt Index: 3
gap> STOP_TEST("test_forms10.tst", 10000 );
