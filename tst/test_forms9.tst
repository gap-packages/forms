gap> START_TEST("Forms: test_forms9.tst");
gap> r := PolynomialRing( GF(16), 4);
GF(2^4)[x_1,x_2,x_3,x_4]
gap> vars := IndeterminatesOfPolynomialRing( r );
[ x_1, x_2, x_3, x_4 ]
gap> z := Z(16);
Z(2^4)
gap> pol := z^5*vars[3]^2+vars[1]*vars[3]+z^8*vars[1]^2 + vars[2]*vars[4];
Z(2^4)^8*x_1^2+x_1*x_3+x_2*x_4+Z(2^2)*x_3^2
gap> form := QuadraticFormByPolynomial(pol,r);
< quadratic form >
gap> B := BaseChangeToCanonical(form);
[ [ Z(2^4)^11, 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ Z(2^4)^7, 0*Z(2), Z(2^4)^4, 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ] ]
gap> mat := form!.matrix;
[ [ Z(2^4)^8, 0*Z(2), Z(2)^0, 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), Z(2^2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> mat2 := B*mat*TransposedMat(B);
[ [ Z(2)^0, Z(2^4)^12, 0*Z(2), 0*Z(2) ], 
  [ Z(2^4)^11, Z(2^4)^3, 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ],
  [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ] ]
gap> Display(mat2);
z = Z(16)
    1 z^12    .    .
 z^11  z^3    .    .
    .    .    .    .
    .    .    1    .
gap> DiscriminantOfForm(form);
"square"
gap> STOP_TEST("test_forms9.tst", 10000 );
