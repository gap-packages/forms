gap> START_TEST("Forms: test_forms6.tst");
gap> r := PolynomialRing( GF(11), 4);
GF(11)[x_1,x_2,x_3,x_4]
gap> vars := IndeterminatesOfPolynomialRing( r );
[ x_1, x_2, x_3, x_4 ]
gap> pol := vars[1]*vars[2]+vars[3]*vars[4];
x_1*x_2+x_3*x_4
gap> form := BilinearFormByPolynomial(pol, r, 4);
< bilinear form >
gap> BaseChangeToCanonical(form);
[ [ Z(11)^2, 0*Z(11), Z(11)^2, 0*Z(11) ], 
  [ 0*Z(11), Z(11)^8, Z(11)^7, 0*Z(11) ], 
  [ Z(11)^4, Z(11)^5, Z(11)^4, Z(11)^0 ], 
  [ 0*Z(11), 0*Z(11), Z(11)^0, 0*Z(11) ] ]
gap> TypeOfForm(form);
1
gap> STOP_TEST("test_forms6.tst", 10000 );
