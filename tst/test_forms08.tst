gap> START_TEST("Forms: test_forms8.tst");
gap> r := PolynomialRing( GF(8), 3);
GF(2^3)[x_1,x_2,x_3]
gap> vars := IndeterminatesOfPolynomialRing( r );
[ x_1, x_2, x_3 ]
gap> pol := vars[1]^2 + vars[3]^2 + vars[2]^2;
x_1^2+x_2^2+x_3^2
gap> form := QuadraticFormByPolynomial(pol, r);
< quadratic form >
gap> BaseChangeToCanonical(form);
[ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ Z(2)^0, Z(2)^0, 0*Z(2) ], 
  [ Z(2)^0, 0*Z(2), Z(2)^0 ] ]
gap> IsDegenerateForm(form);
#I  Testing degeneracy of the *associated bilinear form*
true
gap> RadicalOfForm(form);
<vector space over GF(2^3), with 63 generators>
gap> STOP_TEST("test_forms8.tst", 10000 );
