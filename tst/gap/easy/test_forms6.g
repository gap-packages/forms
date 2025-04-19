#test_forms6: forms by polynomials.
r := PolynomialRing( GF(11), 4);
vars := IndeterminatesOfPolynomialRing( r );
pol := vars[1]*vars[2]+vars[3]*vars[4];
form := BilinearFormByPolynomial(pol, r, 4);
BaseChangeToCanonical(form);
TypeOfForm(form);
quit;
