#test_forms8: forms by polynomials
r := PolynomialRing( GF(8), 3);
vars := IndeterminatesOfPolynomialRing( r );
pol := vars[1]^2 + vars[3]^2 + vars[2]^2;
form := QuadraticFormByPolynomial(pol, r);
BaseChangeToCanonical(form);
IsDegenerateForm(form);
RadicalOfForm(form);
quit;
