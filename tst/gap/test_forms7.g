#test_forms7: forms by polynomials
r := PolynomialRing( GF(7), 6);
vars := IndeterminatesOfPolynomialRing( r );
pol := (Z(7)^4)*vars[1]^2-vars[2]*vars[3]-vars[4]*vars[6]+vars[5]^2;
form := BilinearFormByPolynomial(pol, r);
IsEllipticForm(form);
Display(form);
quit;
