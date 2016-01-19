#Constructing form: BilinearFormByPolynomial
r := PolynomialRing( GF(11), 4);
vars := IndeterminatesOfPolynomialRing( r );
pol := vars[1]*vars[2]+vars[3]*vars[4];
form := BilinearFormByPolynomial(pol, r, 4);
Display(form);
r := PolynomialRing(GF(4),2);
pol := r.1*r.2;
form := BilinearFormByPolynomial(pol,r);
quit;
quit;
