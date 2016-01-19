#Constructing form: HermitianFormByPolynomial
r := PolynomialRing( GF(9), 4);
vars := IndeterminatesOfPolynomialRing( r );
poly := vars[1]*vars[2]^3+vars[1]^3*vars[2]+
             vars[3]*vars[4]^3+vars[3]^3*vars[4];
form := HermitianFormByPolynomial(poly,r);
Display(form);
quit;
