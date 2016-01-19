#Constructing form: QuadraticFormByPolynomial
r := PolynomialRing( GF(8), 3);
poly := r.1^2 + r.2^2 + r.3^2;
form := QuadraticFormByPolynomial(poly, r);
RadicalOfForm(form);
r := PolynomialRing(GF(9),4);
poly := Z(3)^2*r.1^2+r.2^2+r.3*r.4;
qform := QuadraticFormByPolynomial(poly,r);
Display(qform);
quit;
