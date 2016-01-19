#Constructing form: AssociatedBilinearForm
r:= PolynomialRing(GF(121),6);
poly := r.1*r.5-r.2*r.6+r.3*r.4;
form := QuadraticFormByPolynomial(poly,r);
aform := AssociatedBilinearForm(form);
Display(aform);
quit;
