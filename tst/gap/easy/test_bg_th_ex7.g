#Background theory: example 7
r := PolynomialRing(GF(8),4);
poly := r.1*r.2+r.3*r.4;
qform := QuadraticFormByPolynomial(poly, r);
Display(qform);
r := RadicalOfForm(qform);;
Dimension(r);
quit;
