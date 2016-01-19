#Constructing form: BilinearFormByQuadraticForm
r := PolynomialRing(GF(9),4);
poly := -r.1*r.2+Z(3^2)*r.3^2+r.4^2;
qform := QuadraticFormByPolynomial(poly,r);
Display( qform );
form := BilinearFormByQuadraticForm( qform );
Display(form);
Set(GF(9)^4, x -> [x,x]^form = x^qform);
quit;
