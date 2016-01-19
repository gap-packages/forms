#Radical of form: RadicalOfForm
r := PolynomialRing( GF(8), 3 );
poly := r.1^2 + r.2 * r.3;
form := QuadraticFormByPolynomial( poly, r );
r := RadicalOfForm( form );
Dimension(r);
quit;
