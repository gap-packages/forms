gap> START_TEST("Forms: assocbilform.tst");
gap> #Constructing form: AssociatedBilinearForm
gap> r:= PolynomialRing(GF(121),6);
GF(11^2)[x_1,x_2,x_3,x_4,x_5,x_6]
gap> poly := r.1*r.5-r.2*r.6+r.3*r.4;
x_1*x_5-x_2*x_6+x_3*x_4
gap> form := QuadraticFormByPolynomial(poly,r);
< quadratic form >
gap> aform := AssociatedBilinearForm(form);
< bilinear form >
gap> Display(aform);
Bilinear form
Gram Matrix:
  .  .  .  .  1  .
  .  .  .  .  . 10
  .  .  .  1  .  .
  .  .  1  .  .  .
  1  .  .  .  .  .
  . 10  .  .  .  .
gap> STOP_TEST("assocbilform.tst", 10000 );
