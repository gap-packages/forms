gap> START_TEST("Forms: bilformbypoly.tst");
gap> #Constructing form: BilinearFormByPolynomial
gap> r := PolynomialRing( GF(11), 4);
GF(11)[x_1,x_2,x_3,x_4]
gap> vars := IndeterminatesOfPolynomialRing( r );
[ x_1, x_2, x_3, x_4 ]
gap> pol := vars[1]*vars[2]+vars[3]*vars[4];
x_1*x_2+x_3*x_4
gap> form := BilinearFormByPolynomial(pol, r, 4);
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
  .  6  .  .
  6  .  .  .
  .  .  .  6
  .  .  6  .
Polynomial: x_1*x_2+x_3*x_4

gap> r := PolynomialRing(GF(4),2);
GF(2^2)[x_1,x_2]
gap> pol := r.1*r.2;
x_1*x_2
gap> form := BilinearFormByPolynomial(pol,r);
Error, No orthogonal form can be associated with a quadratic polynomial in eve\
n characteristic
gap> STOP_TEST("bilformbypoly.tst", 10000 );
