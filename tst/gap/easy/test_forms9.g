#test_forms9: forms by polynomials
r := PolynomialRing( GF(16), 4);
vars := IndeterminatesOfPolynomialRing( r );
z := Z(16);
pol := z^5*vars[3]^2+vars[1]*vars[3]+z^8*vars[1]^2 + vars[2]*vars[4];
form := QuadraticFormByPolynomial(pol,r);
B := BaseChangeToCanonical(form);
mat := form!.matrix;
mat2 := B*mat*TransposedMat(B);
Display(mat2);
DiscriminantOfForm(form);
quit;
