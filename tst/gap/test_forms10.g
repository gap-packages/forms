#Klein's quadric in PG(5,8)
r := PolynomialRing( GF(8), 6);
vars := IndeterminatesOfPolynomialRing( r );
pol := vars[1]*vars[6]+vars[2]*vars[5]+vars[3]*vars[4];
form := QuadraticFormByPolynomial(pol,r,6);
B := BaseChangeToCanonical(form);
mat := form!.matrix;;
Display(mat);
mat2 := B*mat*TransposedMat(B);;
Display(mat2);
iso := IsometricCanonicalForm(form);
Display(iso);
quit;
