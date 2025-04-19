# test_forms3: some forms in even char, parabolic example
f := GF(8);
mat := [[Z(8),0*Z(2),0*Z(2),0*Z(2),0*Z(2)],[0*Z(2),Z(2)^0,Z(2^3)^5,0*Z(2),0*Z(2)],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),0*Z(2)],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),Z(2)^0],[0*Z(2),0*Z(2),0*Z(2),0*Z(2),0*Z(2)]];
form := QuadraticFormByMatrix(mat,f);
TypeOfForm(form);
BaseChangeToCanonical(form);
iso := IsometricCanonicalForm(form);
Display(form);
Display(iso);
quit;
