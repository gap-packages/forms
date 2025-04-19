#test_forms4
f := GF(8);
mat := [[Z(8),0,0,0],[0,0,Z(8)^4,0],[0,0,0,1],[0,0,0,0]]*Z(8)^0;
form := QuadraticFormByMatrix(mat,f);
IsSingularForm(form);
TypeOfForm(form);
iso := IsometricCanonicalForm(form);
Display(form);
Display(iso);
IsDegenerateForm(iso);
RadicalOfForm(iso);
quit;
