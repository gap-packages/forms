# test_forms5: hermitian forms
f := GF(9);
gram := [[0,0,0,Z(9)^2],[0,0,1,0],[0,1,0,0],[-Z(9)^2,0,0,0]]*Z(9)^0;
form := HermitianFormByMatrix(gram,f);
TypeOfForm(form);
BaseChangeToCanonical(form);
iso := IsometricCanonicalForm(form);
Display(form);
Display(iso);
quit;
