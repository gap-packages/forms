#Background theory: example 2
mat := [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]]*Z(9)^0;
form := BilinearFormByMatrix(mat,GF(9));
Display(form);
IsReflexiveForm(form);
IsSymmetricForm(form);
IsAlternatingForm(form);
r := RadicalOfForm(form);
Dimension(r);
quit;
