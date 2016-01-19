#Background theory: example 3
mat := [[0,0,-2],[0,0,1],[2,-1,0]]*Z(7)^0;
form := BilinearFormByMatrix(mat,GF(7));
Display(form);
IsSymmetricForm(form);
IsAlternatingForm(form);
r := RadicalOfForm(form);
Dimension(r);
quit;
