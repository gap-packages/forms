# Witt index. Also of degenerated forms
mat := [[0,0,1,0,0],[0,0,0,0,0],[-1,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]*Z(7)^0;
form := BilinearFormByMatrix(mat,GF(7));
WittIndex(form);
RadicalOfForm(form);
Dimension(last);
mat := IdentityMat(6,GF(5));
form := QuadraticFormByMatrix(mat,GF(5));
WittIndex(form);
mat := IdentityMat(6,GF(7));
form := QuadraticFormByMatrix(mat,GF(7));
WittIndex(form);
quit;

