#Constructing form: BilinearFormByMatrix
mat := IdentityMat(4, GF(9));
form := BilinearFormByMatrix(mat,GF(9));
Display(form);
mat := [[0*Z(2),Z(16)^12,0*Z(2),Z(4)^2,Z(16)^13],
   [Z(16)^12,0*Z(2),0*Z(2),Z(16)^11,Z(16)],
   [0*Z(2),0*Z(2),0*Z(2),Z(4)^2,Z(16)^3],
   [Z(4)^2,Z(16)^11,Z(4)^2,0*Z(2),Z(16)^3],
   [Z(16)^13,Z(16),Z(16)^3,Z(16)^3,0*Z(2) ]];
form := BilinearFormByMatrix(mat,GF(16));
Display(form);
mat := [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]*Z(7)^0;
form := BilinearFormByMatrix(mat);
WittIndex(form);
form := BilinearFormByMatrix(mat,GF(49));
WittIndex(form);
quit;
