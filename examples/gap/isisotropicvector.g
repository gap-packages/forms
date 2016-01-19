#testing isotropy for vectors.
mat := [[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]]*Z(41)^0;
form := BilinearFormByMatrix(mat);
v := [1,1,0,0]*Z(41)^0;
IsIsotropicVector(form,v);
mat := [[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0;
form := QuadraticFormByMatrix(mat);
v1 := [1,0,0,0,0]*Z(8)^0;
v2 := [0,1,0,0,0]*Z(8)^0;
IsIsotropicVector(form,v1);
IsIsotropicVector(form,v2);
quit;
