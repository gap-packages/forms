#testing singularity for vectors.
mat := [[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0;
form := QuadraticFormByMatrix(mat);
v1 := [1,0,0,0,0]*Z(8)^0;
v2 := [0,1,0,0,0]*Z(8)^0;
IsSingularVector(form,v1);
IsSingularVector(form,v2);
IsIsotropicVector(form,v1);
IsIsotropicVector(form,v2);
quit;
