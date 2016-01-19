#computing subspace orthogonal to given vector
mat := [[0,0,0,-2],[0,0,-3,0],[0,3,0,0],[2,0,0,0]]*Z(7)^0;
form := BilinearFormByMatrix(mat);
v := Random(GF(7)^4);
vperp := OrthogonalSubspaceMat(form,v);
sub := [[1,1,0,0],[0,0,1,2]]*Z(7)^0;
subperp := OrthogonalSubspaceMat(form,sub);
mat := [[1,0,0],[0,0,1],[0,0,0]]*Z(2)^0;
form := QuadraticFormByMatrix(mat);
v := Random(GF(2)^3);
vperp := OrthogonalSubspaceMat(form,v);
sub := [[1,0,1],[1,0,0]]*Z(2)^0;
subperp := OrthogonalSubspaceMat(form,sub);
quit;
