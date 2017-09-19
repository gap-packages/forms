#computing subspace orthogonal to given vector
mat := [[0,0,0,-2],[0,0,-3,0],[0,3,0,0],[2,0,0,0]]*Z(7)^0;
form := BilinearFormByMatrix(mat);
v := [0*Z(7),Z(7)^0,Z(7)^3,Z(7)^5];
vperp := OrthogonalSubspaceMat(form,v);
List(vperp,x->[x,v]^form);
sub := [[1,1,0,0],[0,0,1,2]]*Z(7)^0;
subperp := OrthogonalSubspaceMat(form,sub);
List(subperp,x->List(sub,y->[x,y]^form));
mat := [[1,0,0],[0,0,1],[0,0,0]]*Z(2)^0;
form := QuadraticFormByMatrix(mat);
v := [Z(2)^0,Z(2)^0,0*Z(2)];
vperp := OrthogonalSubspaceMat(form,v);
bil_form := AssociatedBilinearForm(form);
List(vperp,x->[x,v]^bil_form);
sub := [[1,0,1],[1,0,0]]*Z(2)^0;
subperp := OrthogonalSubspaceMat(form,sub);
List(subperp,x->List(sub,y->[x,y]^bil_form));
quit;
