#testing total isotropy for subspaces
mat := [[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]]*Z(7)^0;
form := BilinearFormByMatrix(mat);
sub:= [[Z(7)^0,0*Z(7),Z(7)^0,Z(7)],[0*Z(7),Z(7)^0,Z(7)^0,Z(7)^4]];
IsTotallyIsotropicSubspace(form,sub);
mat := IdentityMat(6,GF(2));
form := HermitianFormByMatrix(mat,GF(4));
sub := [[Z(2)^0,0*Z(2),0*Z(2),Z(2)^0,Z(2)^0,Z(2)^0], 
  [0*Z(2),Z(2)^0,0*Z(2),Z(2^2)^2,Z(2^2),Z(2)^0], 
  [0*Z(2),0*Z(2),Z(2)^0,Z(2)^0,Z(2^2),Z(2^2)^2]];
IsTotallyIsotropicSubspace(form,sub);
quit;
