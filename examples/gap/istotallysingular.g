#testing singularity for subspaces.
mat := [[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0;
form := QuadraticFormByMatrix(mat);
sub := [[Z(2)^0,0*Z(2),Z(2^3)^6,Z(2^3),Z(2^3)^3],
       [0*Z(2),Z(2)^0,Z(2^3)^6,Z(2^3)^2,Z(2^3)]];
IsTotallySingularSubspace(form,sub);
quit;
