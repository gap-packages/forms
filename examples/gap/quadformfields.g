#constructing the same form over different fields
mat := 
[[Z(2)^0,Z(2)^0,0*Z(2),0*Z(2)],[0*Z(2),Z(2)^0,0*Z(2),0*Z(2)], 
 [0*Z(2),0*Z(2),0*Z(2),Z(2)^0],[0*Z(2),0*Z(2),0*Z(2),0*Z(2)]];
form := QuadraticFormByMatrix(mat);
WittIndex(form);
form := QuadraticFormByMatrix(mat,GF(4));
WittIndex(form);
quit;
