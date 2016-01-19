#Constructing form: QuadraticFormByMatrix
mat := [[1,0,0,0],[0,3,0,0],[0,0,0,6],[0,0,6,0]]*Z(7)^0;
form := QuadraticFormByMatrix(mat,GF(7));
Display(form);
gf := GF(2^2);
mat := InvariantQuadraticForm( SO(-1, 4, 4) )!.matrix;
form := QuadraticFormByMatrix( mat, gf );
Display(form);
quit;
