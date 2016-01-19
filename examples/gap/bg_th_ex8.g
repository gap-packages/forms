#Background theory: example 8
mat := [[Z(16)^3,1,0,0],[0,Z(16)^5,0,0],
             [0,0,Z(16)^3,1],[0,0,0,Z(16)^12]]*Z(16)^0;
qform := QuadraticFormByMatrix(mat,GF(16));
Display( qform );
mat2 := [[Z(16)^7,1,0,0],[0,0,0,0],
             [0,0,Z(16)^2,1],[0,0,0,Z(16)^9]]*Z(16)^0;
qform2 := QuadraticFormByMatrix(mat2, GF(16));
Display( qform2 );
biform := AssociatedBilinearForm( qform2 );
Display( biform );
quit;
