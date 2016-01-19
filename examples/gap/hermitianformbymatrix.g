#Constructing form: HermitianFormByMatrix
gf := GF(3^2);
mat := IdentityMat(4, gf);
form := HermitianFormByMatrix( mat, gf );
Display(form);
mat := [[Z(11)^0,0*Z(11),0*Z(11)],[0*Z(11),0*Z(11),Z(11)],
    [0*Z(11),Z(11),0*Z(11)]];
form := HermitianFormByMatrix(mat,GF(121));
Display(form);
quit;
