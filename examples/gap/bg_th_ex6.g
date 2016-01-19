#Background theory: example 6
V := GF(4)^3;                           
mat := [[Z(2^2)^2,Z(2^2),Z(2^2)^2],[Z(2^2)^2,Z(2)^0,Z(2)^0],
        [0*Z(2),Z(2)^0,0*Z(2)]];
qform := QuadraticFormByMatrix(mat, GF(4));
Display( qform );
PolynomialOfForm( qform );
quit;
