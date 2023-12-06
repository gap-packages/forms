gap> START_TEST("Forms: bilformbymatrix.tst");
gap> #Constructing form: BilinearFormByMatrix
gap> mat := IdentityMat(4, GF(9));
[ [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ] ]
gap> form := BilinearFormByMatrix(mat,GF(9));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
 1 . . .
 . 1 . .
 . . 1 .
 . . . 1
gap> mat := [[0*Z(2),Z(16)^12,0*Z(2),Z(4)^2,Z(16)^13],
>    [Z(16)^12,0*Z(2),0*Z(2),Z(16)^11,Z(16)],
>    [0*Z(2),0*Z(2),0*Z(2),Z(4)^2,Z(16)^3],
>    [Z(4)^2,Z(16)^11,Z(4)^2,0*Z(2),Z(16)^3],
>    [Z(16)^13,Z(16),Z(16)^3,Z(16)^3,0*Z(2) ]];
[ [ 0*Z(2), Z(2^4)^12, 0*Z(2), Z(2^2)^2, Z(2^4)^13 ], 
  [ Z(2^4)^12, 0*Z(2), 0*Z(2), Z(2^4)^11, Z(2^4) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^2)^2, Z(2^4)^3 ], 
  [ Z(2^2)^2, Z(2^4)^11, Z(2^2)^2, 0*Z(2), Z(2^4)^3 ], 
  [ Z(2^4)^13, Z(2^4), Z(2^4)^3, Z(2^4)^3, 0*Z(2) ] ]
gap> form := BilinearFormByMatrix(mat,GF(16));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
z = Z(16)
    . z^12    . z^10 z^13
 z^12    .    . z^11  z^1
    .    .    . z^10  z^3
 z^10 z^11 z^10    .  z^3
 z^13  z^1  z^3  z^3    .
gap> mat := [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]*Z(7)^0;
[ [ Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ], [ 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7) ] ]
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> WittIndex(form);
1
gap> form := BilinearFormByMatrix(mat,GF(49));
< bilinear form >
gap> WittIndex(form);
2
gap> STOP_TEST("bilformbymatrix.tst", 10000 );
