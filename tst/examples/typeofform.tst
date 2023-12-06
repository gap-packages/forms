gap> START_TEST("Forms: typeofform.tst");
gap> #Type Of Form
gap> mat := [[0,0,0,-1],[0,0,3,0],[0,-3,0,0],[1,0,0,0]]*Z(25)^0;
[ [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^2 ], [ 0*Z(5), 0*Z(5), Z(5)^3, 0*Z(5) ], 
  [ 0*Z(5), Z(5), 0*Z(5), 0*Z(5) ], [ Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5) ] ]
gap> form := BilinearFormByMatrix(mat,GF(25));
< bilinear form >
gap> IsDegenerateForm(form);
false
gap> TypeOfForm(form);
0
gap> mat := IdentityMat(3,GF(7));
[ [ Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^0, 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), Z(7)^0 ] ]
gap> form := QuadraticFormByMatrix(mat,GF(7));
< quadratic form >
gap> IsSingularForm(form);
false
gap> TypeOfForm(form);
0
gap> mat := [[0,1,0,0],[-1,0,0,0],[0,0,0,0],[0,0,0,0]]*Z(5)^0;
[ [ 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5) ], [ Z(5)^2, 0*Z(5), 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ] ]
gap> form := BilinearFormByMatrix(mat,GF(5));
< bilinear form >
gap> IsDegenerateForm(form);
true
gap> TypeOfForm(form);
0
gap> mat := [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]*Z(7)^0;
[ [ Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ], [ 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7) ] ]
gap> form := BilinearFormByMatrix(mat,GF(7));
< bilinear form >
gap> IsDegenerateForm(form);
false
gap> TypeOfForm(form);
-1
gap> mat := IdentityMat(3,GF(9));
[ [ Z(3)^0, 0*Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0, 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), Z(3)^0 ] ]
gap> form := HermitianFormByMatrix(mat,GF(9));
< hermitian form >
gap> IsDegenerateForm(form);
false
gap> TypeOfForm(form);
-1/2
gap> mat := [[0,0,0,1],[0,1,0,0],[0,0,1,0],[1,0,0,0]]*Z(8)^0;
[ [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], [ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := BilinearFormByMatrix(mat,GF(8));
< bilinear form >
gap> IsDegenerateForm(form);
false
gap> TypeOfForm(form);
Error, <f> is a pseudo form and has no defined type
gap> STOP_TEST("typeofform.tst", 10000 );
