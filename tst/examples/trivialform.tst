gap> START_TEST("Forms: trivialform.tst");
gap> mat := [[0,0,0],[0,0,0],[0,0,0]]*Z(7)^0;
[ [ 0*Z(7), 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7) ] ]
gap> form1 := BilinearFormByMatrix(mat,GF(7));
< trivial form >
gap> form2 := QuadraticFormByMatrix(mat,GF(7));
< trivial form >
gap> form1 = form2;
true
gap> IsQuadraticForm(form1);
false
gap> IsSesquilinearForm(form1);
false
gap> mat := [[0,0],[0,0]]*Z(4)^0;
[ [ 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2) ] ]
gap> form3 := BilinearFormByMatrix(mat,GF(4));
< trivial form >
gap> form3 = form1;
false
gap> STOP_TEST("trivialform.tst", 10000 );
