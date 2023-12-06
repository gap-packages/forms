gap> START_TEST("Forms: trivialform2.tst");
gap> 
gap> #Behaviour of a trivial form
gap> mat := [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]*Z(3)^0;
[ [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ] ]
gap> form := BilinearFormByMatrix(mat,GF(3));
< trivial form >
gap> v := Random(GF(3)^4);
[ Z(3), Z(3), 0*Z(3), Z(3) ]
gap> [v,v]^form;
0*Z(3)
gap> v^form;
0*Z(3)
gap> STOP_TEST("trivialform2.tst", 10000 );
