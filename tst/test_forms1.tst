gap> START_TEST("Forms: test_forms1.tst");
gap> f := GF(3);;
gap> gram := [[0,0,0,0,0,2],[0,0,0,0,2,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,2,0,0,0,0],[2,0,0,0,0,0]]*Z(3)^0;
[ [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3) ], 
  [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
  [ Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ] ]
gap> form := BilinearFormByMatrix(gram,f);
< bilinear form >
gap> Display(BaseChangeToCanonical(form));
 1 . 1 1 . 1
 1 . 2 2 . 1
 1 2 2 1 . 2
 2 . 1 2 2 1
 2 . 1 2 . 1
 . 1 2 1 2 .
gap> Display(form);
Hyperbolic bilinear form
Gram Matrix:
 . . . . . 2
 . . . . 2 .
 . . . 1 . .
 . . 1 . . .
 . 2 . . . .
 2 . . . . .
Witt Index: 3
gap> TypeOfForm(form);
1
gap> STOP_TEST("test_forms1.tst", 10000 );
