gap> START_TEST("Forms: test_forms2.tst");
gap> f := GF(7);
GF(7)
gap> gram := [[-3,0,0,0,0,0],[0,0,3,0,0,0],[0,3,0,0,0,0],[0,0,0,0,0,-1/2],[0,0,0,0,1,0],[0,0,0,-1/2,0,0]]*Z(7)^0;
[ [ Z(7)^4, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7), 0*Z(7), 0*Z(7) ] ]
gap> form := BilinearFormByMatrix(gram,f);
< bilinear form >
gap> IsEllipticForm(form);
true
gap> TypeOfForm(form);
-1
gap> Display(form);
Elliptic bilinear form
Gram Matrix:
 4 . . . . .
 . . 3 . . .
 . 3 . . . .
 . . . . . 3
 . . . . 1 .
 . . . 3 . .
Witt Index: 2
gap> STOP_TEST("test_forms2.tst", 10000 );
