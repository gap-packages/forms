gap> START_TEST("Forms: bg_th_ex1.tst");
gap> #Background theory: example 1
gap> mat := [[1,0,0],[0,1,4],[1,2,1]]*Z(5)^0;
[ [ Z(5)^0, 0*Z(5), 0*Z(5) ], [ 0*Z(5), Z(5)^0, Z(5)^2 ], 
  [ Z(5)^0, Z(5), Z(5)^0 ] ]
gap> form := BilinearFormByMatrix(mat,GF(5));
Error, Invalid Gram matrix
gap> STOP_TEST("bg_th_ex1.tst", 10000 );
