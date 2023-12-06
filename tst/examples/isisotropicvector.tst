gap> START_TEST("Forms: isisotropicvector.tst");
gap> #testing isotropy for vectors.
gap> mat := [[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]]*Z(41)^0;
[ [ Z(41)^0, 0*Z(41), 0*Z(41), 0*Z(41) ], 
  [ 0*Z(41), Z(41)^20, 0*Z(41), 0*Z(41) ], 
  [ 0*Z(41), 0*Z(41), 0*Z(41), Z(41)^0 ], 
  [ 0*Z(41), 0*Z(41), Z(41)^0, 0*Z(41) ] ]
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> v := [1,1,0,0]*Z(41)^0;
[ Z(41)^0, Z(41)^0, 0*Z(41), 0*Z(41) ]
gap> IsIsotropicVector(form,v);
true
gap> mat := [[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0;
[ [ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> v1 := [1,0,0,0,0]*Z(8)^0;
[ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
gap> v2 := [0,1,0,0,0]*Z(8)^0;
[ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ]
gap> IsIsotropicVector(form,v1);
true
gap> IsIsotropicVector(form,v2);
true
gap> STOP_TEST("isisotropicvector.tst", 10000 );
