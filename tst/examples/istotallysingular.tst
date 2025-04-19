gap> START_TEST("Forms: istotallysingular.tst");
gap> mat := [[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0;
[ [ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> sub := [[Z(2)^0,0*Z(2),Z(2^3)^6,Z(2^3),Z(2^3)^3],
>        [0*Z(2),Z(2)^0,Z(2^3)^6,Z(2^3)^2,Z(2^3)]];
[ [ Z(2)^0, 0*Z(2), Z(2^3)^6, Z(2^3), Z(2^3)^3 ], 
  [ 0*Z(2), Z(2)^0, Z(2^3)^6, Z(2^3)^2, Z(2^3) ] ]
gap> IsTotallySingularSubspace(form,sub);
true
gap> STOP_TEST("istotallysingular.tst", 10000 );
