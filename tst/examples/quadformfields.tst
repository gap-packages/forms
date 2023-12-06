gap> START_TEST("Forms: quadformfields.tst");
gap> #constructing the same form over different fields
gap> mat := 
> [[Z(2)^0,Z(2)^0,0*Z(2),0*Z(2)],[0*Z(2),Z(2)^0,0*Z(2),0*Z(2)], 
>  [0*Z(2),0*Z(2),0*Z(2),Z(2)^0],[0*Z(2),0*Z(2),0*Z(2),0*Z(2)]];
[ [ Z(2)^0, Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> WittIndex(form);
1
gap> form := QuadraticFormByMatrix(mat,GF(4));
< quadratic form >
gap> WittIndex(form);
2
gap> STOP_TEST("quadformfields.tst", 10000 );
