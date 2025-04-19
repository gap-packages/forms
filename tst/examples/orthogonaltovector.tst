gap> START_TEST("Forms: orthogonaltovector.tst");
gap> mat := [[0,0,0,-2],[0,0,-3,0],[0,3,0,0],[2,0,0,0]]*Z(7)^0;
[ [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^5 ], [ 0*Z(7), 0*Z(7), Z(7)^4, 0*Z(7) ], 
  [ 0*Z(7), Z(7), 0*Z(7), 0*Z(7) ], [ Z(7)^2, 0*Z(7), 0*Z(7), 0*Z(7) ] ]
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> v := [0*Z(7),Z(7)^0,Z(7)^3,Z(7)^5];
[ 0*Z(7), Z(7)^0, Z(7)^3, Z(7)^5 ]
gap> vperp := OrthogonalSubspaceMat(form,v);
[ [ Z(7)^0, Z(7)^0, 0*Z(7), 0*Z(7) ], [ Z(7)^0, 0*Z(7), Z(7)^0, 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ] ]
gap> List(vperp,x->[x,v]^form);
[ 0*Z(7), 0*Z(7), 0*Z(7) ]
gap> sub := [[1,1,0,0],[0,0,1,2]]*Z(7)^0;
[ [ Z(7)^0, Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), Z(7)^0, Z(7)^2 ] ]
gap> subperp := OrthogonalSubspaceMat(form,sub);
[ [ Z(7)^0, Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), Z(7)^4, Z(7)^0 ] ]
gap> List(subperp,x->List(sub,y->[x,y]^form));
[ [ 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7) ] ]
gap> mat := [[1,0,0],[0,0,1],[0,0,0]]*Z(2)^0;
[ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> v := [Z(2)^0,Z(2)^0,0*Z(2)];
[ Z(2)^0, Z(2)^0, 0*Z(2) ]
gap> vperp := OrthogonalSubspaceMat(form,v);
[ <an immutable GF2 vector of length 3>, <an immutable GF2 vector of length 
    3> ]
gap> bil_form := AssociatedBilinearForm(form);
< bilinear form >
gap> List(vperp,x->[x,v]^bil_form);
[ 0*Z(2), 0*Z(2) ]
gap> sub := [[1,0,1],[1,0,0]]*Z(2)^0;
[ [ Z(2)^0, 0*Z(2), Z(2)^0 ], [ Z(2)^0, 0*Z(2), 0*Z(2) ] ]
gap> subperp := OrthogonalSubspaceMat(form,sub);
[ <an immutable GF2 vector of length 3>, <an immutable GF2 vector of length 
    3> ]
gap> List(subperp,x->List(sub,y->[x,y]^bil_form));
[ [ 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2) ] ]
gap> STOP_TEST("orthogonaltovector.tst", 10000 );
