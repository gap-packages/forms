gap> START_TEST("Forms: matobj/orthogonaltovector.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> ToVecObj := {v, F} -> Vector(F, v);;
gap> mat := ToMatObj([[0,0,0,-2],[0,0,-3,0],[0,3,0,0],[2,0,0,0]]*Z(7)^0, GF(7));
<4x4 plist matrix over GF(7)>
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> v := ToVecObj([0*Z(7),Z(7)^0,Z(7)^3,Z(7)^5], GF(7));;
gap> Unpack(v);
[ 0*Z(7), Z(7)^0, Z(7)^3, Z(7)^5 ]
gap> vperp := OrthogonalSubspaceMat(form,v);;
gap> Unpack(vperp);
[ [ Z(7)^0, Z(7)^0, 0*Z(7), 0*Z(7) ], [ Z(7)^0, 0*Z(7), Z(7)^0, 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ] ]
gap> List(vperp,x->[x,v]^form);
[ 0*Z(7), 0*Z(7), 0*Z(7) ]
gap> sub := ToMatObj([[1,1,0,0],[0,0,1,2]]*Z(7)^0, GF(7));;
gap> Unpack(sub);
[ [ Z(7)^0, Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), Z(7)^0, Z(7)^2 ] ]
gap> subperp := OrthogonalSubspaceMat(form,sub);
[ [ Z(7)^0, Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), Z(7)^4, Z(7)^0 ] ]
gap> List(subperp,x->List(sub,y->[x,y]^form));
[ [ 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7) ] ]
gap> mat := ToMatObj([[1,0,0],[0,0,1],[0,0,0]]*Z(2)^0, GF(2));
<3x3 plist matrix over GF(2)>
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> v := ToVecObj([Z(2)^0,Z(2)^0,0*Z(2)], GF(2));;
gap> Unpack(v);
[ Z(2)^0, Z(2)^0, 0*Z(2) ]
gap> vperp := OrthogonalSubspaceMat(form,v);
<an immutable 2x3 matrix over GF2>
gap> bil_form := AssociatedBilinearForm(form);
< bilinear form >
gap> List(vperp,x->[x,v]^bil_form);
[ 0*Z(2), 0*Z(2) ]
gap> sub := ToMatObj([[1,0,1],[1,0,0]]*Z(2)^0, GF(2));;
gap> Unpack(sub);
[ [ Z(2)^0, 0*Z(2), Z(2)^0 ], [ Z(2)^0, 0*Z(2), 0*Z(2) ] ]
gap> subperp := OrthogonalSubspaceMat(form,sub);
<an immutable 2x3 matrix over GF2>
gap> List(subperp,x->List(sub,y->[x,y]^bil_form));
[ [ 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2) ] ]
gap> STOP_TEST("matobj/orthogonaltovector.tst", 10000 );
