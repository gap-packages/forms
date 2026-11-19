gap> START_TEST("Forms: matobj/isisotropicvector.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> ToVecObj := {v, F} -> Vector(F, v);;
gap> mat := ToMatObj([[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]]*Z(41)^0, GF(41));
<4x4 plist matrix over GF(41)>
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> v := ToVecObj([1,1,0,0]*Z(41)^0, GF(41));
[ Z(41)^0, Z(41)^0, 0*Z(41), 0*Z(41) ]
gap> IsIsotropicVector(form,v);
true
gap> mat := ToMatObj([[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0, GF(8));
<5x5 plist matrix over GF(2^3)>
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> v1 := ToVecObj([1,0,0,0,0]*Z(8)^0, GF(8));
[ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
gap> v2 := ToVecObj([0,1,0,0,0]*Z(8)^0, GF(8));
[ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ]
gap> IsIsotropicVector(form,v1);
true
gap> IsIsotropicVector(form,v2);
true
gap> STOP_TEST("matobj/isisotropicvector.tst", 10000 );
