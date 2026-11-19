gap> START_TEST("Forms: matobj/issingularvector.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> ToVecObj := {v, F} -> Vector(F, v);;
gap> mat := ToMatObj([[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0, GF(8));
<5x5 plist matrix over GF(2^3)>
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> v1 := ToVecObj([1,0,0,0,0]*Z(8)^0, GF(8));
[ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
gap> v2 := ToVecObj([0,1,0,0,0]*Z(8)^0, GF(8));
[ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ]
gap> IsSingularVector(form,v1);
false
gap> IsSingularVector(form,v2);
true
gap> IsIsotropicVector(form,v1);
true
gap> IsIsotropicVector(form,v2);
true
gap> STOP_TEST("matobj/issingularvector.tst", 10000 );
