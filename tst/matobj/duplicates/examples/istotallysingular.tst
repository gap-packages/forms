gap> START_TEST("Forms: matobj/istotallysingular.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0, GF(8));
<5x5 plist matrix over GF(2^3)>
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> sub := ToMatObj([[Z(2)^0,0*Z(2),Z(2^3)^6,Z(2^3),Z(2^3)^3],
>        [0*Z(2),Z(2)^0,Z(2^3)^6,Z(2^3)^2,Z(2^3)]], GF(8));
<2x5 plist matrix over GF(2^3)>
gap> IsTotallySingularSubspace(form,sub);
true
gap> STOP_TEST("matobj/istotallysingular.tst", 10000 );
