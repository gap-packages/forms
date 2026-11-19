gap> START_TEST("Forms: matobj/quadformfields.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := 
> ToMatObj([[Z(2)^0,Z(2)^0,0*Z(2),0*Z(2)],[0*Z(2),Z(2)^0,0*Z(2),0*Z(2)], 
>  [0*Z(2),0*Z(2),0*Z(2),Z(2)^0],[0*Z(2),0*Z(2),0*Z(2),0*Z(2)]], GF(2));
<4x4 plist matrix over GF(2)>
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> WittIndex(form);
1
gap> form := QuadraticFormByMatrix(ChangedBaseDomain(mat, GF(4)),GF(4));
< quadratic form >
gap> WittIndex(form);
2
gap> STOP_TEST("matobj/quadformfields.tst", 10000 );
