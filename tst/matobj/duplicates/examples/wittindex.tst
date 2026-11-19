gap> START_TEST("Forms: matobj/wittindex.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[0,0,1,0,0],[0,0,0,0,0],[-1,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]*Z(7)^0, GF(7));
<5x5 plist matrix over GF(7)>
gap> form := BilinearFormByMatrix(mat,GF(7));
< bilinear form >
gap> WittIndex(form);
1
gap> RadicalOfForm(form);
<vector space of dimension 3 over GF(7)>
gap> Dimension(last);
3
gap> mat := ToMatObj(IdentityMat(6,GF(5)), GF(5));
<6x6 plist matrix over GF(5)>
gap> form := QuadraticFormByMatrix(mat,GF(5));
< quadratic form >
gap> WittIndex(form);
3
gap> mat := ToMatObj(IdentityMat(6,GF(7)), GF(7));
<6x6 plist matrix over GF(7)>
gap> form := QuadraticFormByMatrix(mat,GF(7));
< quadratic form >
gap> WittIndex(form);
2
gap> STOP_TEST("matobj/wittindex.tst", 10000 );
