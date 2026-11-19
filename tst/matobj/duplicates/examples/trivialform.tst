gap> START_TEST("Forms: matobj/trivialform.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[0,0,0],[0,0,0],[0,0,0]]*Z(7)^0, GF(7));
<3x3 plist matrix over GF(7)>
gap> form1 := BilinearFormByMatrix(mat,GF(7));
< trivial form >
gap> form2 := QuadraticFormByMatrix(mat,GF(7));
< trivial form >
gap> form1 = form2;
true
gap> IsQuadraticForm(form1);
false
gap> IsSesquilinearForm(form1);
false
gap> mat := ToMatObj([[0,0],[0,0]]*Z(4)^0, GF(4));
<2x2 plist matrix over GF(2^2)>
gap> form3 := BilinearFormByMatrix(mat,GF(4));
< trivial form >
gap> form3 = form1;
false
gap> STOP_TEST("matobj/trivialform.tst", 10000 );
