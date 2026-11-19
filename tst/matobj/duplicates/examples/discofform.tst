gap> START_TEST("Forms: matobj/discofform.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> gram := ToMatObj(InvariantQuadraticForm(GO(-1,4,5))!.matrix, GF(5));;
gap> qform := QuadraticFormByMatrix(gram, GF(5));
< quadratic form >
gap> DiscriminantOfForm( qform );
"nonsquare"
gap> STOP_TEST("matobj/discofform.tst", 10000 );
