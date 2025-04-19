gap> START_TEST("Forms: discofform.tst");
gap> gram := InvariantQuadraticForm(GO(-1,4,5))!.matrix;;
gap> qform := QuadraticFormByMatrix(gram, GF(5));
< quadratic form >
gap> DiscriminantOfForm( qform );
"nonsquare"
gap> STOP_TEST("discofform.tst", 10000 );
