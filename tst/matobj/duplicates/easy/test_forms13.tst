gap> START_TEST("Forms: matobj/test_forms13.tst");
gap> q := 9;
9
gap> f := GF(q);
GF(3^2)
gap> dim := 3;
3
gap> v := f^dim;
( GF(3^2)^3 )
gap> mat := Matrix(IsPlistMatrixRep, f, IdentityMat(dim,f));
<3x3 plist matrix over GF(3^2)>
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> lines := Subspaces(v,1);
Subspaces( ( GF(3^2)^3 ), 1 )
gap> matrices := List(lines,x->Matrix(BasisVectors(Basis(x)), mat));; # have to convert to matrix objects as list times matrix obj multiplication is dubious at this point in time. See https://github.com/gap-system/gap/issues/6635 This is in a lot of the test files, namely test_forms13, test_forms14, test_forms15, test_forms16
gap> vectors := List(matrices,x->RowsOfMatrix(x)[1]);;
gap> results := Collected(List(vectors,x->EvaluateForm(form,x)));;
gap> [Zero(f),(q^(dim-1)-1)/(q-1)] in results;
true
gap> results := Collected(List(matrices,x->x^form));;
gap> [ZeroMatrix(1,1,mat),(q^(dim-1)-1)/(q-1)] in results;
true
gap> Number(vectors,x->IsSingularVector(form,x))=(q^(dim-1)-1)/(q-1);
true
gap> Number(matrices,x->IsTotallySingularSubspace(form,x))=(q^(dim-1)-1)/(q-1);
true
gap> STOP_TEST("matobj/test_forms13.tst", 10000 );
