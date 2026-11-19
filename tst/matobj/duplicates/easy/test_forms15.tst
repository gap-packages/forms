gap> START_TEST("Forms: matobj/test_forms15.tst");
gap> q := 5;
5
gap> f := GF(q);
GF(5)
gap> dim := 5;
5
gap> v := f^dim;
( GF(5)^5 )
gap> mat := Matrix(IsPlistMatrixRep, f, IdentityMat(dim,f));
<5x5 plist matrix over GF(5)>
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> lines := Subspaces(v,1);
Subspaces( ( GF(5)^5 ), 1 )
gap> matrices := List(lines,x->Matrix(BasisVectors(Basis(x)), mat));;
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
gap> planes := Subspaces(v,2);
Subspaces( ( GF(5)^5 ), 2 )
gap> matrices := List(planes,x->Matrix(BasisVectors(Basis(x)), mat));;
gap> results := Collected(List(matrices,x->x^form));;
gap> [ZeroMatrix(2,2,mat),(q^(dim-1)-1)/(q-1)] in results;
true
gap> Number(matrices,x->IsTotallySingularSubspace(form,x))=(q^(dim-1)-1)/(q-1);
true
gap> STOP_TEST("matobj/test_forms15.tst", 10000 );
