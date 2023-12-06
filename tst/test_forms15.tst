gap> START_TEST("Forms: test_forms15.tst");
gap> q := 5;
5
gap> f := GF(q);
GF(5)
gap> dim := 5;
5
gap> v := f^dim;
( GF(5)^5 )
gap> mat := IdentityMat(dim,f);
[ [ Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0 ] ]
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> lines := Subspaces(v,1);
Subspaces( ( GF(5)^5 ), 1 )
gap> matrices := List(lines,x->BasisVectors(Basis(x)));;
gap> vectors := List(matrices,x->x[1]);;
gap> results := Collected(List(vectors,x->EvaluateForm(form,x)));;
gap> [Zero(f),(q^(dim-1)-1)/(q-1)] in results;
true
gap> results := Collected(List(matrices,x->x^form));;
gap> [[[Zero(f)]],(q^(dim-1)-1)/(q-1)] in results;
true
gap> Number(vectors,x->IsSingularVector(form,x))=(q^(dim-1)-1)/(q-1);
true
gap> Number(matrices,x->IsTotallySingularSubspace(form,x))=(q^(dim-1)-1)/(q-1);
true
gap> planes := Subspaces(v,2);
Subspaces( ( GF(5)^5 ), 2 )
gap> matrices := List(planes,x->BasisVectors(Basis(x)));;
gap> results := Collected(List(matrices,x->x^form));;
gap> z := Zero(f);
0*Z(5)
gap> [[[z,z],[z,z]],(q^(dim-1)-1)/(q-1)] in results;
true
gap> Number(matrices,x->IsTotallySingularSubspace(form,x))=(q^(dim-1)-1)/(q-1);
true
gap> STOP_TEST("test_forms15.tst", 10000 );
