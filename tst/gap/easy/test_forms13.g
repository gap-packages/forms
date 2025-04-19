#isotropic subspaces. quadratic forms
q := 9;
f := GF(q);
dim := 3;
v := f^dim;
mat := IdentityMat(dim,f);
form := QuadraticFormByMatrix(mat,f);
lines := Subspaces(v,1);
matrices := List(lines,x->BasisVectors(Basis(x)));;
vectors := List(matrices,x->x[1]);;
results := Collected(List(vectors,x->EvaluateForm(form,x)));;
[Zero(f),(q^(dim-1)-1)/(q-1)] in results;
results := Collected(List(matrices,x->x^form));;
[[[Zero(f)]],(q^(dim-1)-1)/(q-1)] in results;
Number(vectors,x->IsSingularVector(form,x))=(q^(dim-1)-1)/(q-1);
Number(matrices,x->IsTotallySingularSubspace(form,x))=(q^(dim-1)-1)/(q-1);
quit;
