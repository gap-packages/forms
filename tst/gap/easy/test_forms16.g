#EvaluateForm and Isotropic subspaces: hermitian form
q := 5;
f := GF(q^4);
mat := NullMat(3,3,f);
mat[1][2] := Z(q^4);
mat[2][1] := Z(q^4)^(q^2);
mat[3][3] := Z(q)^0;
form := HermitianFormByMatrix(mat,GF(q^4));
v := f^3;
lines := Subspaces(v,1);
matrices := List(lines,x->BasisVectors(Basis(x)));;
vectors := List(matrices,x->x[1]);;
results := Collected(List(vectors,x->EvaluateForm(form,x,x)));;
[Zero(f),q^6+1] in results;
results := Collected(List(matrices,x->EvaluateForm(form,x,x)));;
[[[Zero(f)]],q^6+1] in results;
Number(vectors,x->IsIsotropicVector(form,x))=q^6+1;
Number(matrices,x->IsTotallyIsotropicSubspace(form,x))=q^6+1;
quit;
