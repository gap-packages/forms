#@local filt, f, dim, ToMatObj, v, w, frob, M, bil, her, quad, triv, id, her2,  MM, bil2, quad2, triv2, vv, ww
gap> START_TEST("Forms: /matobj/evaluate_form_respects_matrix_obj.tst"); # TODO i am not sure entirely this operation should be allowed...
gap> ToMatObj := m -> Matrix(IsGenericMatrixRep, GF(9), m);;
gap> MM := [[Z(9), Z(9)^2, Z(9)^3], [One(GF(9)), Z(9)^5, Zero(GF(9))]];;
gap> M := ToMatObj(MM);;
gap> v := RowsOfMatrix(M)[1];; w := RowsOfMatrix(M)[2];;
gap> vv := Unpack(v);; ww := Unpack(w);;
gap> frob := FrobeniusAutomorphism(GF(9));; id := IdentityMapping(GF(9));;
gap> bil := BilinearFormByMatrix(ToMatObj(IdentityMat(3, GF(9))), GF(9));;
gap> bil2 := BilinearFormByMatrix(IdentityMat(3, GF(9)), GF(9));;
gap> her := HermitianFormByMatrix(ToMatObj(IdentityMat(3, GF(9))), GF(9));;
gap> her2 := HermitianFormByMatrix(IdentityMat(3, GF(9)), GF(9));;
gap> quad := QuadraticFormByMatrix(ToMatObj(IdentityMat(3, GF(9))), GF(9));;
gap> quad2 := QuadraticFormByMatrix(IdentityMat(3, GF(9)), GF(9));;
gap> triv := BilinearFormByMatrix(ToMatObj(NullMat(3, 3, GF(9))), GF(9));;
gap> triv2 := BilinearFormByMatrix(NullMat(3, 3, GF(9)), GF(9));;
gap> IsTrivialForm(triv);
true
gap> v^frob;
<plist vector over GF(3^2) of length 3>
gap> ConstructingFilter(v^frob) = ConstructingFilter(v);
true
gap> v^id;
<immutable plist vector over GF(3^2) of length 3>
gap> v^id = v;
true
gap> M^frob = ToMatObj(MM^frob);
true
gap> ConstructingFilter(M^frob) = ConstructingFilter(M);
true
gap> M^id = M;
true
gap> [v, w]^bil;
Z(3^2)^6
gap> [M, M]^bil = ToMatObj([MM, MM]^bil2); 
true
gap> [v, w]^her;
Z(3^2)^5
gap> [M, M]^her = ToMatObj([MM, MM]^her2);
true
gap> [v, w]^triv;
0*Z(3)
gap> v^quad;
Z(3)
gap> M^quad = ToMatObj(MM^quad2);
true
gap> v^triv;
0*Z(3)
gap> ConstructingFilter([M, M]^bil) = ConstructingFilter(M);
true
gap> ConstructingFilter([M, M]^her) = ConstructingFilter(M);
true
gap> ConstructingFilter(M^quad) = ConstructingFilter(M);
true
gap> ConstructingFilter(M^id) = ConstructingFilter(M);
true
gap> ConstructingFilter(v^id) = ConstructingFilter(v);
true
gap> [v, w]^bil = [vv, ww]^bil2;
true
gap> [v, w]^her = [vv, ww]^her2;
true
gap> v^quad = vv^quad2;
true
gap> ([M, M]^bil)[1, 2] = [v, w]^bil;
true
gap> ([M, M]^her)[1, 2] = [v, w]^her;
true
gap> (M^quad)[1, 1] = v^quad;
true
gap> (v^frob)^frob = v and v^(frob^2) = v;
true
gap> [v, w]^her = ([w, v]^her)^frob;
true
gap> STOP_TEST("Forms: /matobj/evaluate_form_respects_matrix_obj.tst");
