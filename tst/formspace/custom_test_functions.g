TestPolyEval := function(benchmark)
    local iters, n, F, mat, coeffs, f, frob, eval, i, time_start, time_average_frob, time_average_normal, normal_eval, x;
    iters := 50;
    time_average_frob := 0;
    time_average_normal := 0;
    for i in [1..iters] do
        n := Random([1..200]);
        F := GF(Random([2, 3, 5])^(Random(1, 3)));
        mat := RandomMat(n, n, F);
        coeffs := List([1..Random(1, 300)], x -> Random(F));
        f := UnivariatePolynomial(F, coeffs);
        time_start := NanosecondsSinceEpoch();
        frob := FrobeniusNormalForm(mat);
        eval := FORMS_EvaluatePolynomialWithFrobenius(coeffs, Matrix(mat), frob, Inverse(frob[2]), F, n);
        time_average_frob := time_average_frob + (NanosecondsSinceEpoch() - time_start);

        time_start := NanosecondsSinceEpoch();
        normal_eval := f(mat);
        time_average_normal := time_average_normal + (NanosecondsSinceEpoch() - time_start);
        if eval <> normal_eval then
            Error("Polynomial: ", f, " Matrix: ", mat, " Frobenius Normal Form: ", frob, "\n");
        fi;
    od;
    time_average_frob := Float(time_average_frob) / Float((iters * 1000000)); # convert to ms
    time_average_normal := Float(time_average_normal) / Float((iters  * 1000000));
    if benchmark then
        Print("Frob: ", time_average_frob, " Normal: ", time_average_normal, "\n");
    fi;
    return "Ok";
end;

TestMatricesAreForms := function(G, Lambdas, unitary, forms)
    local p_exponent, hom, F, f, g, i, Gens, n, j;
    F := DefaultFieldOfMatrixGroup(G);
    Gens := GeneratorsOfGroup(G);
    n := DimensionsMat(Gens[1])[1];
    p_exponent := DegreeOverPrimeField(F);
    hom := fail;
    if unitary then
        if p_exponent mod 2 <> 0 then
            if Size(forms) = 0 then
                return "Ok";
            else 
                Error("Claims to have found unitary form although the field does not admit a field automorphism of order two!");
            fi;
        fi;
        hom := FrobeniusAutomorphism(F)^(p_exponent/2);
    fi;

    for i in [1..Size(forms)] do
        f := forms[i];
        for j in [1..Size(Gens)] do
            g := Gens[j];
            if g * f * FORMS_CalculateAdjoint(g, unitary, hom, n, F) <> Lambdas[j] * f then
                Error("Computed non formspace element ", f, "group " , G, " unitary ", unitary, " Lambdas ", Lambdas);
            fi;
        od;
    od;
    return "Ok";
end;

TestMatricesAreForms2 := function(G, forms)
    local F, n, lambdas, Gens, i;
    F := DefaultFieldOfMatrixGroup(G);
    Gens := GeneratorsOfGroup(G);
    n := DimensionsMat(Gens[1])[1];
    lambdas := [];
    for i in [1..Size(Gens)] do
        Add(lambdas, One(F));
    od;
    return [TestMatricesAreForms(G, lambdas, false, forms[1]), TestMatricesAreForms(G, lambdas, true, forms[2])];
end;

# computes F by solving big system of linear equations.
# this is a good idea for testing to see that the presered formspace function actually computes the entire formspace not just a subspace. (this is the function that provides the formsace solutions used in other tests)
# however maybe we should test the test? this seems a little stupid
TestComputeFormspaceBruteForce := function(G, Lambdas, unitary)
    local Gens, i, j, F, n, base, b, g, eqs, p_exponent, ComputeMatrixVectorSpaceBasis, VectorToMatrix, MatrixToVector, sol, out, s, hom, v;
    Gens := GeneratorsOfGroup(G);
    n := DimensionsMat(Gens[1])[1];
    F := DefaultFieldOfMatrixGroup(G);

    ComputeMatrixVectorSpaceBasis := function(F, n)
        local O, i, j, m;
        O := [];
        for i in [1..n] do
            for j in [1..n] do
                m := NullMat(n, n, F);
                m[i][j] := One(F);
                Add(O, m);
            od;
        od;
        return O;
    end;

    VectorToMatrix := function(vec, F, n)
        local i, j, m;
        m := NullMat(n, n, F);
        for i in [1..n] do
            for j in [1..n] do
                m[i][j] := vec[(i-1)*n + j];
            od;
        od;
        return m;
    end;

    MatrixToVector := function(mat, F, n)
        local vec, i, j;
        vec := ZeroVector(F, n^2);
        for i in [1..n] do
            for j in [1..n] do
                vec[(i-1)*n + j] := mat[i][j];
            od;
        od;
        return vec;
    end;

    p_exponent := DegreeOverPrimeField(F);
    hom := fail;
    if unitary then
        if p_exponent mod 2 <> 0 then
            return [];
        fi;
        hom := FrobeniusAutomorphism(F)^(p_exponent/2);
    fi;
    base := ComputeMatrixVectorSpaceBasis(F, n);
    eqs := [];
   
    for j in [1..Size(base)] do  
        for i in [1..Size(Gens)] do
            b := base[j];  
            v :=  MatrixToVector(Gens[i] * b * FORMS_CalculateAdjoint(Gens[i], unitary, hom, n, F) - Lambdas[i] * b, F, n);
            if i = 1 then
                Add(eqs, v);
            else
                eqs[j] := Concatenation(List(eqs[j]), List(v));
            fi;
        od;
    od;
    # Print(aaa);
    sol := NullspaceMat(eqs);
    out := [];
    for s in sol do
        Add(out, VectorToMatrix(s, F, n));
    od;
    # Print(aaa);
    return out;
end;
