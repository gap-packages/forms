RandomVector:=function(F, n)
    local randvec, i;
    randvec := EmptyPlist(n);
    for i in [1..n] do
            randvec[i] := PseudoRandom(F);
    od;
    return randvec;
end;

RandomMatrix := function(n, F)
     local M, i, j;
     while true do
        # M := ZeroMatrix(F, n, n);
          M := NullMat(n, n, F);
          for i in [1..n] do
               for j in [1..n] do
                    M[i][j] := PseudoRandom(F);
               od;
          od;
          ConvertToMatrixRep(M, F);
          return M;
     od;
end;

TestPolyEval := function(benchmark)
    local iters, n, F, mat, coeffs, f, frob, eval, i, time_start, time_average_frob, time_average_normal, normal_eval;
    iters := 50;
    time_average_frob := 0;
    time_average_normal := 0;
    for i in [1..iters] do
        n := PseudoRandom([1..200]);
        F := GF(PseudoRandom([2, 3, 5])^(PseudoRandom([1, 2, 3])));
        mat := RandomMatrix(n, F);
        coeffs := RandomVector(F, PseudoRandom([0..300]));
        f := UnivariatePolynomial(F, coeffs);
        time_start := NanosecondsSinceEpoch();
        frob := FrobeniusNormalForm(mat);
        eval := __FORMSPACE__INTERNAL__EvaluatePolynomialWithFrobenius(coeffs, Matrix(mat), frob, Inverse(frob[2]), F, n);
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