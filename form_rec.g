#Print("Dont forget to load M. Gecks NoFoMa Package!\n");

# IN
# g : Zyklische Matrix
# v : Vektor der zyklisch ist
# OUT
# Basis (v gv g^2v ... g^(n-1)v)
Spin := function(g, v) 
    local B, n, i;
    n := DimensionsMat(g)[1];
    B := NullMat(n, n, Zero(Field(v[1])));
    B[1] := v;

    for i in [2..n] do
        v := v*g;
        B[i] := v;
    od;
    # return CMat(B);
    ConvertToMatrixRep(B, DefaultFieldOfMatrix(B));
    return B;
end;

RandomVector:=function(F, n)
    local randvec, i;
    randvec := EmptyPlist(n);
    for i in [1..n] do
            randvec[i] := PseudoRandom(F);
    od;
    return randvec;
end;
## TODO use predefined func!!
ZeroVectorFunc:=function(F, n)
    local i, v;
    v := EmptyPlist(n);
    for i in [1..n] do
        v[i] := Zero(F);
    od;
    return v;
end;

IsZeroVec := function(v)
    local i;
    for i in [1..Length(v)] do
        if not IsZero(v[i]) then
            return false;
        fi;
    od;
    return true;
end;

EvaluateMatrixPolynomialWithVector := function(F, n, g, v, coeffs)
    local res, i, deg;
    deg := Size(coeffs);
    if deg = 0 then
        return ZeroVector(F, n);
    fi;
    if deg = 1 then
        return v * coeffs[1];
    fi;
    
    res := v * g * coeffs[deg];
    for i in [1..deg-2] do
        res := res + coeffs[deg - i]*v;
        res := res * g;
    od;
    res := res + coeffs[1]*v;
    return res;
end;

# mode = False \iff mat* := mat^tr, mode = True \iff mat* = komplex_konjugiert(mat^tr)
CalculateAdjoint := function(mat, mode, hom, n)
    local transposed, i, j;
    transposed := TransposedMat(mat);
    if mode = false then
        return transposed;
    fi;
    if mode then
        return transposed^hom;
    fi;
    return fail;
end;

# follows forms package recognition.gi (See : https://github.com/gap-packages/forms/blob/master/lib/recognition_new.gi)
# the reason for own implementation is to reuse the minpol computation (maybe just ignore this and use the one from forms package anyways???)
CalculateLambdas := function(F, n, g, g_star_inv, mode, hom)
    local t, t_star_inv, p, p_deg, cyc, I, as, gcd_rep, l, a;
    t := Trace(g);
    t_star_inv := Trace(g_star_inv);
    cyc := fail;
    if t = Zero(F) and t_star_inv <> Zero(F) then
        return fail;
    fi;
    if t <> Zero(F) and t_star_inv = Zero(F) then
        return fail;
    fi;
    if t <> Zero(F) and t_star_inv <> Zero(F) then 
        # TODO
        return [t * Inverse(t_star_inv), fail];
    fi;
    # diser ganze ansatz hat den nachteil, das die suche nach dem minimalpolynom eigentlich einen zyklischen vektor liefert. Am besten geben wir diesen zurück statt true/false/fail. Dafür müsste aber die methode für das minpol angepasst werden.
    p := MinimalPolynomial(g);
    
    # the reasom
    p_deg := Degree(p);
    as := CoefficientsOfUnivariatePolynomial(p);
    I := Filtered([0..p_deg], x -> as[x+1] <> Zero(F));
    
    # muss symmetrisch sein, sonst wird keine non-deg form invariant gelassen.
    if ForAny(I, x -> not (p_deg-x) in I) then
        return fail;
    fi;

    gcd_rep := GcdRepresentation(I);
    
    if mode = false then
        l:=List([1..Length(I)-1], x ->((as[1])*as[p_deg-I[x]+1]/(as[I[x]+1])));
    else
        l:=List([1..Length(I)-1], x ->((as[1]^hom)*as[p_deg-I[x]+1]/(as[I[x]+1]^hom)));
    fi;
    # todo a is not guranateed to be Lambdas here, only Lambdas^gcd = a!
    a:= Product([1..Length(I)-1], x->l[x]^gcd_rep[x]);

    if Degree(p) = n then
        cyc := true;
    else 
        cyc := false;
    fi;
    #TODO komplizierte analyse um Lambdas zu finden..
    return [a, cyc];
end;

RandomMatrixInvertible := function(n, q)
    local F, M, i, j;
    F := GF(q);
    while true do
        # M := ZeroMatrix(F, n, n);
        M := NullMat(n, n, F);
        for i in [1..n] do
            for j in [1..n] do
                M[i][j] := PseudoRandom(F);
            od;
        od;
        ConvertToMatrixRep(M, F);
        if RankMat(M) = n then
            return M;
        fi;
    od;
end;

##### --------- uses functions from Geck "nofoma" package

# Avoids the use of PseudoRandomElement(G), since it is slow and the actual "randomness" is irrelevant.
#returns g cyclic, scalar that belongs to g, 
# has some issues!!!!!!! TODO FIX ME
FindCyclicGroupElementAndScalars := function(Gens, Lambdas)
    local cur_group_element, cur_scalar, known_elements, i, mode, n, res, j, known_scalars, best_known_element_index, best_known_res, best_known_length, mod_elem, g, e;

    known_elements := ShallowCopy(Gens);
    known_scalars := ShallowCopy(Lambdas);

    n := DimensionsMat(Gens[1])[1];
    i := Size(known_elements);
    mod_elem := 1;
    best_known_length := n + 1;
    while i < 25 do
        # no accidental identity mat
        if false then
            mode := 1;
        else
            mode := PseudoRandom([1, 2]);
        fi;
        j := PseudoRandom([1..Size(known_elements)]);
        cur_group_element := known_elements[j];
        cur_scalar := known_scalars[j];

        if mode = 1 then
            j := PseudoRandom([1..Size(known_elements)]);
            cur_group_element := cur_group_element * known_elements[j];
            cur_scalar := cur_scalar * known_scalars[j];
        elif mode = 2 then
            j := PseudoRandom([1..Size(known_elements)]);
            cur_group_element := cur_group_element * Inverse(known_elements[j]);
            cur_scalar := cur_scalar * Inverse(known_scalars[j]);
        fi;
        if i mod mod_elem = 0 then
            res := FrobeniusNormalForm(cur_group_element);
            if Size(res[3]) = 1 then
                return Concatenation([cur_group_element, cur_scalar], res, [i]); 
            fi;
            if Size(res[3]) <= best_known_length and not (cur_group_element in Gens) then
                best_known_element_index := i + 1;
                best_known_res := res;
                best_known_length := Size(res[3]);
            fi;
        fi;
        Add(known_elements, cur_group_element);
        Add(known_scalars, cur_scalar);
        i := i + 1;
    od;
    ## This fixes all the issues for some ungodly reason!!! why :(
    # best_known_element_index := PseudoRandom(Group(Gens));
    # return Concatenation([best_known_element_index, 1], FrobeniusNormalForm(best_known_element_index), [-1]);
    Print("only found element of length ", best_known_length, "\n");
    # TODO: möglichst kurze zyklische modul basis usw....
    return Concatenation([known_elements[best_known_element_index], known_scalars[best_known_element_index]], best_known_res, [i]);
    # return fail;
end;

# turns the jn vector into a j times n matrix
VectorReorganize := function(vec, j, F, n)
    local A, i;
    A := NullMat(j, n, F);
    ConvertToMatrixRep(A, F);
    for i in [1..j] do
        A[i] := vec{[((i - 1) * n + 1)..(i*n)]};
    od;
    return A;
end;

# p is given as a list of coefficients.
EvaluatePolynomialWithFrobenius := function(p, g, frob_base, frob_base_inv, F, n)
    local ws, C, i, end_pos, j, k;
    ws := [];
    j := Size(frob_base[3]);
    for k in [1..j] do
    #  function(F, n, g, v, coeffs)
        Add(ws, EvaluateMatrixPolynomialWithVector(F, n, g, frob_base[2][frob_base[3][k]]{[1..n]}, p));
    od;
    C := NullMat(n, n, F);
    ConvertToMatrixRep(C, F);
    for k in [1..j] do
        end_pos := 0;
        if j = k then
            end_pos := n;
        else 
            end_pos := frob_base[3][k + 1];
        fi;
        for i in [0..(end_pos - frob_base[3][k])] do
            if i = 0 then
                C[frob_base[3][k] + i] := ws[k];
            else
                C[frob_base[3][k] + i] := C[frob_base[3][k] + i - 1]*g;
            fi;
        od;
    od;
    ## FrobBasSpin needed?
    return frob_base_inv * C;
end;

ComputeConditionMatrixFrob := function(u, h, h_star, scalar_h, g_star_inv_scaled, frob_base, frob_base_inv, frob_base_inv_star, F, n)
    local coeffs_c, coeffs_f, Ps, i, j, b_end, b_start, cpol, fpol;
    coeffs_c := (u * h) * frob_base_inv;
    coeffs_f := (u * frob_base_inv) * scalar_h;
    # Display(coeffs_c);
    j := Size(frob_base[3]);
    Ps := NullMat(n * j, n, F);
    ConvertToMatrixRep(Ps, F);
    # Print(frob_base[3], "\n");
    for i in [1..j] do
        if i = j then
            b_end := n;
        else
            b_end := frob_base[3][i + 1] - 1;
        fi;
        # Print("[", ((i - 1)*n + 1), ",", (i*n), "\n");
        #Print("->", frob_base[3][i],",",b_end, "\n");
        Ps{[((i - 1)*n + 1)..(i*n)]}{[1..n]} :=
            EvaluatePolynomialWithFrobenius(coeffs_c{[frob_base[3][i]..b_end]}, g_star_inv_scaled, frob_base, frob_base_inv_star, F, n) * h_star - EvaluatePolynomialWithFrobenius(coeffs_f{[frob_base[3][i]..b_end]}, g_star_inv_scaled, frob_base, frob_base_inv_star, F, n);
        # cpol := UnivariatePolynomial(F, coeffs_c);
        # fpol := UnivariatePolynomial(F, coeffs_f);
        # Ps{[((i - 1)*n + 1)..(i*n)]}{[1..n]} :=
        #     cpol(g_star_inv_scaled) * h_star - fpol(g_star_inv_scaled);
        # Print(aua);
    od;
    return Ps;
end;

FrobSpin := function(Images, spin_elem, frob_base_blocks, n, F)
    local A, j, i, k, end_pos;
    j := Size(frob_base_blocks);
    A := NullMat(n, n, F);
    ConvertToMatrixRep(A, F);
    for i in [1..j] do
        if i = j then
            end_pos := n;
        else
            end_pos := frob_base_blocks[i + 1] - 1;
        fi;
        A[frob_base_blocks[i]] := Images[i];
        for k in [(frob_base_blocks[i] + 1)..end_pos] do
            A[k] := A[k - 1]*spin_elem;
        od;
    od;
    return A;
end;

FindFormspaceInternal := function(Gens, Lambdas, unitary, hom, g_res, g_inv_frob, g_star_inv_scaled, g_star_inv_scaled_frob, frob_base_inv_star, d, F, n)
    local k, i, W, first, j, Cond, Conds, h_star, h, w, O, needs_checking, A, failed_check, nspace, vec;
    first := true;
    needs_checking := false; ## maybe remove?

    for i in [1..d] do
        if needs_checking then
            break;
        fi;
        h := Gens[i];
        h_star := CalculateAdjoint(h, unitary, hom, n);
        for j in [1..n] do
            # vec := RandomVector(F, n);
            vec := g_res[4][j];
            # Display(vec);
            Conds := 
                ComputeConditionMatrixFrob(vec, h, h_star, Lambdas[i], g_star_inv_scaled, g_star_inv_scaled_frob, g_inv_frob, frob_base_inv_star, F, n);
            
            ConvertToMatrixRep(Conds, F);
            if not first then
                nspace := NullspaceMat(W * Conds);
                ConvertToMatrixRep(nspace, F);
                if Size(nspace) = 0 then 
                    # Print("empty");
                    return [];
                fi;
                W := nspace * W;
            fi;
            if first then
                nspace := NullspaceMat(Conds);
                ConvertToMatrixRep(nspace, F);
                if Size(nspace) = 0 then 
                    # Print("empty");
                    return [];
                fi;
                W := nspace;
                first := false;
            fi;
            if Size(nspace) = 1 then
                needs_checking := true;
                break;
            fi;
        od;
    od;
    O := [];
    for w in W do
        A := g_inv_frob * FrobSpin(VectorReorganize(w, Size(g_res[5]), F, n), g_star_inv_scaled, g_res[5], n, F);
        if needs_checking then
            failed_check := false;
            for i in [1..d] do
                if not failed_check and Gens[i] * A * CalculateAdjoint(Gens[i], unitary, hom, n) <> Lambdas[i] * A then
                    failed_check := true;
                fi;
            od;
            if not failed_check then
                Add(O, A);
            fi;
        fi;
        if not needs_checking then
            Add(O, A);
        fi;
    od;

    return O;

end;

FindFormspace := function(G, Lambdas, unitary) 
    local F, p_exponent, Gens, n, d, first, hom, g_res, g_inv_frob, g_star_inv_unscaled, g_star_inv_scaled_frob, frob_base_inv_star;
    if not IsMatrixGroup(G) then
        # Print("The method only works for matrix groups. \n");
        return;
    fi;
    F := DefaultFieldOfMatrixGroup(G);
    p_exponent := DegreeOverPrimeField(F);
    if unitary and (p_exponent mod 2 <> 0) then
        Print("Field does not admit field automorphism of order two!");
        return;
    fi;
    # Prüfen ob es sich um einen endlichen körper handelt??
    Gens := GeneratorsOfGroup(G);
    # F := DefaultFieldOfMatrix(Gens[1]);
    n := DimensionsMat(Gens[1])[1];
    d := Size(Gens);
    hom := fail;

    if d = 1 then
        # Todo.... !!
    fi;
    
    if unitary then
        hom := FrobeniusAutomorphism(F)^(p_exponent/2);
    fi;
    #contains  group element, scalar, (factors of minopol), Basis change to Frobenius (v, vg, vg^2, ...), Frobenius block lengths, number of iterations to compute
    g_res := FindCyclicGroupElementAndScalars(Gens, Lambdas);
    g_inv_frob := Inverse(g_res[4]);
    #CalculateAdjoint := function(mat, mode, hom, n)
    g_star_inv_unscaled := TransposedMat(Inverse(g_res[1]));
    if unitary then
        g_star_inv_unscaled := g_star_inv_unscaled^hom;
    fi;
    #* g_res[2]; 
    # Todo the computation of this can probably be sped up by using the known information about g!!!
    g_star_inv_scaled_frob := FrobeniusNormalForm(g_star_inv_unscaled * g_res[2]);
    frob_base_inv_star := Inverse(g_star_inv_scaled_frob[2]);

    return FindFormspaceInternal(Gens, Lambdas, unitary, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n);
end;

# todo scalars
FindForms := function(G)
    local F, p_exponent, Gens, n, d, first, hom, g_res, g_inv_frob, g_star_inv_unscaled, g_star_inv_scaled_frob, frob_base_inv_star, Lambdas, Out, i;
    Out := [];
    if not IsMatrixGroup(G) then
        # Print("The method only works for matrix groups. \n");
        return;
    fi;
    F := DefaultFieldOfMatrixGroup(G);
    p_exponent := DegreeOverPrimeField(F);
    # Prüfen ob es sich um einen endlichen körper handelt??
    Gens := GeneratorsOfGroup(G);
    # F := DefaultFieldOfMatrix(Gens[1]);
    n := DimensionsMat(Gens[1])[1];
    d := Size(Gens);
    hom := fail;

    if d = 1 then
        # Todo.... !!
    fi;
  
    #contains  group element, scalar, (factors of minopol), Basis change to Frobenius (v, vg, vg^2, ...), Frobenius block lengths, number of iterations to compute
    Lambdas := [];
    for i in [1..d] do
        Add(Lambdas, One(F));
    od;
    g_res := FindCyclicGroupElementAndScalars(Gens, Lambdas);
    ConvertToMatrixRep(g_res[1], F);
    ConvertToMatrixRep(g_res[4], F);
    g_inv_frob := Inverse(g_res[4]);
    #CalculateAdjoint := function(mat, mode, hom, n)
    g_star_inv_unscaled := TransposedMat(Inverse(g_res[1]));
    #* g_res[2]; 
    # Todo the computation of this can probably be sped up by using the known information about g!!!
    g_star_inv_scaled_frob := FrobeniusNormalForm(g_star_inv_unscaled); # hier ist eventuell noch ein bug mit den skalaren
    frob_base_inv_star := Inverse(g_star_inv_scaled_frob[2]);

    Add(Out, FindFormspaceInternal(Gens, Lambdas, false, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n));
    if p_exponent mod 2 = 0 then
        hom := FrobeniusAutomorphism(F)^(p_exponent/2);
        g_star_inv_unscaled := g_star_inv_unscaled^hom;
        g_star_inv_scaled_frob[2] := g_star_inv_scaled_frob[2]^hom;
        frob_base_inv_star := frob_base_inv_star^hom;
        Add(Out, FindFormspaceInternal(Gens, Lambdas, true, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n));
    fi;
    return Out;
end;