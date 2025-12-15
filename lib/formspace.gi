#! @Chapter Formspace
#! For a matrix Group $G$ consisting of $n\times n$ matrices over the field $K$. Let $\Lambda : G \to K^\times $ be a group homomorphism. Suppose $h : K \to K$ is a field automorphism of $K$. For $g\in G$ we define $g^* := h(g^T)$ where the automorphism $h$ is applied entrywise.
#! We define the Formspace $\mathcal{F}_h(G, \Lambda) := \{A\in K^{n\times n} \mid gAg^* = \Lambda(g)A \forall g\in G\}$.
#! @Section Computing the Formspace
#! For details on how the form space is computed SEE: somewhere that does not exist yet

# underlying field F, g is in F^{n\times n}, v is in F^n, coeffs are the coefficients of a polynomial p in F[X]. Returns vp(g).
__FORMSPACE__INTERNAL__EvaluateMatrixPolynomialWithVector := function(F, n, g, v, coeffs)
    local res, i, deg;
    deg := Size(coeffs);
    v := Vector(v);
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


# Given mat in F^{n\times n} F Field, n \in N, mode says whether to apply hom or not. Returns (mat^T)^hom.
# mode = False \iff mat* := mat^tr, mode = True \iff mat* = komplex_konjugiert(mat^tr)
__FORMSPACE__INTERNAL__CalculateAdjoint := function(mat, mode, hom, n, F)
    local transposed, i, j;
    transposed := TransposedMat(mat);
    ConvertToMatrixRep(transposed, F);
    if mode = false then
        return transposed; # this is apparently very slow????? what why?
    fi;
    if mode then
        return transposed^hom;
    fi;
    return fail;
end;

# tries to find a element g \in <Gens> such that the Frobenius Normal form of g has as few blocks as possible. Lambdas describes a group homomorphism induced by Phi : Gens[i] \mapsto Lambdas[i]
# Returns [g, Phi(g), FrobeniusNormalForm(g), nrOfTries]
# nrOfTries contains the number of random elements tested before g was found.
__FORMSPACE__INTERNAL__FindCyclicGroupElementAndScalars := function(Gens, Lambdas)
    local cur_group_element, cur_scalar, known_elements, i, mode, n, res, j, known_scalars, best_known_element_index, best_known_res, best_known_length, mod_elem, g, e;

    known_elements := ShallowCopy(Gens);
    known_scalars := ShallowCopy(Lambdas);

    n := DimensionsMat(Gens[1])[1];
    i := Size(known_elements);
    mod_elem := 1;
    best_known_length := n + 1;
    while i < 25 do # 25 is a magic number. maybe investigate a good number here
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
    # Print("only found element of length ", best_known_length, "\n");
    # TODO: möglichst kurze zyklische modul basis usw....
    if best_known_length = n + 1 then
        return Concatenation([Gens[1], Lambdas[1]], FrobeniusNormalForm(Gens[1]), [-1]);
    fi;
    return Concatenation([known_elements[best_known_element_index], known_scalars[best_known_element_index]], best_known_res, [i]);
end;

# turns the (jn) vector vec in F^{jn} and returns a F^{j times n} matrix
__FORMSPACE__INTERNAL__VectorReorganize := function(vec, j, F, n)
    local A, i;
    A := NullMat(j, n, F);
    ConvertToMatrixRep(A, F);
    for i in [1..j] do
        A[i] := vec{[((i - 1) * n + 1)..(i*n)]};
    od;
    return A;
end;

# turns the F^{j times n} matrix into F^{jn} vector
__FORMSPACE__INTERNAL__MatrixReorganize := function(mat, j, F, n)
    local vec, i;
    vec := ZeroVector(F, j*n);
    for i in [1..j] do
        vec{[((i - 1) * n + 1)..(i*n)]} := mat[i];
    od;
    return vec;
end;

# Evaluates the univariate polynomial p (given as coefficients) in the matrix g \in F^{n\times n}. frob_base = FrobeniusNormalForm(g) must be satisfied and frob_base_inv = Inverse(FrobeniusNormalForm(g)[2]). The reason these two parameters are given, and not computed in the function itself is to not compute FrobeniusNormalForm(g) multiple times when evaluating multiple polynomials in g.
__FORMSPACE__INTERNAL__EvaluatePolynomialWithFrobenius := function(p, g, frob_base, frob_base_inv, F, n) # frob_base_inv useless, right now this is not actually evaluating the polynomial frob_base_inv missing
    local ws, C, i, end_pos, j, k;
    ws := [];
    j := Size(frob_base[3]);
    for k in [1..j] do
    #  function(F, n, g, v, coeffs)
        Add(ws, __FORMSPACE__INTERNAL__EvaluateMatrixPolynomialWithVector(F, n, g, frob_base[2][frob_base[3][k]]{[1..n]}, p));
    od;
    C := NullMat(n, n, F);
    #ConvertToMatrixRep(C, F);
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
    # this mulitiplication might not be needed here 
    # just multiply the found basis in the end when all the condition matrices null spaces got found already (i.e. just return C here)#
    # this could potentially save up to n*d*B matrix matrix multiplications so quite significant!!! definetly investigate!
    ConvertToMatrixRep(C, F);
    return frob_base_inv * C; #to actually evaluate the polynomial
    
    # return frob_base_inv * C;
end;

# This computes (and returns) the matrix \mathcal{P}_{h, u} from the bachelors thesis.
# Here we have u, h \in F^{n\times n}, h_star = h^*, scalar_h = \lambda_h, g_star_inv_scaled = g^{-*}*lambda_g, frob_base = FrobeniusNormalForm(g^{-*}), frob_base_inv = Inverse(FrobeniusNormalForm(g)[2]), frob_base_inv_star = Inverse(frob_base[2]). The reason we to give ten billion parameters is to avoid computing the same matrices multiple times. TODO: better names!!!!!!
__FORMSPACE__INTERNAL__ComputeConditionMatrixFrob := function(u, h, h_star, scalar_h, g_star_inv_scaled, frob_base, frob_base_inv, frob_base_inv_star, F, n)
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
            __FORMSPACE__INTERNAL__EvaluatePolynomialWithFrobenius(coeffs_c{[frob_base[3][i]..b_end]}, g_star_inv_scaled, frob_base, frob_base_inv_star, F, n) * h_star -  __FORMSPACE__INTERNAL__EvaluatePolynomialWithFrobenius(coeffs_f{[frob_base[3][i]..b_end]}, g_star_inv_scaled, frob_base, frob_base_inv_star, F, n);
        # cpol := UnivariatePolynomial(F, coeffs_c);
        # fpol := UnivariatePolynomial(F, coeffs_f);
        # Ps{[((i - 1)*n + 1)..(i*n)]}{[1..n]} :=
        #     cpol(g_star_inv_scaled) * h_star - fpol(g_star_inv_scaled);
        # Print(aua);
    od;
    return Ps;
end;

# Spins the stuff
__FORMSPACE__INTERNAL__FrobSpin := function(Images, spin_elem, frob_base_blocks, n, F)
    local A, j, i, k, end_pos;
    j := Size(frob_base_blocks);
    A := NullMat(n, n, F);
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
    ConvertToMatrixRep(A, F);
    return A;
end;

__FORMSPACE__INTERNAL__FrobSpinAtBlock := function(Image, spin_elem, frob_base_blocks, block_index, n, F)
    local A, j, i, k, end_pos;
    j := Size(frob_base_blocks);
    A := NullMat(n, n, F);
    ConvertToMatrixRep(A, F);
    if block_index = j then
        end_pos := n;
    else
        end_pos := frob_base_blocks[block_index + 1] - 1;
    fi;
    A[frob_base_blocks[block_index]] := Image;
    for k in [(frob_base_blocks[block_index] + 1)..end_pos] do
        A[k] := A[k - 1]*spin_elem;
    od;
    return A;
end;

# checks if (Form^T)^hom = c Form for some c in F and returns d such that d * Form = ((d Form)^T)^hom if possible or if not possibe fail
__FORMSPACE__INTERNAL__ScalarFormIdentifyCheck := function(Form, F, n, hom, p, q)
    local lambda, i, j;
    lambda := fail;
    for i in [1..n] do
        for j in [1..n] do
            if Form[i][j] <> Zero(F) then 
                if Form[j][i] = Zero(F) then
                    return fail;
                fi;
                if lambda <> fail and Form[i][j] * lambda <> hom(Form[j][i]) then
                    return fail;
                fi;
                if lambda = fail then
                    lambda := hom(Form[j][i]) * Inverse(Form[i][j]);
                fi;
            fi;
        od;
    od;
    return RootFFE(F, Inverse(lambda), q - 1);
end;

# find symplectic and symmetric matrices in Forms
# returns bases [[symmetric forms], [symplectic forms]]
__FORMSPACE__INTERNAL__FilterBilinearForms := function(Forms, F, n)
    local transposed_equal_result, form, symmetric_base, symplectic_base, transposed_form, TransposedEqual, symmetric_base_vecs, symplectic_base_vecs, char_2_eqs, sol, out;

    # computes if f = f^T or f = -f^T or none
    # return 0 if f = -f^T
    # return 1 if f = f^T
    # return 2 if f <> f^T and f <> -f^T
    TransposedEqual := function(f, F, n) 
        local i, j, symmetric_possible, symplectic_possible;
        
        if Characteristic(F) mod 2 = 0 then
        # - = + now
            for i in [1..n] do
                for j in [1..n] do
                    if f[i][j] <> f[j][i] then
                        return 2;
                    fi;
                od;
            od;
            return 1;
        fi;
        symmetric_possible := true;
        symplectic_possible := true;
        for i in [1..n] do
            for j in [1..n] do
                if symmetric_possible and f[i][j] <> f[j][i] then
                    symmetric_possible := false;
                fi;
                if symplectic_possible and f[i][j] <> -f[j][i] then
                    symplectic_possible := false;
                fi;
            od;
        od;
        if symmetric_possible then
            return 1;
        fi;
        if symplectic_possible then
            return 0;
        fi;
        return 2;
    end;

    if Size(Forms) = 0 then
        return [];
    fi;
    if Size(Forms) = 1 then
        transposed_equal_result := TransposedEqual(Forms[1], F, n);
        if transposed_equal_result = 2 then
            return [];
        fi;
        if transposed_equal_result = 1 then
            return [[Forms[1]], []];
        fi;
        if transposed_equal_result = 2 then
            return [[], [Forms[1]]];
        fi;

    fi;

    if Characteristic(F) mod 2 <> 0 then
        symmetric_base := MutableBasis(F, [NullMat(n, n, F)]);
        symplectic_base := MutableBasis(F, [NullMat(n, n, F)]);
        for form in Forms do
            transposed_form := TransposedMat(form);
            CloseMutableBasis(symmetric_base, transposed_form + form);
            CloseMutableBasis(symplectic_base, form - transposed_form);
        od;

        # this does not create compressed matrices, rather it returns lists consisting of compressed vectors... TODO: investigate how MutableBasis works on the inside, and maybe change it to use compressed matrices since they are faster to deal with
        symmetric_base_vecs := BasisVectors(ImmutableBasis(symmetric_base));
        symplectic_base_vecs := BasisVectors(ImmutableBasis(symplectic_base));

        if Size(symmetric_base_vecs) + Size(symplectic_base_vecs) <> Size(Forms) then 
            Error("This should not have happend!! there are supposedly ", Size(symmetric_base_vecs), " linearly independent symmetric forms and ", Size(symplectic_base_vecs), " linearly independent symplectic forms, yet the dimension of the formspace is ", Size(Forms), "\n");
        fi;
        return [symmetric_base_vecs, symplectic_base_vecs];

    fi;
    # TODO: ignore the diagonal entries for the symmetric matrices, for symplectic these need to be zero.
    # TODO: solve this in some efficient way that does not formulate this by solving sets of linear equations of matrices. Instead it would be better to gradually consider entries of the matrices. Then use some heuristic to determine we are done and just check wether the resulting matrices are infact symmetric. this should be faster because now worst case we are solving a system of linear equations that consists of n^2 equations and Size(Forms) indeterminates.
    # TODO: maybe do these todos for all characteristics the nice thing about the current algorithm for char F <> 2 is that we do not have to solve many systems of equations, rather just check wether we have the zero matrix

    char_2_eqs := [];
    for form in Forms do
        Add(char_2_eqs, __FORMSPACE__INTERNAL__MatrixReorganize(form - TransposedMat(form), n, F, n));
    od;
    sol := NullspaceMatDestructive(char_2_eqs);
    out := [];
    for form in sol do
        Add(out, __FORMSPACE__INTERNAL__VectorReorganize(form, n, F, n));
    od;
    return out;
end;

# tries to filter the F = GF(q^2) vectorspace generated by <Forms> and return the GF(q) vector space A such that A = <Forms> \cap B where B = {A \in F^{n\times n}, A* = A} TODO: THIS NEEDS SOME FURTHER INVESTIGATION
# there must be a better way to compute these matrices..
__FORMSPACE__INTERNAL__FilterUnitaryForms := function(Forms, F, n, hom)
    local i, j, ent, q, half, O, l, FF, p, tr_form, Base, baseVecs, gf_base, hgf_base;
    if Size(Forms) = 0 then
        return [];
    fi;
    p := Characteristic(F);
    q := p^(DegreeOverPrimeField(F)/2);
    
    if Size(Forms) = 1 then
        #checks if A = A* or A = cA* if A = A* return A, if A = cA* we want to return scalar multiples of A, namely lA for l such that c = l^(1-q) iff c^-1 = l^(q-1)
        # all the solutions then are lA*GF(q) is this correct?? i am not sure if this are indeed all the possible solutions, but it certanly are solutions.
        # Print(ff);
        l := __FORMSPACE__INTERNAL__ScalarFormIdentifyCheck(Forms[1], F, n, hom, p, q);
        if l = fail then
            return [];
        fi;
        
        return [Forms[1] * l];
    fi;
    # this is where it gets interesting
    # Print("ahhh this needs work!\n");
    # kind of okay case?
       
    ## TODO proof if gAg^* = cA for some c it must hold that c \in GF(q) (not in GF(q^2))

    ## we use (A + A*) is a hermitian form if gAg* = A the problem here is that for A, B such that gA*g = A and gB*g = B we may loose information namely it may be the case that 1/2 (A + A*) = c/2 (B + B*) this is annoying.... one thing one could do is check whether these matrices <1/2 (A + A*), 1/2 (B + B*), ...> are lineraly independent. if that is the case, we know that all forms must have been found (but do we???). but what if not? then there might be another form... this is annoying. Then we may add matrices A, B and so on such that <C1,C2, ., D1, D2..> is a basis of F where C1 and so on are hermitian and D1, D2 and so on are not. We may write D1 = A + B where A is hermitian and B is not. we can then try to write B in terms of the other matrices??? does this help... idk :(
    ## for char = 2 this can be a bad idea as it can make the diagonal disaapear..   oh well
    if Characteristic(F) mod 2 <> 0 then
        Base := MutableBasis(GF(q), [NullMat(n, n, F)]);
        gf_base := BasisVectors(Basis(GF(GF(q), 2)))[2];
        hgf_base := hom(gf_base);
        for FF in Forms do
            # l := __FORMSPACE__INTERNAL__ScalarFormIdentifyCheck(FF, F, n, hom, p, q);
            # if l = fail then
            tr_form := TransposedMat(FF^hom);
            CloseMutableBasis(Base, FF + tr_form);
            CloseMutableBasis(Base, gf_base * FF + hgf_base*tr_form);
            # else
                # CloseMutableBasis(Base, l * FF);
            # fi;
        od;

        O := [];
        baseVecs := ImmutableBasis(Base);
        # for FF in baseVecs do
        #     Add(O, HermitianFormByMatrix(FF, F));
        # od;

        # if Size(O) = Size(Forms) then
        #     return O;
        # fi;
        # if Size(baseVecs) = Size(Forms) then
        #     return baseVecs;
        # fi;

        return BasisVectors(baseVecs);
    fi;
    Print("char 2 case is missing!!");
    # # TODO this function definetly needs work!!!
    # Print("Could not find a basis of hermitian Forms. Returned hermitian Forms and a Bigger space of matrices that contains all possible hermitian forms. \n");
    # return [baseVecs, Forms];
end;



# Compute the formspace for cyclic matrix group. 
# TODO: optimizations:
# this function (sometimes) yields a very large formspace 
# to better recognize forms in this case it would be good to add a function that does not do this (since we only care about non degenerate classical forms). to find a bilinear/symplectic/hermitian form we can just compute a invertibe matrix S such that gS = Sg^{-*} (with frobenius normal form) and hope that S + S^*, S - S^* are also inevertible. Then we have found symmetric/symplectic non degenrate forms This seems like a good idea? idk
__FORMSPACE__INTERNAL__CyclicGroupCase := function(Gen, Gen_adjoint_inv_scaled, Lambda, unitary, hom, frob, frob_inv_star_scaled, frob_inv_star_base_change, frob_inv_base_change, F, n)
    # maybe recoginize the trivial group here as a special case
    local p, mat, outspace, i, j, w, OutForms;

    outspace := [];
    for p in frob[1] do#function(p, g, frob_base, frob_base_inv, F, n)
        mat := __FORMSPACE__INTERNAL__EvaluatePolynomialWithFrobenius(CoefficientsOfUnivariatePolynomial(p), Gen_adjoint_inv_scaled * Lambda, frob_inv_star_scaled, frob_inv_star_base_change, F, n);
        Add(outspace, NullspaceMatDestructive(mat));
    od;
    OutForms := [];
    #(Image, spin_elem, frob_base_blocks, block_index, n, F)
    
    for i in [1..Size(outspace)] do
        for w in outspace[i] do
            Add(OutForms, frob_inv_base_change * __FORMSPACE__INTERNAL__FrobSpinAtBlock(w, Gen_adjoint_inv_scaled, frob[3], i, n, F));
        od;
    od;
    return OutForms;
end;

# Returns formspace preserved by the group <Gens> modulo Lambdas. Unitary says wheter to look for unitary forms or not. hom can be the Field Automorphism of order two. g_res = [g, Lambda_g, ## The elements of ## FrobeniusNormalForm(g)]. Where g is randomly determenied. g_inv_frob = Inverse(FrobeniusNormalForm(g)[2]). g_star_inv_scaled = g^{-*} * Lambda_g, g_star_inv_scaled_frob = FrobeniusNormalForm(g^{-*} * Lambda_g), frob_base_inv_star = Inverse(g_star_inv_scaled_frob[2]). d = Size(Gens), F is base field and n is the matrix dimension. 
__FORMSPACE__INTERNAL__FormspaceInternal := function(Gens, Lambdas, unitary, hom, g_res, g_inv_frob, g_star_inv_scaled, g_star_inv_scaled_frob, frob_base_inv_star, d, F, n)

    local k, i, W, first, j, Cond, Conds, h_star, h, w, O, needs_checking, A, failed_check, nspace, vec;
    first := true;
    needs_checking := false; ## maybe remove?

    for i in [1..d] do
        if needs_checking then
            break;
        fi;
        h := Gens[i];
        h_star := __FORMSPACE__INTERNAL__CalculateAdjoint(h, unitary, hom, n, F);
        for j in [1..n] do
            # vec := RandomVector(F, n);
            vec := g_res[4][j];
            # Display(vec);
            Conds := 
                __FORMSPACE__INTERNAL__ComputeConditionMatrixFrob(vec, h, h_star, Lambdas[i], g_star_inv_scaled, g_star_inv_scaled_frob, g_inv_frob, frob_base_inv_star, F, n);
            
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
                # technically we still need checking here..
            fi;
            if Size(nspace) = 1 then
                needs_checking := true;
                break;
            fi;
        od;
    od;
    O := [];
    # W := frob_base_inv_star * W;
    # Print(WWW);
    for w in W do
        A := g_inv_frob * __FORMSPACE__INTERNAL__FrobSpin(__FORMSPACE__INTERNAL__VectorReorganize(w, Size(g_res[5]), F, n), g_star_inv_scaled, g_res[5], n, F);
        ## This whole checking should be removed, it is stupid to always check!!!
        if needs_checking then
            failed_check := false;
            for i in [1..d] do
                if not failed_check and Gens[i] * A * __FORMSPACE__INTERNAL__CalculateAdjoint(Gens[i], unitary, hom, n, F) <> Lambdas[i] * A then
                    failed_check := true;
                fi;
            od;
            if not failed_check then
                Add(O, A);
                # if not unitary then
                #     Add(O, BilinearFormByMatrix(A, F));
                # else
                #     Add(O, A);
                # fi;
            fi;
        fi;
        if not needs_checking then
            # if not unitary then
            #     # throws errors since the forms package only accepts symmetric bilinear forms... i think this should be changed..
            #     # or should i solve for symmetric matrices as a system of linear equations..?
            #     Add(O, BilinearFormByMatrix(A, F));
            # else
            #     Add(O, A);
            # fi;
            Add(O, A);
        fi;
    od;
    # if unitary then
    #     return __FORMSPACE__INTERNAL__FilterUnitaryForms(O, F, n, hom);
    # fi;
    return O;
end;


#! @Arguments G, L, unitary
#! @Returns a basis of $\mathcal{F}_h(G, \Lambda)$
#! @Description
#!  <A>G</A> is a finitely generated matrix group over a finite field $K$. 
#!  L is a list of scalars corresponding to the group generators. I.e. $\Lambda(G_i)$ = <A>L[i]</A> for $i = 1, \dots, $<A>Size(GeneratorsOfGroup(G))</A>. 
#!  If unitary is true, the function uses $h(x) = x$ for all $x \in K$. 
#!  If unitary is false, the function tries to use a field automorphism of order two for $h$. If no such field automorphism exists, it will return [].
InstallMethod(PreservedFormspace,
    "for matrix group over finite field, with given scalars, and search for unitary forms", [IsMatrixGroup, IsVector and IsFFECollection, IsBool],
    function(G, Lambdas, unitary) 
        local F, p_exponent, Gens, n, d, first, hom, g_res, g_inv_frob, g_star_inv_unscaled, g_star_inv_scaled_frob, frob_base_inv_star, frob, frob_inv_star, Gen, Gen_adjoint, Gen_adjoint_inv_scaled, frob_inv_star_scaled, frob_inv_star_base_change;
        F := DefaultFieldOfMatrixGroup(G);
        p_exponent := DegreeOverPrimeField(F);
        if unitary and (p_exponent mod 2 <> 0) then
            # Print("Field does not admit field automorphism of order two!");
            return [];
        fi;
        # Prüfen ob es sich um einen endlichen körper handelt??
        Gens := GeneratorsOfGroup(G);
        # F := DefaultFieldOfMatrix(Gens[1]);
        n := DimensionsMat(Gens[1])[1];
        d := Size(Gens);
        hom := fail;

        
        if unitary then
            hom := FrobeniusAutomorphism(F)^(p_exponent/2);
        fi;
        
        if d = 1 then
            Gen := Gens[1];
            Gen_adjoint := __FORMSPACE__INTERNAL__CalculateAdjoint(Gen, unitary, hom, n, F);
            Gen_adjoint_inv_scaled := Lambda[1]*Inverse(Gen_adjoint);
            frob := FrobeniusNormalForm(Gen);
            frob_inv_star_scaled := FrobeniusNormalForm(Gen_adjoint_inv_scaled);
            frob_inv_star_base_change := Inverse(frob_inv_star_scaled[2]);

            return __FORMSPACE__INTERNAL__CyclicGroupCase(Gen, Gen_adjoint_inv_scaled, Lambda[1], unitary, hom, frob, frob_inv_star_scaled, frob_inv_star_base_change, Inverse(frob[2]), F, n);
        fi;
        #contains  group element, scalar, (factors of minopol), Basis change to Frobenius (v, vg, vg^2, ...), Frobenius block lengths, number of iterations to compute
        g_res := __FORMSPACE__INTERNAL__FindCyclicGroupElementAndScalars(Gens, Lambdas);
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

        return __FORMSPACE__INTERNAL__FormspaceInternal(Gens, Lambdas, unitary, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n);
    end
);

#! @Arguments G,
#! @Returns a basis of $\mathcal{F}_{id}(G, \Lambda_1)$ and $\mathcal{F}_h(G, \Lambda_1)$
#! @Description
#!  <A>G</A> is a finitely generated Matrix group over a finite field $K$. 
#! Here $\Lambda_1(g) := 1$ for all $g\in G$ and $id(x) := x$ for all $x \in K$. Furthermore $h$ is a field automorphism of order two, admitted by $h$. If no such field automorphism exists, the second returned basis is empty.
InstallMethod(PreservedFormspace,
    "for matrix group (finds bilinear/unitary forms modulo 1)",
    [IsMatrixGroup],
    function(G)
        local F, p_exponent, Gens, n, d, first, hom, g_res, g_inv_frob, g_star_inv_unscaled, g_star_inv_scaled_frob, frob_base_inv_star, Lambdas, Out, i, Gen, Gen_adjoint, Gen_adjoint_inv_scaled, frob, frob_inv_star_scaled, frob_inv_star_base_change;
        Out := [];
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
            Gen := Gens[1];
            Gen_adjoint := __FORMSPACE__INTERNAL__CalculateAdjoint(Gen, false, fail, n, F);
            Gen_adjoint_inv_scaled := One(F)*Inverse(Gen_adjoint); # scaling happens here!!
            frob := FrobeniusNormalForm(Gen);
            frob_inv_star_scaled := FrobeniusNormalForm(Gen_adjoint_inv_scaled);
            frob_inv_star_base_change := Inverse(frob_inv_star_scaled[2]);
            g_inv_frob := Inverse(frob[2]);

            Add(Out, __FORMSPACE__INTERNAL__CyclicGroupCase(Gen, Gen_adjoint_inv_scaled, One(F), false, fail, frob, frob_inv_star_scaled, frob_inv_star_base_change, g_inv_frob, F, n));

            if p_exponent mod 2 = 0 then
                hom := FrobeniusAutomorphism(F)^(p_exponent/2);
                Gen_adjoint_inv_scaled := Gen_adjoint_inv_scaled^hom;
                frob_inv_star_base_change := frob_inv_star_base_change^hom;
                frob_inv_star_scaled[2] := frob_inv_star_scaled[2]^hom;

                Add(Out, __FORMSPACE__INTERNAL__CyclicGroupCase(Gen, Gen_adjoint_inv_scaled, One(F), true, hom, frob, frob_inv_star_scaled, frob_inv_star_base_change, g_inv_frob, F, n));
            else 
                Add(Out, []);
            fi;
            return Out;
        fi;
    
        Lambdas := [];
        for i in [1..d] do
            Add(Lambdas, One(F));
        od;
        #contains  group element, scalar, (factors of minopol), Basis change to Frobenius (v, vg, vg^2, ...), Frobenius block lengths, number of iterations to compute
        g_res := __FORMSPACE__INTERNAL__FindCyclicGroupElementAndScalars(Gens, Lambdas);
        ConvertToMatrixRep(g_res[1], F);
        ConvertToMatrixRep(g_res[4], F);
        g_inv_frob := Inverse(g_res[4]);
        #CalculateAdjoint := function(mat, mode, hom, n)
        g_star_inv_unscaled := TransposedMat(Inverse(g_res[1]));
        #* g_res[2]; 
        # Todo the computation of this can probably be sped up by using the known information about g!!!
        g_star_inv_scaled_frob := FrobeniusNormalForm(g_star_inv_unscaled); # hier ist eventuell noch ein bug mit den skalaren
        frob_base_inv_star := Inverse(g_star_inv_scaled_frob[2]);

        Add(Out, __FORMSPACE__INTERNAL__FormspaceInternal(Gens, Lambdas, false, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n));
        if p_exponent mod 2 = 0 then
            # is ^hom actually cheaper than computing the frobenius normal form?? investigate! certainly makes the code ugly... oh well
            hom := FrobeniusAutomorphism(F)^(p_exponent/2);
            g_star_inv_unscaled := g_star_inv_unscaled^hom;
            g_star_inv_scaled_frob[2] := g_star_inv_scaled_frob[2]^hom;
            frob_base_inv_star := frob_base_inv_star^hom;
            Add(Out, __FORMSPACE__INTERNAL__FormspaceInternal(Gens, Lambdas, true, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n));
        else
            Add(Out, []);
        fi;
        return Out;
    end
);

#! @Arguments Forms, Field, unitary
#! @Returns basis of the spaces of symmetric/symplectic matrices or a basis of the hermitian matrices contained in Forms
#! @Description
#! In the case where unitary is false, this will return a list that contains to lists of matrices. The first is a basis of the symmetric matrices, the second is a basis of the symplectic matrices. Be carefull: If the characteristic of the field is 2, then symplectic matrices and symmetric matrices are the same. Hence only one basis will be returned. If unitary is false this will return a basis of the hermitian matrices.
#! The reason the the field must be given as an argument, is so that the Field automorphims of order two, which is used to compute the adjoined, is specified.
InstallMethod(FilterFormspace, "for list of matrices, F finite field, bool hermitian", [IsList, IsFinite and IsField, IsBool], function(Forms, F, unitary)
    local n, hom, p_exponent;
    if Size(Forms) = 0 then
        return [];
    fi;
    n := NrRows(Forms[1]);

    if Size(Forms) = n^2 then
        # TODO: special case where we should just return a precomputed basis... 
    fi;

    if not unitary then
        return __FORMSPACE__INTERNAL__FilterBilinearForms(Forms, F, n);
    else
        p_exponent := DegreeOverPrimeField(F);
        if p_exponent mod 2 <> 0 then 
            Error("The given Field ", F, " must admit a field automorphism of order two for unitary = true\n");
            return [];
        fi;
        hom := FrobeniusAutomorphism(F)^(p_exponent/2);
        return __FORMSPACE__INTERNAL__FilterUnitaryForms(Forms, F, n, hom);
    fi;
end);