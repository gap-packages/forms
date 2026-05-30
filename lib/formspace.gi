#! @Chapter Formspace
#! For a matrix Group $G$ consisting of $n\times n$ matrices over the field $K$. Let $\Lambda : G \to K^\times $ be a group homomorphism. Suppose $h : K \to K$ is a field automorphism of $K$. For $g\in G$ we define $g^* := h(g^T)$ where the automorphism $h$ is applied entrywise.
#! We define the Formspace $\mathcal{F}_h(G, \Lambda) := \{A\in K^{n\times n} \mid gAg^* = \Lambda(g)A \forall g\in G\}$.
#! @Section Computing the Formspace
#! For details on how the form space is computed SEE: somewhere that does not exist yet

# underlying field F, g is in F^{n\times n}, v is in F^n, coeffs are the coefficients of a polynomial p in F[X]. Returns vp(g).
FORMS_EvaluateMatrixPolynomialWithVector := function(F, n, g, v, coeffs)
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


# Given mat in F^{n\times n} F Field, n \in N, mode says whether to apply hom or not. Returns (mat^T)^hom.
# mode = False \iff mat* := mat^tr, mode = True \iff mat* = komplex_konjugiert(mat^tr)
FORMS_CalculateAdjoint := function(mat, mode, hom, n, F)
    local transposed, i, j;
    transposed := TransposedMat(mat);
    # ConvertToMatrixRep(transposed, F);
    if mode = false then
        return transposed;
    fi;
    if mode then
        return transposed^hom;
    fi;
    return fail;
end;

# tries to find a element g \in <Gens> such that the Frobenius Normal form of g has as few blocks as possible. Lambdas describes a group homomorphism induced by Phi : Gens[i] \mapsto Lambdas[i]
# Returns [g, Phi(g), FrobeniusNormalForm(g), nrOfTries]
# nrOfTries contains the number of random elements tested before g was found.
FORMS_FindCyclicGroupElementAndScalars := function(Gens, Lambdas, n)
    local known_elements, known_scalars, i, j, run, runs_max, best_known_length, best_known_res, best_known_scalar, best_known_element, g, lamb, s, d, res, interval;

    d := Size(Gens);

    # no point to random group generators in cyclic groups.
    if d = 1 then 
        return Concatenation([Gens[1], Lambdas[1]], FrobeniusNormalForm(Gens[1]), [-1]); 
    fi;
    runs_max := 3*d + 15; # pretty random idk what to put here

    known_elements := ShallowCopy(Gens);
    known_scalars := ShallowCopy(Lambdas);

    best_known_length := n + 1;

    for run in [1..runs_max] do
        i := PseudoRandom([1..d]);
        interval := [1..d];
        Remove(interval, i);
        j := PseudoRandom(interval);
        # TODO: add inverses? seems fine without
        g := known_elements[i]*known_elements[j];
        lamb := known_scalars[i]*known_scalars[j];
        res := FrobeniusNormalForm(g);
        s := Size(res[3]);
        if s = 1 then
            return Concatenation([g, lamb], res, [run]);
        fi;
        if s < best_known_length then
            best_known_length := s;
            best_known_element := g;
            best_known_scalar := lamb;
            best_known_res := res;
        fi;
        known_elements[i] := g;
        known_scalars[i] := lamb;
    od;
   
    if best_known_length = n + 1 then
        return Concatenation([Gens[1], Lambdas[1]], FrobeniusNormalForm(Gens[1]), [-1]);
    fi;
    return Concatenation([best_known_element, best_known_scalar], best_known_res, [runs_max]);
end;

# turns the (jn) vector vec in F^{jn} and returns a F^{j times n} matrix
FORMS_VectorReorganize := function(vec, j, F, n)
    local A, i;
    A := NullMat(j, n, F);
    # A := ZeroMatrix(F, j, n);
    # ConvertToMatrixRep(A, F);
    for i in [1..j] do
        A[i] := vec{[((i - 1) * n + 1)..(i*n)]};
    od;
    return A;
end;

# turns the F^{j times n} matrix into F^{jn} vector
FORMS_MatrixReorganize := function(mat, j, F, n)
    local vec, i;
    vec := ZeroVector(F, j*n);
    for i in [1..j] do
        vec{[((i - 1) * n + 1)..(i*n)]} := mat[i];
    od;
    return vec;
end;

# given j times n matrices over field F, this tries to compute a linear combination of the matrices such that their sum is zero.
# TODO: maybe make this funciton more efficient by only considereing some vectors, a better optimatizion however, would be to not even generate the unused equations
FORMS_SolveMatrixSystem := function(mats, j, n, F)
    local eqs, mat, sol, out;
    eqs := [];
    for mat in mats do 
        Add(eqs, FORMS_MatrixReorganize(mat, j, F, n));
    od;
    return NullspaceMatDestructive(eqs);
end;

# Spins a block in the frobenius normal form, with the images.
FORMS_FrobSpin := function(Images, spin_elem, frob_base_blocks, n, F)
    local A, j, i, k, end_pos;
    j := Size(frob_base_blocks);
    A := NullMat(n, n, F);
    # A := ZeroMatrix(F, n, n);
    CopySubMatrix(Images, A, [1..j], frob_base_blocks, [1..n], [1..n]);
    for i in [1..j] do
        if i = j then
            end_pos := n;
        else
            end_pos := frob_base_blocks[i + 1] - 1;
        fi;
        # causes error if A is matrix obj and Images is plain list... :(( 
        # A[frob_base_blocks[i]] := Images[i];
        for k in [(frob_base_blocks[i] + 1)..end_pos] do
            CopySubMatrix([A[k - 1]*spin_elem], A, [1], [k], [1..n], [1..n]);
            # A[k] := A[k - 1]*spin_elem;
        od;
    od;
    return A;
end;

# Evaluates the univariate polynomial p (given as coefficients) in the matrix g \in F^{n\times n}. frob_base = FrobeniusNormalForm(g) must be satisfied and frob_base_inv = Inverse(FrobeniusNormalForm(g)[2]). The reason these two parameters are given, and not computed in the function itself is to not compute FrobeniusNormalForm(g) multiple times when evaluating multiple polynomials in g.
FORMS_EvaluatePolynomialWithFrobenius := function(p, g, frob_base, frob_base_inv, F, n) 
    local ws, C, i, end_pos, j, k;
    ws := [];
    j := Size(frob_base[3]);
    for k in [1..j] do
        Add(ws, FORMS_EvaluateMatrixPolynomialWithVector(F, n, g, frob_base[2][frob_base[3][k]]{[1..n]}, p));
    od;
    C := FORMS_FrobSpin(ws, g, frob_base[3], n, F);
    # TODO: this function is used to build the condition matrices, however since we are only interested in the kernels of these matrices and C is an invertible matrix we can probably eliminate this multiplication and do it at a later stage?
    return frob_base_inv * C;
end;

# This computes (and returns) the matrix \mathcal{P}_{h, u} from the bachelors thesis.
# Here we have u, h \in F^{n\times n}, h_star = h^*, scalar_h = \lambda_h, g_star_inv_scaled = g^{-*}*lambda_g, frob_base = FrobeniusNormalForm(g^{-*}), frob_base_inv = Inverse(FrobeniusNormalForm(g)[2]), frob_base_inv_star = Inverse(frob_base[2]). The reason we to give so many parameters is to avoid computing the same matrices multiple times. TODO: better names!!!!!!
FORMS_ComputeConditionMatrixFrob := function(u, h, h_star, scalar_h, g_star_inv_scaled, frob_base, frob_base_inv, frob_base_inv_star, F, n)
    local coeffs_c, coeffs_f, Ps, i, j, b_end, b_start, cpol, fpol;
    coeffs_c := (u * h) * frob_base_inv;
    coeffs_f := (u * frob_base_inv) * scalar_h;
    j := Size(frob_base[3]);
    Ps := NullMat(n * j, n, F);
    # Ps := ZeroMatrix(F, n*j, n);
    for i in [1..j] do
        if i = j then
            b_end := n;
        else
            b_end := frob_base[3][i + 1] - 1;
        fi;
        Ps{[((i - 1)*n + 1)..(i*n)]}{[1..n]} :=
            FORMS_EvaluatePolynomialWithFrobenius(coeffs_c{[frob_base[3][i]..b_end]}, g_star_inv_scaled, frob_base, frob_base_inv_star, F, n) * h_star - FORMS_EvaluatePolynomialWithFrobenius(coeffs_f{[frob_base[3][i]..b_end]}, g_star_inv_scaled, frob_base, frob_base_inv_star, F, n);
    od;
    return Ps;
end;

## similar to FORMS_FrobSpin but only spins one block namely that of block_index
FORMS_FrobSpinAtBlock := function(Image, spin_elem, frob_base_blocks, block_index, n, F)
    local A, j, i, k, end_pos;
    j := Size(frob_base_blocks);
    A := NullMat(n, n, F);
    # A := ZeroMatrix(F, n, n);
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

# find symplectic and symmetric matrices in Forms
# returns bases [[symmetric forms], [symplectic forms]] for char <> 2 and [symmetric forms] for char = 2
FORMS_FilterBilinearForms := function(Forms, F, n)
    local transposed_equal_result, form, symmetric_base, symplectic_base, transposed_form, TransposedEqual, symmetric_base_vecs, symplectic_base_vecs, char_2_eqs, sol, out, mat, tmp, i, s;

    if Size(Forms) = 0 then
        return [];
    fi;
    
    if Characteristic(F) <> 2 then
        if Size(Forms) = 1 then
            if FORMS_IsSymplecticMatrix(Forms[1], F) then
                return [[], [Forms[1]]];
            fi;
            if FORMS_IsSymmetricMatrix(Forms[1]) then
                return [[Forms[1]], []];
            fi;
            return [];
        fi;


        # maybe not use these mutable bases
        symmetric_base := MutableBasis(F, [NullMat(n, n, F)]);
        symplectic_base := MutableBasis(F, [NullMat(n, n, F)]);
        # symmetric_base := MutableBasis(F, [ZeroMatrix(F, n, n)]);
        # symplectic_base := MutableBasis(F, [ZeroMatrix(F, n, n)]);
        # symmetric_base := MutableBasis(F, [], ZeroVector(F, n));
        # symplectic_base := MutableBasis(F, [], ZeroVector(F, n));
        for form in Forms do
            transposed_form := TransposedMat(form);
            CloseMutableBasis(symmetric_base, transposed_form + form);
            CloseMutableBasis(symplectic_base, form - transposed_form);
        od;

        symmetric_base_vecs := BasisVectors(ImmutableBasis(symmetric_base));
        symplectic_base_vecs := BasisVectors(ImmutableBasis(symplectic_base));

        if Size(symmetric_base_vecs) + Size(symplectic_base_vecs) <> Size(Forms) then 
            Error("This should not have happend!! there are supposedly ", Size(symmetric_base_vecs), " linearly independent symmetric forms and ", Size(symplectic_base_vecs), " linearly independent symplectic forms, yet the dimension of the formspace is ", Size(Forms), "\n");
        fi;
        return [symmetric_base_vecs, symplectic_base_vecs];
    fi;
    # TODO: solve this in some efficient way that does not formulate this by solving sets of linear equations of matrices. Instead it would be better to gradually consider entries of the matrices. Then use some heuristic to determine we are done and just check wether the resulting matrices are infact symmetric. this should be faster because now worst case we are solving a system of linear equations that consists of n^2 equations and Size(Forms) indeterminates.
    # TODO: maybe do these todos for all characteristics?

    if Size(Forms) = 1 then
        if FORMS_IsSymplecticMatrix(Forms[1], F) then
            return [Forms[1]];
        else 
            return [];
        fi;
    fi;
    
    char_2_eqs := [];
    for form in Forms do
        Add(char_2_eqs, form - TransposedMat(form));
    od;
    out := [];
    sol := FORMS_SolveMatrixSystem(char_2_eqs, n, n, F);
    for s in sol do
        tmp := Forms[1]*s[1];
        for i in [2..Size(s)] do
            tmp := tmp + s[i]*Forms[i];
        od;
        Add(out, tmp);
    od;
    return out;
end;

# tries to filter the F = GF(q^2) vectorspace generated by <Forms> and return the GF(q) vector space A such that A = <Forms> \cap B where B = {A \in F^{n\times n}, A* = A}
FORMS_FilterUnitaryForms := function(Forms, F, n, hom)
    local i, j, ent, q, half, O, l, FF, p, tr_form, Base, baseVecs, gf_base, hgf_base, mat, small_field, field_aut_mat, gf_base_vecs, big_aut_mat, to_smaller_field_matrix, transpose_weird_mat, apply_aut_to_mat, eqs, s, form, form_changed, sol, out, tmp, big_field, changed_forms, apply_coefficient_and_get_smaller_matrix, square_base_entry_rep, multiply_with_scalar, form_changed_mult, form_changed_mult_star, form_changed_star, FORMS_ScalarFormIdentifyCheck;
    if Size(Forms) = 0 then
        return [];
    fi;
    p := Characteristic(F);
    q := p^(DegreeOverPrimeField(F)/2);
    
    # checks if (Form^T)^hom = c Form for some c in F and returns d such that d * Form = ((d Form)^T)^hom if possible or if not possibe fail
    FORMS_ScalarFormIdentifyCheck := function(Form, F, n, hom, p, q)
        local lambda, i, j;
        lambda := fail;
        for i in [1..n] do
            for j in [1..n] do
                if not IsZero(Form[i, j]) then 
                    if IsZero(Form[j, i]) then
                        return fail;
                    fi;
                    if lambda <> fail and Form[i, j] * lambda <> hom(Form[j, i]) then
                        return fail;
                    fi;
                    if lambda = fail then
                        lambda := hom(Form[j, i]) * Inverse(Form[i, j]);
                    fi;
                fi;
            od;
        od;
        return RootFFE(F, Inverse(lambda), q - 1);
    end;

    if Size(Forms) = 1 then
        #checks if A = A* or A = cA* if A = A* return A, if A = cA* we want to return scalar multiples of A, namely lA for l such that c = l^(1-q) iff c^-1 = l^(q-1)
        # all the solutions then are lA*GF(q) is this correct?? i am not sure if this are indeed all the possible solutions, but it certanly are solutions.
        l := FORMS_ScalarFormIdentifyCheck(Forms[1], F, n, hom, p, q);
        if l = fail then
            return [];
        fi;
        
        return [Forms[1] * l];
    fi;
       
    ## Lemma: If A = A* is a form preserved modulo scalar c in GF(q^2) i.e. gAg* = cA for group elements g and A is not 0. then c in GF(q).

    ## for char = 2 this can not yield all hermitian matrices as it deletes the diagonal.
    if p <> 2 then
        Base := MutableBasis(GF(q), [NullMat(n, n, GF(q))]);
        # Base := MutableBasis(GF(q), [ZeroMatrix(GF(q), n, n)]);
        
        # Base := MutableBasis(GF(q), [], ZeroVector(GF(q), n));
        gf_base := BasisVectors(Basis(GF(GF(q), 2)))[2];
        hgf_base := hom(gf_base);
        for FF in Forms do
            tr_form := TransposedMat(FF^hom);
            CloseMutableBasis(Base, FF + tr_form);
            CloseMutableBasis(Base, gf_base * FF + hgf_base*tr_form);
        od;

        O := [];
        baseVecs := BasisVectors(ImmutableBasis(Base));

        return baseVecs;
    fi;

    # the idea for char 2 is to solve the semiliner system of equations. we take a GF(q) basis of GF(q^2) namely <1, delta> and express n times n matrices with this basis. 

    # TODO: this function is potentially really slow, it would be much better two only take a few equations and constain the problem instead of taking the entire n times n matrix. One such example would be the form space preserved by G := Group(SU(200, 2^2).1) then Size(Forms) = 38420 consisting of 200x200 matrices. For cases such as this it might also be good to just have some function that tries to find one random non-deg form instead of generating random ones

    small_field := GF(q);
    big_field := GF(q^2);
    gf_base := Basis(GF(small_field, 2));
    gf_base_vecs :=  BasisVectors(gf_base);
    # expresses n times n matrices as 2n times n matrices where 2 times 1 collumn vectors contain the coefficients ascoiciated with the basis gf_base
    to_smaller_field_matrix := function(basis_of_field, small_field_, mat, n, c)
        local outmat, i, j;
        outmat := NullMat(2 * n, c, small_field_);
        # outmat := ZeroMatrix(small_field_, 2 * n, c);
        for i in [1..n] do
            for j in [1..c] do
                outmat{[(2*i - 1)..(2*i)]}[j] := Coefficients(basis_of_field, mat[i][j]);
            od;
        od;
        return outmat;
    end;
    # this matrix applies the field automorphism hom to a element expressed in the basis over the smaller field
    field_aut_mat := to_smaller_field_matrix(gf_base, small_field, [[hom(gf_base_vecs[1]), hom(gf_base_vecs[2])]], 1, 2); 
    apply_aut_to_mat := function(mat, n, aut_mat)
        local outmat, i, j;
        outmat := 1*mat; # ugly hack to force gap to make a new copy of the matrix.
        for i in [1..n] do
            for j in [1..n] do
                outmat{[(2*i - 1)..(2*i)]}[j] := aut_mat * outmat{[(2* i - 1)..(2*i)]}[j];
            od;
        od;
        return outmat;
    end;

    # transposes 2n times n matrix (as if it were a n times n matrix)
    transpose_weird_mat := function(mat, n)
        local outmat, i, j;
        outmat := 1*mat; # hack to force gap to copy the matrix.
        for i in [1..n] do
            for j in [1..n] do
                outmat{[(2*i-1)..(2*i)]}[j] := mat{[(2*j-1)..(2*j)]}[i];
            od;
        od;
        return outmat;
    end;

    square_base_entry_rep := Coefficients(gf_base, gf_base_vecs[2]^2);
    # multiplies a 2n times n matrix over F_q with a scalar from F_q^2
    multiply_with_scalar := function(scalar,gf_base, square_base_entry_rep, mat, n)
        local smat, srep;
        srep := Coefficients(gf_base, scalar);
        smat := [[srep[1], srep[2] * square_base_entry_rep[1]], 
        [srep[2], srep[1] + srep[2] * square_base_entry_rep[2]]];
        return apply_aut_to_mat(mat, n, smat);
    end;

    eqs := [];
    # build linear equations that can now be solved over F_q
    for form in Forms do
        form_changed := to_smaller_field_matrix(gf_base, small_field, form, n, n);
        form_changed_star := apply_aut_to_mat(transpose_weird_mat(form_changed, n), n, field_aut_mat);
        form_changed_mult := multiply_with_scalar(gf_base_vecs[2], gf_base, square_base_entry_rep, form_changed, n);
        form_changed_mult_star := multiply_with_scalar(hom(gf_base_vecs[2]), gf_base, square_base_entry_rep, form_changed, n);

        Add(eqs, form_changed - form_changed_star);
        Add(eqs, form_changed_mult - form_changed_mult_star);
    od;

    sol := FORMS_SolveMatrixSystem(eqs, 2*n, n, small_field);
    out := [];
    for s in sol do
        tmp := (s[1]*gf_base_vecs[1] + s[2]*gf_base_vecs[2])*Forms[1];
        for i in [2..Size(Forms)] do
            tmp := tmp + (s[2*i - 1]*gf_base_vecs[1] + s[2*i]*gf_base_vecs[2])*Forms[i];
        od;
        Add(out, tmp);
    od;
    return out;
end;

# Compute the formspace for cyclic matrix group. 
# TODO: optimizations:
# this function (sometimes) yields a very large formspace 
# to better recognize forms in this case it would be good to add a function that does not do this (since we only care about non degenerate classical forms). to find a bilinear/symplectic/hermitian form we can just compute a invertibe matrix S such that gS = Sg^{-*} (with frobenius normal form) and hope that S + S^*, S - S^* are also inevertible. Then we have found symmetric/symplectic non degenrate forms This seems like a good idea? idk
FORMS_CyclicGroupCase := function(Gen, Gen_adjoint_inv_scaled, Lambda, unitary, hom, frob, frob_inv_star_scaled, frob_inv_star_base_change, frob_inv_base_change, F, n)
    # maybe recoginize the trivial group here as a special case
    local p, mat, outspace, i, j, w, OutForms;

    outspace := [];
    for p in frob[1] do
        mat := FORMS_EvaluatePolynomialWithFrobenius(CoefficientsOfUnivariatePolynomial(p), Gen_adjoint_inv_scaled * Lambda, frob_inv_star_scaled, frob_inv_star_base_change, F, n);
        Add(outspace, NullspaceMatDestructive(mat));
    od;
    OutForms := [];
    
    for i in [1..Size(outspace)] do
        for w in outspace[i] do
            Add(OutForms, frob_inv_base_change * FORMS_FrobSpinAtBlock(w, Gen_adjoint_inv_scaled, frob[3], i, n, F)); # this can be used to bound the rank very cheaply, it is smaller than the lenght of the ith frobenius block, furthermore we use the different frobenius blocks, to build a non-deg form we should use a form from each block and add them??
        od;
    od;
    return OutForms;
end;


# builds the forms from image vectors and frobenius normal forms and so on, might needs to check if the forms are actually preserved forms
FORMS_ReturnFormspace := function(needs_checking, W, g_res, g_star_inv_scaled, Lambdas, unitary, hom, Gens, d, F, g_inv_frob, n)
    local O, w, A, i;
    # no kernel, return empty
    O := [];
    for w in W do
        A := g_inv_frob * FORMS_FrobSpin(FORMS_VectorReorganize(w, Size(g_res[5]), F, n), g_star_inv_scaled, g_res[5], n, F);
        if needs_checking then
            for i in [1..d] do
                #todo: store adjoints so they do not have to be reused!
                if Gens[i] * A * FORMS_CalculateAdjoint(Gens[i], unitary, hom, n, F) <> Lambdas[i] * A then
                    if Size(W) = 1 then
                        return []; #one dimensional case no further computation needed
                    fi;
                    return false;
                fi;
            od;
            Add(O, A);
        else
            Add(O, A);
        fi;
    od;
    return O;
end;

# Returns formspace preserved by the group <Gens> modulo Lambdas. Unitary says wheter to look for unitary forms or not. hom can be the Field Automorphism of order two. g_res = [g, Lambda_g, ## The elements of ## FrobeniusNormalForm(g)]. Where g is randomly determenied. g_inv_frob = Inverse(FrobeniusNormalForm(g)[2]). g_star_inv_scaled = g^{-*} * Lambda_g, g_star_inv_scaled_frob = FrobeniusNormalForm(g^{-*} * Lambda_g), frob_base_inv_star = Inverse(g_star_inv_scaled_frob[2]). d = Size(Gens), F is base field and n is the matrix dimension. 
FORMS_FormspaceInternal := function(Gens, Lambdas, unitary, hom, g_res, g_inv_frob, g_star_inv_scaled, g_star_inv_scaled_frob, frob_base_inv_star, d, F, n)
    local k, i, W, first, j, Cond, Conds, h_star, h, w, O, A, failed_check, nspace, vec, stagnatiton, stagnation_max, old_kernel_size, nspace_size;
    first := true;
   
    ## after computing 5 kernels and the formspace not getting smaller, assume that all forms have been found.
    old_kernel_size := n + 1;
    stagnatiton := 0;
    stagnation_max := 2; # test three matrices
    for i in [1..d] do
        h := Gens[i];
        h_star := FORMS_CalculateAdjoint(h, unitary, hom, n, F);
        for j in [1..n] do
            # vec := RandomVector(F, n); # does not seem to be better
            vec := g_res[4][j];
            Conds := 
                FORMS_ComputeConditionMatrixFrob(vec, h, h_star, Lambdas[i], g_star_inv_scaled, g_star_inv_scaled_frob, g_inv_frob, frob_base_inv_star, F, n);
            
            if not first then
                nspace := NullspaceMat(W * Conds);
                nspace_size := Size(nspace);
                if nspace_size = 0 then 
                    return [];
                fi;
                W := nspace * W;
            fi;
            if first then
                nspace := NullspaceMat(Conds);
                nspace_size := Size(nspace);
                if nspace_size = 0 then 
                    return [];
                fi;
                W := nspace;
                first := false;
            fi;
            if nspace_size = 1 then
                break;
            fi;
            if old_kernel_size = nspace_size then
                stagnatiton := stagnatiton + 1;
            else 
                stagnatiton := 0;
            fi;
            old_kernel_size := nspace_size;
            # try exit early
            if stagnatiton = stagnation_max then
                stagnatiton := 0;
                # maybe not do this computation too often..?
                # stagnation_max := stagnation_max * 2;
                O := FORMS_ReturnFormspace(true, W, g_res, g_star_inv_scaled, Lambdas, unitary, hom, Gens, d, F, g_inv_frob, n);
                if O <> false then
                    return O;
                fi;
            fi;
        od;
    od;
    
    # TODO: only check if W one dimensional or cyclic group element condition is not met, i think this can be determined from frob normal forms if we even need to do this.
    O := FORMS_ReturnFormspace(true, W, g_res, g_star_inv_scaled, Lambdas, unitary, hom, Gens, d, F, g_inv_frob, n);

    if O = false then
        # TODO: cyclic matrix conditition and example Group where this is needed if such a group exists. Can only happen for a group that does not preserve non-deg form but preserves deg form i think.
        # TODO: one could also do this instead of checking if forms are preserved however this seems slower
    fi;
    return O;
end;

# ideas for scalars : 
## non-degenerate 
# - do the things already in the forms package with trace/ analyze minimal and characteristic polynomial DONE! :)
## degenerate : 
# - do the trick where a group generator with scalar guranteed to be one gets added and compute the eigenvalues of this one these are then all possible choices of scalars, we can even do this a few times to arrow down the possibilities. 
    ## structure : 
        # 1. find a random element in the group of form ghg^{-1}h^{-1} (i.e.) it has scalar one with as few frobenius blocks (best case cyclic) in the group.
        # 2. we compute the eigenvalues of one part of the condition matrix. TODO: specify
        # 3. repeat the first two steps until scalars are sufficiently constricted
        # 4. compute formspaces of these scalars.

        ## If there is a frobenius hom wi need to do the above twice. once with the normal conjugation and once with frobenius conjugation?
## commutator group? idk if this is a good approach, ideally i want it to be able to handle every case, even if the commutator subgroup generated is trivial. the issue with this is that the resulting forms of the commutator group may not be forms of the bigger group. In this case we might need linear combinations and are still missing scalars

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
        # check if field is finite?
        Gens := GeneratorsOfGroup(G);
        # F := DefaultFieldOfMatrix(Gens[1]);
        n := NrRows(Gens[1]);
        d := Size(Gens);
        hom := fail;

        
        if unitary then
            hom := FrobeniusAutomorphism(F)^(p_exponent/2);
        fi;
        
        if d = 1 then
            Gen := Gens[1];
            Gen_adjoint := FORMS_CalculateAdjoint(Gen, unitary, hom, n, F);
            Gen_adjoint_inv_scaled := Lambda[1]*Inverse(Gen_adjoint);
            frob := FrobeniusNormalForm(Gen);
            frob_inv_star_scaled := FrobeniusNormalForm(Gen_adjoint_inv_scaled);
            frob_inv_star_base_change := Inverse(frob_inv_star_scaled[2]);

            return FORMS_CyclicGroupCase(Gen, Gen_adjoint_inv_scaled, Lambda[1], unitary, hom, frob, frob_inv_star_scaled, frob_inv_star_base_change, Inverse(frob[2]), F, n);
        fi;
        #contains  group element, scalar, (factors of minopol), Basis change to Frobenius (v, vg, vg^2, ...), Frobenius block lengths, number of iterations to compute
        g_res := FORMS_FindCyclicGroupElementAndScalars(Gens, Lambdas, n);
        g_inv_frob := Inverse(g_res[4]);
        g_star_inv_unscaled := TransposedMat(Inverse(g_res[1]));
        if unitary then
            g_star_inv_unscaled := g_star_inv_unscaled^hom;
        fi;
        g_star_inv_scaled_frob := FrobeniusNormalForm(g_star_inv_unscaled * g_res[2]);
        frob_base_inv_star := Inverse(g_star_inv_scaled_frob[2]);

        return FORMS_FormspaceInternal(Gens, Lambdas, unitary, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n);
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
        # Todo check if F finite?
        Gens := GeneratorsOfGroup(G);
        n := NrRows(Gens[1]);
        d := Size(Gens);
        hom := fail;

        if d = 1 then
            Gen := Gens[1];
            Gen_adjoint := FORMS_CalculateAdjoint(Gen, false, fail, n, F);
            Gen_adjoint_inv_scaled := One(F)*Inverse(Gen_adjoint); # scaling happens here!!
            frob := FrobeniusNormalForm(Gen);
            frob_inv_star_scaled := FrobeniusNormalForm(Gen_adjoint_inv_scaled);
            frob_inv_star_base_change := Inverse(frob_inv_star_scaled[2]);
            g_inv_frob := Inverse(frob[2]);

            Add(Out, FORMS_CyclicGroupCase(Gen, Gen_adjoint_inv_scaled, One(F), false, fail, frob, frob_inv_star_scaled, frob_inv_star_base_change, g_inv_frob, F, n));

            if p_exponent mod 2 = 0 then
                hom := FrobeniusAutomorphism(F)^(p_exponent/2);
                Gen_adjoint_inv_scaled := Gen_adjoint_inv_scaled^hom;
                frob_inv_star_base_change := frob_inv_star_base_change^hom;
                frob_inv_star_scaled[2] := frob_inv_star_scaled[2]^hom;

                Add(Out, FORMS_CyclicGroupCase(Gen, Gen_adjoint_inv_scaled, One(F), true, hom, frob, frob_inv_star_scaled, frob_inv_star_base_change, g_inv_frob, F, n));
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
        g_res := FORMS_FindCyclicGroupElementAndScalars(Gens, Lambdas, n);
        g_inv_frob := Inverse(g_res[4]);
        g_star_inv_unscaled := TransposedMat(Inverse(g_res[1]));
        g_star_inv_scaled_frob := FrobeniusNormalForm(g_star_inv_unscaled); # maybe this causes a bug with scalars
        frob_base_inv_star := Inverse(g_star_inv_scaled_frob[2]);

        Add(Out, FORMS_FormspaceInternal(Gens, Lambdas, false, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n));
        if p_exponent mod 2 = 0 then
            # is ^hom actually cheaper than computing the frobenius normal form?? investigate! certainly makes the code ugly... oh well
            hom := FrobeniusAutomorphism(F)^(p_exponent/2);
            g_star_inv_unscaled := g_star_inv_unscaled^hom;
            g_star_inv_scaled_frob[2] := g_star_inv_scaled_frob[2]^hom;
            frob_base_inv_star := frob_base_inv_star^hom;
            Add(Out, FORMS_FormspaceInternal(Gens, Lambdas, true, hom, g_res, g_inv_frob, g_star_inv_unscaled * g_res[2], g_star_inv_scaled_frob, frob_base_inv_star, d, F, n));
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
        return FORMS_FilterBilinearForms(Forms, F, n);
    else
        p_exponent := DegreeOverPrimeField(F);
        if p_exponent mod 2 <> 0 then 
            Error("The given Field ", F, " must admit a field automorphism of order two for unitary = true\n");
            return [];
        fi;
        hom := FrobeniusAutomorphism(F)^(p_exponent/2);
        return FORMS_FilterUnitaryForms(Forms, F, n, hom);
    fi;
end);


# computes non degenerate forms and scalars, TODO: random_non_deg should specify whether the function should use random methods to make sure that the returned forms are non deg, this is a todo and not necessary for irreducible groups.
FORMS_PreservedNonDegFormsWithScalarsOp := function(G, random_non_deg, sesquilinear)
    local F, Gens, n, p_exponent, hom, newGens, forms_bil, forms_herm, filtered, space, forms_out, filt_bil, filt_herm, Form, quad, x, biglist;
    F := DefaultFieldOfMatrixGroup(G);
    Gens := GeneratorsOfGroup(G);
    n := NrRows(Gens[1]);
    p_exponent := DegreeOverPrimeField(F);
    forms_bil := [];
    forms_herm := [];

    hom := FrobeniusAutomorphism(F)^0; # better way to obtain identity mapping?
    newGens := ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(G, hom);
    filt_bil := [[], []];
    if newGens <> false then
        # this should be improved, many computations from PreservedFormspace could probably be recycled for example the cyclic group element, and so on... Does not seem to be such a big issue, this is only inefficient for groups that actually have multiple possible scalars which seems to be very rare. Infact TODO: test case with (irreducible) group which preserves forms modulo different sets of scalars. Does such a group even exist?
        biglist := Cartesian(newGens[2]);
        # Print(biglist);
        for x in biglist do
            forms_bil := PreservedFormspace(Group(newGens[1]), x, false);
            filtered := FilterFormspace(forms_bil, F, false);
            if forms_bil = [] then
                # do nothing .. :)
            elif Characteristic(F) = 2 then
                for Form in filtered do
                    if not sesquilinear then
                        quad := ClassicalForms_QuadraticForm2(F, Form, newGens[1], x);
                        if quad = false then
                            Add(filt_bil[2], Form);
                        else
                            Add(filt_bil[1], quad);
                        fi;
                    else
                        Add(filt_bil[2], Form);
                    fi;
                od;
            else 
                filt_bil[2] := Concatenation(filt_bil[2], filtered[2]);
                for Form in filtered[1] do
                    if not sesquilinear then
                        quad := ClassicalForms_QuadraticForm(F, Form);
                        Add(filt_bil[1], quad);
                    else
                        Add(filt_bil[2], Form);
                    fi;
                od;
            fi;
        od;
    fi;
   
    filt_herm := [];
    if p_exponent mod 2 = 0 then
        hom := FrobeniusAutomorphism(F)^(p_exponent/2);
        newGens := ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(G, hom);
        if newGens <> false then
            biglist := Cartesian(newGens[2]);
            for x in biglist do
                Add(forms_herm, PreservedFormspace(Group(newGens[1]), x, true));
            od;
        fi;
        for space in forms_herm do
            filt_herm := Concatenation(filt_herm, FilterFormspace(space, F, true));
        od;
    fi;
    
    ## process the forms in forms_bil, forms_herm
    forms_out := [];
    for Form in filt_bil[1] do
        Add(forms_out, QuadraticFormByMatrix(Form, F));             
    od;
    for Form in filt_bil[2] do
        Add(forms_out, BilinearFormByMatrix(Form, F));
    od;
    for Form in filt_herm do
        Add(forms_out, HermitianFormByMatrix(Form, F));
    od;
    return forms_out;
end;


#! @Arguments G
#! @Returns List of preserved Forms by G
#! @Description
#! <A>G</A> is a finitely generated matrix group over a finite field $K$. 
#! This function computes scalars $\lambda$ and forms $\Phi$ such that <A>G</A> preserves $\Phi$ modulo $\lambda$. The function does only find scalars that belong to non-degenerate forms but may return degenerate forms if the group action is reducible. 
InstallMethod(PreservedFormsWithScalars, "for matrix groups", [IsMatrixGroup], 
function(G)
    return FORMS_PreservedNonDegFormsWithScalarsOp(G, true, false);
end);

#! @Arguments G
#! @Returns List of preserved sesquilinear Forms by G
#! @Description 
#! Very similar to <A>PreservedFormsWithScalars</A> where <A>G</A> is a finitely generated matrix group over a finite field $K$. 
#! This function computes scalars $\lambda$ and forms $\Phi$ such that <A>G</A> preserves $\Phi$ modulo $\lambda$. The function does only find scalars that belong to non-degenerate forms but may return degenerate forms if the group action is reducible.
InstallMethod(PreservedSesquilinearFormsWithScalars, "for matrix groups", [IsMatrixGroup], 
function(G)
    return FORMS_PreservedNonDegFormsWithScalarsOp(G, true, true);
end);
