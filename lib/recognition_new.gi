#############################################################################
##
##  recognition_new.gi              'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##                                                              Frank Celler
##  Copyright 2024, Vrije Universiteit Brussel
##  Copyright 2024, The University of Western Austalia
##  Copyright (C) 2024,  Lehrstuhl D fuer Mathematik, RWTH Aachen, Germany
##
##  This file contains functions to compute sesquilinear forms invariant under
##  matrix groups. It is work in progress. This file contains the additions
##  and corrected version of some functions in recognition.gi
##
##  *** Bamberg and De Beule are very grateful to Frank Celler for
##  generously providing the bulk of this code.  ***
##

#############################################################################
##
#F  ClassicalForms_PossibleScalarsSesquilinear( <field>, <mat>, <frob> )
# This function returns [i0,a] such that if <mat> preserves a sesquilinear form B
# modulo scalars then the scalar lambda for which mat*B*(mat^T)^frob = lambda*B
# satisfies lambda^i0 = a.
# Very important: this function is meant to be called *only* for
# sesquilinear forms, that means, the "bar map" is the argument <frob>,
# which is the identity if we are looking for a bilinear form, and
# the involutory field automorphism in case of a hermitian form.
##
InstallGlobalFunction( ClassicalForms_PossibleScalarsSesquilinear,
 function( F, M, frob )
    local cpol, d, c, z, I, t, a, l, q, i0, Minv;

    # compute the characteristic polynomial of <M>
    cpol := CharacteristicPolynomial(M);

    # get the position of the non-zero coefficients
    d := Degree(cpol);
    c := CoefficientsOfUnivariatePolynomial(cpol);
    z := Zero(F);
    q := Size(F);
    I := Filtered( [ 0 .. d ],  x -> c[x+1] <> z );

    #Lemma: Trace(M) <> z implies that Trace(M^*) <> z.
    # Use TraceMatrix to support Matrix Objs
    Minv := M^-1;
    if TraceMatrix(M) = z and TraceMatrix(Minv) <> z or
       TraceMatrix(M) <> z and TraceMatrix(Minv) = z then
        return false;
    fi;
    
    # make sure that <d> and <d>-i both occur, i.e. check that the support of cpol is symmetric
    
    if ForAny( I, x -> not (d-x) in I )  then
        return false;
    fi;

    # we need gcd one in order to get alpha exactly (ignoring +-)
    Add(I,q-1);
    t := GcdRepresentation(I);
    i0 := I*t;

    #Lemma 1: if g preserves a form modulo lambda then c[d-I](c_0^frob) = (c_i^frob) lambda^i.
    #cpol = c_0 +c_1 X + \ldots + c_d X^d
    #l will be the list of the lambda^i's.
    l:=List([1..Length(I)-1], x ->((c[1]^frob)*c[d-I[x]+1]/(c[I[x]+1]^frob))); #here is the line with identity field automorphism!

    a:=Product([1..Length(I)-1], x->l[x]^t[x]);
    # Now the scalar $\lambda$ satisfies $\lambda^{i_0}=a$

    # check: $\forall_i: c_{d-i}c_0^frob=c_i^frob\lambda^i
    # if any of the conditions from Lemma 1 fail, g does not preserve the form modulo this lambda.
    if ForAny([1..Length(I)-1],x->(l[x]<>a^QuoInt(I[x],i0))) then
        Info( InfoForms, 1, "Characteristic polynomial proves that there is no scalar\n" );
      return false;
    fi;

    return [i0,a];
  end );


#############################################################################
##
#F  ClassicalForms_GeneratorsWithBetterScalarsSesquilinear( grp, frob )
# Given a group <grp>, this function tries to compute a list of generators
# for <grp> that preserve a form modulo the returned scalars.
# Note that the user determines whether looking for a preserved bilinear, respectively
# hermitian form, by choosing frob to be the identiy, respectively the involutory
# field automorhipsm. This is actually just a massage of the generating set of the group.
##
InstallGlobalFunction( ClassicalForms_GeneratorsWithBetterScalarsSesquilinear,
  function( grp, frob )
    local tries, gens, field, m1, a1, new, i, scalars, improvegenerator, res, newgens, champion, len, qq, q;

    #field := FieldOfMatrixGroup(grp); #this causes a problem if one has a maximal subgroup of e.g. SU(3,7^2). There are examples of which de FieldOf is GF(7), while DefaultFieldOF is GF(7^2). Then this causes problems.
    field := DefaultFieldOfMatrixGroup(grp);
    qq := Size(field);
    if not IsOne(frob) then
        q := Sqrt(qq);
    else
        q := 1;
    fi;

    # the next function returns a matrix with a list of possible scalars for this matrix.

    improvegenerator := function(m1,i,count,len)
        local a1, s, j, k, scalars, r, w;
        
        a1 := ClassicalForms_PossibleScalarsSesquilinear(field,m1,frob);
        #Recall that the scalars satisfy the equation lambda^a[1] = a[2]
        if a1 = false then
            return false; #The group does not preserve a form modulo scalars
        fi;
                
        if IsList(a1) then
            if a1[1] = 1 then # the matrix m1 has scalar a1[2]
                if LogFFE(a1[2],PrimitiveRoot(field)) mod (q+1) = 0 then
                    return [m1/NthRoot(field,a1[2],q+1),[One(field)]];
                fi;
                return [m1,[a1[2]]]; #originally, the three lines above this return were not there. Those three lines make sure scalar becomes 1 if possible (basically if there is a q+1-st root).
            else
                # All solutions of lambda^a1[1] = a1[2] are r*w^k, where r is one solution and w is a primitive a1[1]-th root of unity (a1[1] divides Size(field) - 1).
                # Previously this was AsList(Group(NthRoot(field,a1[2],a1[1]))), the group generated by a single root. That is not the solution set and lost scalars, e.g. for lambda^3 = 1 it only gave [1], so forms preserved modulo other scalars were not found. An example of this is the group G22 in tst/intersting_groups.g
                r := NthRoot(field,a1[2],a1[1]);
                if r = fail then
                    return false; # no solution, so the group does not preserve a form modulo scalars
                fi;
                w := PrimitiveRoot(field)^((Size(field) - 1) / a1[1]);
                scalars := List([0 .. a1[1] - 1], k -> r * w^k); # add all possible scalars for m1
                if count = 0 then
                    return champion; # add all possible scalars for m1
                elif Length(scalars) < len then
                    champion := [m1,scalars];
                    len := Length(scalars);
                fi;
                k := Random(Difference([1..Length(gens)],[i])); #m1 could not be improved, so try to change m1 by multiplying with another generator.
                m1 := m1*newgens[k]; # this should be newgens, such that the improved generators still generate the group instead of a subgroup.
                return improvegenerator(m1,i,count-1,len);
            fi;
        fi;
    end;


    # start with 2 random elements,  at most 10 tries
    
    # this logic for cyclic groups is not helpfull to the formspace function, as it will return a group with more than one generator, this is bad for the preserved formspace funcitons, as they have extra logic for dealing with groups with one generator.
    # tries := 0;
    gens  := ShallowCopy(GeneratorsOfGroup(grp));
    # if Length(gens) = 1 then
    #     Add(gens,PseudoRandom(grp));
    # fi;
    #We will randomize the generating set in the hope that we obtain generators preserving fewer scalars.
     
    scalars := [];
 
    newgens := ShallowCopy(gens);
        
    for i in [1..Length(gens)] do
        champion := [gens[i],AsList(Group(PrimitiveRoot(field)))];
        len := Length(champion[2]);
        res := improvegenerator(gens[i],i,10,len);
        if res = false then
            return false; #The group does not preserve a bilinear form modulo scalars
        fi;
        newgens[i] := res[1];
        scalars[i] := res[2];
    od;

    return [newgens,scalars];

  end );

#This function tests wheter the group grp preserves a form given by the matrix form modulo scalars if so it returns the scalar for each generator. <form> must be a form object.
#############################################################################
##
#O  ScalarsOfPreservedForm( <grp>, <form> )
##  given <grp> and <form>, check whether <grp> preserves <forms> modulo scalars
# and returns scalars if so, false otherwise.
##
InstallGlobalFunction(ScalarsOfPreservedForm,
    function(grp,form)

    local i, r, c, gens, scalars, a, m, gram, zero, z, gf, frob;
    gens := GeneratorsOfGroup(grp);
    gf := DefaultFieldOfMatrixGroup(grp);
    scalars := [];
    gram := GramMatrix(form);
    z := 0*gram[1][1];
    zero := List([1..Length(gram)],i->z);
    r := First([1..Length(gram)],x->gram[x] <> zero);
    if r = fail then
        # the trivial form is preserved modulo any scalars
        return List(gens, g -> One(gf));
    fi;
    c := First([1..Length(gram)],x->gram[r][x] <> z);
    for i in [1..Length(gens)] do
        if IsQuadraticForm(form) then
            m := gens[i] * gram * TransposedMat(gens[i]);
            m := Forms_RESET(m,Length(gram),gf);
        else
            frob := CompanionAutomorphism(form);
            if IsOne(frob) then
                m := gens[i] * gram * TransposedMat(gens[i]);
            else
                m := gens[i] * gram * TransposedMat(gens[i])^frob;
            fi;
        fi;
        a := m[r][c] * gram[r][c]^-1;
        if m <> a * gram then
            return false;
        fi;
        scalars[i] := a;
    od;
    return scalars;
end );

#############################################################################
##
#O  TestPreservedSesquilinearForms( <grp>, <form> )
##  given <grp> and a list of forms, check whether <grp> preserves <forms> modulo scalars.
##  This function is meant to return true or the position in the list for which the
##  test fails.
##
InstallGlobalFunction(TestPreservedSesquilinearForms,
    function(grp,list)

    local form, mat, result, i, fails, gens, aut;
    gens := GeneratorsOfGroup(grp);
    fails := [];
    for i in [1..Length(list)] do
        form := list[i];
        mat := GramMatrix(form);
        aut := CompanionAutomorphism(form);
        result := List(gens,x->_IsEqualModScalars(x*mat*(TransposedMat(x)^aut),mat));
        if false in result then
            Add(fails,i);
        fi;
    od;
    if fails = [] then
        return true;
    else
        return fails;
    fi;
end );
