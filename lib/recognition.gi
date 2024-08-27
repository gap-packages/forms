#############################################################################
##
##  recognition.gi              'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##                                                              Frank Celler
##  Copyright 2024, Vrije Universiteit Brussel
##  Copyright 2024, The University of Western Austalia
##  Copyright (C) 2024,  Lehrstuhl D fuer Mathematik, RWTH Aachen, Germany
##
##  This file contains functions to compute sesquilinear forms invariant under
##  matrix groups.
##
##  *** Bamberg and De Beule are very grateful to Frank Celler for
##  generously providing the bulk of this code.  ***
##

#############################################################################
##
#F  ClassicalForms_PossibleScalarsSesquilinear( <field>, <mat>, <frob> )
# Very important: this function is meant to be called *only* for
# sesquilinear forms, that means, the "bar map" is the argument <frob>,
# which is the identity if we are looking for a bilinear form, and
# the involutory field automorphism in case of a hermitian form.
# This function returns [i0,a] such that if <mat> preserves a sesquilinear form B
# modulo scalars then the scalar lambda for which mat*B*(mat^T)^frob = lambda*B
# satisfies lambda^i0 = a
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

    Minv := M^-1;
    if Trace(M) = z and Trace(Minv) <> z or
       Trace(M) <> z and Trace(Minv) = z then
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
# Note that the user determines whether looking for a preserved bilinear, respectively
# hermitian form, by choosing frob to be the identiy, respectively the involutory
# field automorhipsm. This is actually just a massage of the generating set of the group.
#
##
InstallGlobalFunction( ClassicalForms_GeneratorsWithBetterScalarsSesquilinear,
  function( grp, frob )
    local tries, gens, field, m1, a1, new, i, scalars, root, improvegenerator, res, newgens, champion, len, qq, q;

    # the aim of this function is to replace the matrix m1 by a
    # matrix that has as few solutions to the scalar equation
    # lambda^a1[1] = a1[2] as possible. It checks first if a1[1] = 1,
    # since then lambda is determined. Next we check if a1[2] is a
    # square. And then we replace m1 by a matrix that leaves the
    # bilinear form invariant modula the scalar 1. If none of these
    # are possible, we try to replace m1 by a matrix that has fewer
    # solutions to the scalar equation.

    #field := FieldOfMatrixGroup(grp); #this causes a problem if one has a maximal subgroup of e.g. SU(3,7^2). There are examples of which de FieldOf is GF(7), while DefaultFieldOF is GF(7^2). Then this causes problems.
    field := DefaultFieldOfMatrixGroup(grp);
    qq := Size(field);
    if not IsOne(frob) then
        q := Sqrt(qq);
    else
        q := 1;
    fi;

    # the next function returns a matrix with a list of possible scalars for this matrix.
    # if it is possible to change the matrix to multiple that preserves up to scalar one, this is achieved.

    improvegenerator := function(m1,i,count,len)
        local a1, s, j, k, scalars, root; #I think root should be declare locally here!
        
        a1 := ClassicalForms_PossibleScalarsSesquilinear(field,m1,frob);
        #Recall that the scalars satisfy the equation lambda^a[1] = a[2]
        if a1 = false then
            return false; #The group does not preserve a form modulo scalars
        fi;
                
        if IsList(a1) then
            root := NthRoot(field,a1[2],a1[1]);
            if a1[1] = 1 then # the matrix m1 has scalar a1[2]
                #if a1[2] has a square root, we can replace m1 with m1*sqrt{a1[2]};
                if LogFFE(a1[2],PrimitiveElement(field)) mod 2 = 0 then
                    return [m1/NthRoot(field,a1[2],2),[One(field)]];
                fi;
                return [m1,[a1[2]]]; #originally, the three lines above this return were not there. Those three lines make sure scalar becomes 1 if possible (basicaly if there is a sqrt).
            elif LogFFE(root,PrimitiveElement(field)) mod (q+1) = 0 then
                return [m1/NthRoot(field,root,q+1),[One(field)]]; #either frob = id, then q+1 = 2, or frob is not trivial, then we take q+1-st root. In both cases, modify m1 to a matrix that has scalar one.
            else
                scalars := AsList(Group(NthRoot(field,a1[2],a1[1]))); # add all possible scalars for m1
                if count = 0 then
                    return champion; # add all possible scalars for m1
                elif Length(scalars) < len then
                    champion := [m1,scalars];
                    len := Length(scalars);
                fi;
                k := Random(Difference([1..Length(gens)],[i])); #m1 could not be improved, so try to change m1 by multiplying with another generator.
                m1 := m1*gens[k];
                return improvegenerator(m1,i,count-1,len);
            fi;
        fi;
    end;

    # start with 2 random elements,  at most 10 tries
    tries := 0;
    gens  := ShallowCopy(GeneratorsOfGroup(grp));
    if Length(gens) = 1 then
        Add(gens,PseudoRandom(grp));
    fi;
    #We will randomize the generating set in the hope that we obtain generators preserving fewer
    #scalars.
     
    scalars := [];
 
    newgens := ShallowCopy(gens);
        
    for i in [1..Length(gens)] do
        champion := [gens[i],AsList(Group(PrimitiveElement(field)))];
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


#############################################################################
##
#F  ClassicalForms_Signum2( <field>, <form>, <quad> )
# This function computes a base change, seemingly sufficient to see from the
# changed form what its type is, all for characteristic 2.
##
InstallGlobalFunction( ClassicalForms_Signum2,
  function( field, form, quad )
    local dim, base, avoid, i, d, j, c, k, x, sgn, pol;

    # compute a new basis,  such that the symmetric form is standard
    #Error("here in ClassicalForms_Signum2");
    base :=OneOp(form);
    form:=MutableCopyMat(form);
    avoid := [];
    dim := NrRows(form);
    for i  in [ 1 .. dim-1 ]  do

        # find first non zero entry
        d:=1;
        while d in avoid or IsZero(form[i,d])  do
            d:=d+1;
        od;
        Add(avoid, d);

        # clear all other entries in this row & column
        for j  in [d+1..dim]  do
            c := form[i,j]/form[i,d];
            if c <> Zero(field)  then
                for k  in [ i .. dim ]  do
                    form[k,j] := form[k,j] - c*form[k,d];
                od;
                AddRowVector(form[j],form[d],c);
                AddRowVector(base[j],base[d],c);
            fi;
        od;
    od;

    # reshuffle base
    c := [];
    j := [];
    for i  in [ 1 .. dim ]  do
        if not i in j  then
            k := form[i,avoid[i]];
            Add( c, base[i]/k );
            Add( c, base[avoid[i]] );
            Add( j, avoid[i] );
        fi;
    od;
    base := c;

    # and try to fix the quadratic form (this is not really necessary)
    x   := X(field);
    sgn := 1;
    for i  in [ 1, 3 .. dim-1 ]  do
        c := base[i] * quad * base[i];
        if IsZero(c)  then
            c := base[i+1] * quad * base[i+1];
            if not IsZero(c) then
                AddRowVector(base[i+1],base[i],-c);
            fi;
        else
            j := base[i+1] * quad * base[i+1];
            if IsZero(j)  then
                AddRowVector(base[i],base[i+1],-c);
            else
                pol := Factors(x^2 + x/j + c/j);
                if Length(pol) = 2  then
                    pol:=List(pol,x->-CoefficientsOfUnivariatePolynomial(x)[1]);
                else
                    sgn := -sgn;
                fi;
            fi;
        fi;
    od;
    return sgn;
  end );


#############################################################################
##
#F  ClassicalForms_Signum( <field>, <form>, <quad> )
##
InstallGlobalFunction( ClassicalForms_Signum,
  ##
  ## This could be replaced by operations in "Forms".
  ## We are, however, not convinced any more, as Determinant might be faster
  ## than computing the base changes.
  ##
  function( field, form, quad )
    local sgn, det, sqr;

    # if dimension is odd,  the signum must be 0
    if NrRows(form) mod 2 = 1  then
        return [ 0 ];

    # hard case: characteristic is 2
    elif Characteristic(field) = 2  then
        Error( "characteristic must be odd" );
    fi;

    # easy case
    det := DeterminantMat(form);
    sqr := LogFFE( det, PrimitiveRoot(field)) mod 2 = 0;
    if (NrRows(form)*(Size(field)-1)/4) mod 2 = 0  then
        if sqr  then
            sgn := +1;
        else
            sgn := -1;
        fi;
    else
        if sqr  then
            sgn := -1;
        else
            sgn := +1;
        fi;
    fi;

    return [ sgn, sqr ];
end );


#############################################################################
##
#F  ClassicalForms_QuadraticForm2( <field>, <form>, <gens>, <scalars> )
##  <form> is a given bilinear form, preserved by <gens> modulo <scalars>.
##  This function computes a quadratic form which has <form> as polar form,
##  and the quadratic form is preserved by gens (if that is all possible).
##  Characteristic is two. Note that since the polar form <form> is given,
##  the quadratic form is already determined above the diagonal. So only
##  the elements on the diagonal have to be computed.
##
##  Because of the optimizations in improvegenerators, we may assume that generators
##  indeed preserve the sesquilinear form up to scalar one, since in characteristic
##  2 we can always achieve this.
##
InstallGlobalFunction( ClassicalForms_QuadraticForm2,
 function( field, form, gens, scalars )
    local dim, H, i, j, e, b, y, x, r, l;

    # raise an error if char is not two
    if Characteristic(field) <> 2  then
        Error( "characteristic must be two" );
    fi;

    # construct the upper half of the quadratic form
    H := ZeroOp(form);
    dim := NrRows(form);
    for i  in [ 1 .. dim ]  do
        for j  in [ i+1 .. dim ]  do
            H[i,j] := form[i,j];
        od;
    od;

    # store the linear equations in <e>
    e := [];

    # loop over all generators
    b := [];
    for y  in [ 1 .. Length(gens) ]  do

        # remove scalars
        x := gens[y]*scalars[y]^-1;

        # first the right hand size
        r := x*H*TransposedMat(x)+H;

        # check <r>
        # here we use that all scalars are one
        # observe that x*H*TransposedMat(x)+H will be symmetric now
        # because H is the upper triangle part of form, hence whatever
        # change is made to H by x*H*TransposedMat(x), if you add H, the
        # result must be symmetric.
        # In fact, we now know mathematically, that scalars should always be one
        # if the given form indeed comes from a quadratic form.
        for i  in [ 1 .. dim ]  do
            for j  in [ i+1 .. dim ]  do
                if not IsZero(r[i,j]+r[j,i]) then
                    Print("returning false at symmetry check\n");
                    return false;
                fi;
            od;
        od;

        # and now the diagonals
        for i  in [ 1 .. dim  ]  do
            l := [];
            for j  in [ 1 .. dim ]  do
                l[j] := x[i,j]^2;
            od;
            l[i] := l[i]+1;
            Add( b, r[i,i] );
            Add( e, l );
        od;
    od;

    # and return a solution
    Print("Rank of e:",RankMat(e),"\n");
    e := SolutionMat( TransposedMat(e), b );
    if e <> fail  then
        for i  in [ 1 .. dim ]  do
            H[i,i] := e[i];
        od;
        return ImmutableMatrix(field,H);
    else
        Print("returning false at no solution found\n");
        return false;
    fi;
  end );


#############################################################################
##
#F  ClassicalForms_QuadraticForm( <field>, <form> )
##
InstallGlobalFunction( ClassicalForms_QuadraticForm,
  function( field, form )
    local H, i, j;

    # special case if <p> = 2
    if Characteristic(field) = 2  then
        Error( "characteristic must be odd" );
    fi;

    # use upper half
    H := ZeroOp(form);
    for i  in [ 1 .. NrRows(form) ]  do
        H[i,i] := form[i,i]/2;
        for j  in [ i+1 .. NrRows(form) ]  do
            H[i,j] := form[i,j];
        od;
    od;
    return H;
  end );

InstallGlobalFunction( DualFrobeniusGModule,
    function(module,frob)
    #make sure frob is involution.
    local   mats,  dmats;

    if SMTX.IsZeroGens(module) then
        return GModuleByMats([],module.dimension,SMTX.Field(module));
    else
        mats  := MTX.Generators(module);
        dmats := List(mats,i->TransposedMat(i^frob)^-1);
        return GModuleByMats(dmats,MTX.Field(module));
    fi;
end );

#############################################################################
##
#F  ClassicalForms_InvariantForms( <gens_scalars> )
#   We generate the module by generators stored in gens_scalars and we know
#   that this module is the module of the original group. Then we construct
#   all possible dual modules according to the list stored in gens_scalars.
#   Then we check which ones yield bilinear forms.
##
InstallGlobalFunction( ClassicalForms_InvariantForms,
    function( gens_scalars, frob )
    local   hom,  scalars,  form,  iform,  identity,  field,  root,
            q,  i,  m,  a,  quad,  sgn, gmodule, forms, biglist, x,
            dmodule, gens, scale_gens, module, output, isform, pair, transposedform,
            check_scalar_matrix, lambda, mu;

    # <dmodule> acts absolutely irreducible without scalars
    gens := gens_scalars[1];
    field := FieldOfMatrixGroup(Group(gens));
    gmodule := GModuleByMats(gens,field);
    forms := [];
    biglist := Cartesian(gens_scalars[2]);
    # we create all possible dual modules that correspond to the scalars we computed.
    for x in biglist do
        scale_gens := List([1..Length(gens)],i->gens[i] / x[i]);
        #Don't compute the dual module multiple times just because the scalars change... to be improved!
        if IsOne(frob) then
            dmodule := DualGModule(GModuleByMats(scale_gens,field));
        else
            dmodule := DualFrobeniusGModule(GModuleByMats(scale_gens,field),frob);
        fi;
        hom := MTX.Homomorphisms( gmodule, dmodule);
        Info( InfoForms, 1, "found homomorphism between V and V^*\n" );
        Append(forms,List(hom,i->[i,x]));
    od;
    #TO DO: get rid of these checks in the final version.
    # make sure that the forms commute with the generators of <module>
    # forms is a list now. We should run the following checks for each of the elements in forms.
    
    output := [];

    check_scalar_matrix := function(T,F)
        local i,lambda;
        i := PositionNonZero(T[1]);
        if i <> PositionNonZero(F[1]) then
            return false;
        fi;
        lambda := T[1][i] / F[1][i];
        if T <> lambda*F then
            return false;
        else
            return lambda;
        fi;
    end;

    for pair in forms do
        form := pair[1];
        scalars  := pair[2];
        transposedform := TransposedMat(form);

        # check the type of form, replace this afterwards with default forms constructors. Then also join the two for loops.
            
        if not IsOne(frob) then
        q := Sqrt(Size(field));
        #check if the form is hermitian. Note that the condition should be TransposedMat(form) = lambda form^frob
        # for some scalar lambda. If so, then mu^(q-1) = \lambda^q. Then mu = RootFFE(GF(q^2),lambda^q,q-1);
            #if TransposedMat(form) <> form^frob then
            lambda := check_scalar_matrix(transposedform,form^frob);
            if lambda <> false then # the form is actually hermitian, we need to change its gram matrix.
                if lambda <> One(field) then
                    mu := RootFFE(field,lambda,q-1);
                    Info( InfoForms, 1,"unknown form made into hermitian\n" );
                else mu := lambda;
                fi;
                #Error("here is now a mistake");
                Add(output, [ "unitary", mu*form, scalars ]);
            else
                #Error("make sure we see what happens");
                Add(output, [ "unknown", "hermitian", form, scalars ]);
            fi;
        else
        #now we know the form is not hermitian.
 
             if TransposedMat(form) = -form  then
                Info(InfoForms, 1, "form is symplectic\n" );
                if Characteristic(field) = 2  then
                    quad := ClassicalForms_QuadraticForm2(field, form, gens, scalars );
                    if quad = false  then
                        Add(output, [ "symplectic", form, scalars ]);
                    elif MTX.Dimension(gmodule) mod 2 = 1  then
                        Error( "no quadratic form but odd dimension" );
                    elif ClassicalForms_Signum2( field, form, quad ) = -1  then
                        Add(output, [ "orthogonalminus", form, scalars, quad ]);
                    else
                        Add(output, [ "orthogonalplus", form, scalars, quad ]);
                    fi;
                else
                    Add(output,[ "symplectic", form, scalars ]);
                fi;
            elif TransposedMat(form) = form  then
                Info(InfoForms, 1, "form is symmetric\n" );
                quad := ClassicalForms_QuadraticForm( field, form );
                if MTX.Dimension(gmodule) mod 2 = 1  then
                    Add(output, [ "orthogonalcircle", form, scalars, quad ]);
                else
                    sgn := ClassicalForms_Signum( field, form, quad );
                    if sgn[1] = -1  then
                        Add(output, [ "orthogonalminus", form, scalars, quad, sgn[2] ]);
                    else
                        Add(output, [ "orthogonalplus", form, scalars, quad, sgn[2] ]);
                    fi;
                fi;
            else
                Info( InfoForms, 1,"unknown form\n" );
                Add(output, [ "unknown", "dual", form, scalars ]);
            fi;
        fi;
    od;
    return output;
end );

#############################################################################
##
#O  ScalarOfSimilarity( <grp> )
##
##    returns a scalar. Is this somewhere used?
##
## To do: see if the check function can get this name and/or be merged with this function.
InstallMethod( ScalarOfSimilarity, [IsMatrix, IsSesquilinearForm],
  function( g, form )

    ## Recall that a similarity of a form f on V, is a linear transformation g
    ## of V where there exists some nonzero scalar c such that for all v,w in V
    ##         f(u^g,v^g) = c f(u,v).
    ## This operation finds for a particular matrix g, giving rise to a similarity
    ## of "form", the scalar c.

    local gram, m, pos, scalar;

    ## check that g and form are compatible in dimension and there fields are OK
    gram := GramMatrix( form );
    if NrRows(g) <> NrRows( gram ) then
       Error("dimensions are incompatible.");
    fi;

    if not ForAll(Flat(g), t -> t in form!.basefield) then
       Error("fields are incompatible");
    fi;

    ## check that g is invertible
    if IsZero(Determinant(g)) then
       Error("g must be invertible");
    fi;

    ## now check to see if g is a similarity, and then
    ## determine the scalar

    m := EvaluateForm(form, g, g);
    pos := PositionNonZero( m[1] );
    scalar := m[1,pos] / gram[1,pos];
    return scalar;
  end );

#############################################################################
##
#O  PreservedFormsOp( <grp> )
##  return a record containing information on preserved forms. This operation
##  is not intended for the user. Its output will be processed in different
##  operations.
##
InstallMethod( PreservedFormsOp, [ IsMatrixGroup ],
  function( grp )
    local   forms, field, i, g, qq, c, module,  invariantforms,
            gens_scalars, fmodule, form, y, newform, newforms, bilinearforms, frob, hermitianforms;
    
    field := DefaultFieldOfMatrixGroup(grp);

    # remember to rename Dual to Bilinaer and Frobenius to Hermitian.

    forms := rec();
    forms.field := field;
    forms.invariantforms := [];
    forms.maybeDual      := false;
    forms.maybeFrobenius := false;

    frob := FrobeniusAutomorphism(field);
    
    # set up the module and other information
    module := GModuleByMats(GeneratorsOfGroup(grp), field);

    # set the possibilities
    forms.maybeDual      := true;
    forms.maybeFrobenius := IsEvenInt(DegreeOverPrimeField(field));

    if forms.maybeFrobenius  then
        qq := Characteristic(field)^(DegreeOverPrimeField(field)/2);
        frob := frob^(DegreeOverPrimeField(field)/2);
    fi;

    # <grp> must act irreducibly
    #Is this really necessary? Can we not simply delete it?
    if not MTX.IsIrreducible(module)  then
        #Error("Currently the use of MeatAxe requires the module to be absolutely irreducible");
        Info( InfoForms, 1,  "group is not irreducible and therefore it does not preserve non-degenerate forms\n" );
        return [];
    fi;

    # <grp> must act absolutely irreducibly
    #Is this really necessary? Can we not simply delete it?
    if not MTX.IsAbsolutelyIrreducible(module)  then
        #Error("Currently the use of MeatAxe requires the module to be absolutely irreducible");
        Info( InfoForms, 1,  "grp not absolutely irreducible\n" );
        #return [];
    fi;

    # try to find generators without scalars
    if forms.maybeDual then
        gens_scalars := ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(grp,frob^0);
        if gens_scalars = false then #happens if the group does not preserve a bilinear form modulo scalars
            forms.maybeDual := false;
        fi;
    fi;
    
    # now try to find an invariant form
    if forms.maybeDual  then
        bilinearforms := ClassicalForms_InvariantForms(gens_scalars,frob^0);
        if bilinearforms <> false  then
            Append( forms.invariantforms, bilinearforms );
        else
            forms.maybeDual := false;
        fi;
    fi;

    if forms.maybeFrobenius  then
        gens_scalars := ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(grp,frob);
        if gens_scalars = false  then
            forms.maybeFrobenius := false;
        fi;
    fi;

    if forms.maybeFrobenius  then
        hermitianforms := ClassicalForms_InvariantForms(gens_scalars,frob);
        if hermitianforms <> false  then
            Append( forms.invariantforms, hermitianforms );
        else
            forms.maybeFrobenius := false;
        fi;
    fi;
    # if all forms are excluded then we are finished
    if not forms.maybeDual and not forms.maybeFrobenius  then
            Append( forms.invariantforms, [ "linear" ] );
    fi;

    return forms;
end );

#############################################################################
##
#O  PreservedForms( <grp> )
##    returns (i) quadratic form(s) if it has one,
##            (ii) a sesquilinear form(s) otherwise
##  it basically converts the information given by Frank Cellers operation
## to output we want...
##
InstallMethod( PreservedForms,
    "for a matrix group over a finite field",
    [ IsMatrixGroup ],
    function( grp )
    local newforms, forms, y, newform, i, field;
    newforms := [];
    field := DefaultFieldOfMatrixGroup(grp);
    forms := PreservedFormsOp(grp);
    if IsList(forms) and IsEmpty(forms) then
        return forms;
    fi;
    for y in forms!.invariantforms do
       if y[1] in ["orthogonalplus", "orthogonalminus", "orthogonalcircle"] then
          newform := QuadraticFormByMatrix(y[4], field);
          Info(InfoForms, 1, Concatenation("preserved up to the following scalars: ", String(y[3])) );
          Info(InfoForms, 1, y[1] );
          Add( newforms, newform );
       elif y[1] = "symplectic" then
          newform := BilinearFormByMatrix(y[2], field);
          Add( newforms, newform );
       elif y[1] = "unitary" then
          newform := HermitianFormByMatrix(y[2], field);
          Add( newforms, newform );
       elif y[1] = "linear" then
          i := NrRows(One(grp));
          newform := BilinearFormByMatrix( NullMat(i,i,field), field );
          Add( newforms, newform );
       fi;
    od;
    return newforms;
end );

#############################################################################
##
#O  PreservedSesquilinearForms( <grp> )
##    returns a sesquilinear form(s) if it has one
##  it basically converts the information given by Frank Cellers operation
## to output we want...
## to do: check whether this is still usefull (and even correctness)
##
InstallMethod( PreservedSesquilinearForms,
    "for a matrix group over a finite field",
    [ IsMatrixGroup ],
    function( grp )
    local newforms, forms, y, newform, i, field;
    newforms := [];
    field := DefaultFieldOfMatrixGroup(grp);
    forms := PreservedFormsOp(grp);
    if IsList(forms) and IsEmpty(forms) then
        return forms;
    fi;
    for y in forms!.invariantforms do
       if y[1] in ["symplectic", "orthogonalplus",
                   "orthogonalminus", "orthogonalcircle"] then
          newform := BilinearFormByMatrix(y[2], field);
          Add( newforms, newform );
       elif y[1] = "unitary" then
          newform := HermitianFormByMatrix(y[2], field);
          Add( newforms, newform );
       elif y[1] = "linear" then
          i := NrRows(One(grp));
          newform := BilinearFormByMatrix( NullMat(i,i,field), field );
          Add( newforms, newform );
       fi;
    od;
    return newforms;
end );

#############################################################################
##
#O  PreservedQuadraticForms( <grp> )
##    returns  quadratic form(s) if it has one.
##  it basically converts the information given by Frank Cellers operation
## to output we want...
## to do: check whether this is still usefull (and even correctness)
##
InstallMethod( PreservedQuadraticForms,
    "for a matrix group over a finite field",
    [ IsMatrixGroup ],
    function( grp )
    local newforms, forms, y, newform, i, field;
    newforms := [];
    field := DefaultFieldOfMatrixGroup(grp);
    forms := PreservedFormsOp(grp);
    for y in forms!.invariantforms do
       if y[1] in ["orthogonalplus", "orthogonalminus", "orthogonalcircle"] then
          newform := QuadraticFormByMatrix(y[4], field);
          Info(InfoForms, 1, Concatenation("preserved up to the following scalars: ", String(y[3])) );
          Info(InfoForms, 1, y[1] );
          Add( newforms, newform );
        fi;
    od;
    return newforms;
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

#For compatibility reasons with recog
#############################################################################
##
#F  ClassicalForms_InvariantFormDual( <module>, <dmodule> )
##
ClassicalForms_InvariantFormDual := function( module, dmodule )
    local   hom,  scalars,  form,  iform,  identity,  field,  root,
            q,  i,  m,  a,  quad,  sgn;

    # <dmodule> acts absolutely irreducible without scalars
    hom := MTX.Homomorphisms( dmodule, DualGModule(dmodule) );
    if 0 = Length(hom)  then
        return false;
    elif 1 < Length(hom)  then
        Error( "module acts absolutely irreducibly but two forms found" );
    fi;
    Info( InfoForms, 1, "found homomorphism between V and V^*\n" );

    # make sure that the forms commute with the generators of <module>
    scalars  := [];
    form     := hom[1];
    iform    := form^-1;
    identity := One(form);
    field   :=  MTX.Field(module);
    root    :=  PrimitiveRoot(field);
    q        := Size(field);
    for i  in MTX.Generators(module)  do
        m := i * form * TransposedMat(i) * iform;
        a := m[1,1];
        if m <> a*identity  then
            Info(InfoForms, 1,
                "form is not invariant under all generators\n" );
            return false;
        fi;
        a := NthRoot(field,a,2);
        Add( scalars, a );
    od;

    # check the type of form
    if TransposedMat(form) = -form  then
        Info(InfoForms, 1, "form is symplectic\n" );
        if Characteristic(field) = 2  then
            quad := ClassicalForms_QuadraticForm2(
                field, form, MTX.Generators(module), scalars );
            if quad = false  then
                return [ "symplectic", form, scalars ];
            elif MTX.Dimension(module) mod 2 = 1  then
                Error( "no quadratic form but odd dimension" );
            elif ClassicalForms_Signum2( field, form, quad ) = -1  then
                return [ "orthogonalminus", form, scalars, quad ];
            else
                return [ "orthogonalplus", form, scalars, quad ];
            fi;
        else
            return [ "symplectic", form, scalars ];
        fi;
    elif TransposedMat(form) = form  then
        Info(InfoForms, 1, "form is symmetric\n" );
        quad := ClassicalForms_QuadraticForm( field, form );
        if MTX.Dimension(module) mod 2 = 1  then
            return [ "orthogonalcircle", form, scalars, quad ];
        else
            sgn := ClassicalForms_Signum( field, form, quad );
            if sgn[1] = -1  then
                return [ "orthogonalminus", form, scalars, quad, sgn[2] ];
            else
                return [ "orthogonalplus", form, scalars, quad, sgn[2] ];
            fi;
        fi;
    else
        Info( InfoForms, 1,"unknown form\n" );
        return [ "unknown", "dual", form, scalars ];
    fi;
end;

#For compatibility reasons with recog
#############################################################################
##
#F  ClassicalForms_InvariantFormFrobenius( <module>, <fmodule> )
##
TransposedFrobeniusMat := function( mat, qq )
    local   i,  j;
    mat:=MutableTransposedMat(mat);
    for i  in [ 1 .. NrRows(mat) ]  do
        for j  in [ 1 .. NrCols(mat) ]  do
            mat[i,j] := mat[i,j]^qq;
        od;
    od;
    return mat;
end;

#For compatibility reasons with recog
#############################################################################
##
#F ClassicalForms_InvariantFormFrobenius( module, fmodule )
##
ClassicalForms_InvariantFormFrobenius := function( module, fmodule )
    local   fro,  hom,  form,  q,  qq,  k,  a,  scalars,  iform,
            identity,  field,  root,  i,  m,  j;

    # <fmodule> acts absolutely irreducible without scalars
    fro := DualFrobeniusGModule(fmodule);
    hom := MTX.Homomorphisms(fmodule, fro);
    if 0 = Length(hom)  then
        return false;
    elif 1 < Length(hom)  then
        Error( "module acts absolutely irreducibly but two form found" );
    fi;
    Info( InfoForms, 1,"found homomorphism between V and (V^*)^frob\n" );

    # invariant form might return a scalar multiple of our form
    field    := MTX.Field(module);
    form := hom[1];
    q  := Size(field);
    qq := Characteristic(field)^(DegreeOverPrimeField(field)/2);
    k  := PositionNonZero(form[1]);
    a  := form[1,k] / form[k,1]^qq;
    a := NthRoot(field,a,(1-qq) mod (q-1));
    if a = fail then
      return false;
    fi;
    form := form * a^-1;


    # make sure that the forms commute with the generators of <module>
    scalars  := [];
    iform    := form^-1;
    identity := One(form);
    root     := PrimitiveRoot(field);
    for i  in MTX.Generators(module)  do
        m := i * form * TransposedFrobeniusMat(i,qq) * iform;
        a := m[1,1];
        if m <> a*identity  then
            Info(InfoForms, 1,
                 "form is not invariant under all generators\n" );
            return false;
        fi;
        a:=NthRoot(field,a,qq+1);
        Add( scalars, a );
    od;

    # check the type of form
    for i  in [ 1 .. NrRows(form) ]  do
        for j  in [ 1 .. NrRows(form) ]  do
            if form[i,j]^qq <> form[j,i]  then
                Info(InfoForms, 1, "unknown form\n" );
                return [ "unknown", "frobenius", form, scalars ];
            fi;
        od;
    od;
    return [ "unitary", form, scalars ];

end;
