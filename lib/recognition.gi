#############################################################################
##
##  recognition.gi              'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##                                                              Frank Celler
##  Copyright 2017, Vrije Universiteit Brussel
##  Copyright 2017, The University of Western Austalia
##  Copyright (C) 1996,  Lehrstuhl D fuer Mathematik, RWTH Aachen, Germany
##
##  This file contains functions to compute sesquilinear forms invariant under
##  matrix groups.
##
##  *** Bamberg and De Beule are very grateful to Frank Celler for
##  generously providing the bulk of this code.  ***
##
#############################################################################
##
#F  ClassicalForms_ScalarMultipleFrobenius( <field>, <mat> )
##
InstallGlobalFunction( ClassicalForms_ScalarMultipleFrobenius,
 function( F, M )
    local mpol, d, c, z, I, q, qq, t, a, l, i0, Minv, tM, tMi;

    # compute the characteristic polynomial of <M>
    mpol := CharacteristicPolynomial(M);

    # get the position of the non-zero coefficients
    d := Degree(mpol);
    c := CoefficientsOfUnivariatePolynomial(mpol);
    z := Zero(F);

    I := Filtered( [0 .. d],  x -> not IsZero(c[x+1]) );
    q := Size(F);
    qq := Characteristic(F)^(DegreeOverPrimeField(F)/2);

    Minv := M^-1;
    tM := Trace(M); tMi := Trace(Minv)^qq; # Frobenius of trace
    if (IsZero(tM) and not IsZero(tMi)) or
       (not IsZero(tM) and IsZero(tMi)) then
        return false;
    fi;

    # make sure that <d> and <d>-i both occur
    if ForAny( I, x ->  IsZero(c[d-x+1]) )  then
        return false;
    fi;

    Add(I,q-1);
    # we need gcd one in order to get alpha exactly (ignoring +-)
    t := GcdRepresentation(I);
    i0 := I*t; #note that this is just Gcd(I)

    if not IsZero(tM) and not IsZero(tMi) then
        return [i0, tM/tMi];
    fi;

    # compute gcd representation
    a:=c[1];
    l:=List([1..Length(I)-1], x ->(a*c[d-I[x]+1]^qq/c[I[x]+1]));

    a:=Product([1..Length(I)-1], x->l[x]^t[x]);
    # Now the scalar $\lambda$ satisfies $\lambda^{i_0}=a$

    # check: $\forall_i: \bar c_{d-i}c_0=c_i\lambda^i$
    if ForAny([1..Length(I)-1],x->l[x]<>a^QuoInt(I[x],i0)) then
        Info( InfoForms, 1,
         "characteristic polynomial does not reveal scalar\n" );
      return false;
    fi;

    # compute a square root of <alpha>
    a:=NthRoot(F,a,(qq+1)*i0);
    if a=fail then
        Info( InfoForms, 1,
        "characteristic polynomial does not reveal scalar\n" );
      return false;
    fi;
    return [i0,a]; #We will compute all solutions in the new version of the code for the next function!
  end );


#############################################################################
##
#F  ClassicalForms_GeneratorsWithoutScalarsFrobenius( grp )
##
InstallGlobalFunction( ClassicalForms_GeneratorsWithoutScalarsFrobenius,
 function( grp )
    local tries, gens, field, m1, a1, new, i;

    # start with 2 random elements,  at most 10 tries
    tries := 0;
    gens  := [];
    field := FieldOfMatrixGroup(grp);
    while Length(gens) < 2  do
        tries := tries + 1;
        if tries = 11  then return false;  fi;
        m1 := PseudoRandom(grp);
        a1 := ClassicalForms_ScalarMultipleFrobenius(field,m1);
        if IsList(a1) and a1[1]=1 then
            a1:=a1[2];
            Add(gens, m1*a1^-1);
        fi;
    od;
    new := GModuleByMats( gens, field );

    # the module must act absolutely irreducibly
    while not MTX.IsAbsolutelyIrreducible(new)  do
        for i  in [ 1 .. 2 ]  do
            repeat
                tries := tries + 1;
                if tries > 10  then return false;  fi;
                m1 := PseudoRandom(grp);
                a1 := ClassicalForms_ScalarMultipleFrobenius(field,m1);
            until IsList(a1) and a1[1]=1;
            a1:=a1[2];
            Add( gens, m1*a1^-1 );
        od;

        new := GModuleByMats( gens, field );
    od;

    return new;
 end );


ClassicalForms_GeneratorsWithoutScalarsFrobeniusNew := 
 function( grp )
    local tries, gens, field, m1, a1, new, i;

    # start with 2 random elements,  at most 10 tries
    tries := 0;
    gens  := [];
    field := FieldOfMatrixGroup(grp);
    #This function will generate all possible Frobenius dual modules that need to be checked.
    for m1 in GeneratorsOfGroup(grp) do 
    #while Length(gens) < 2  do
        #tries := tries + 1;
        #if tries = 11  then return false;  fi;
        #m1 := PseudoRandom(grp);
        a1 := ClassicalForms_ScalarMultipleFrobenius(field,m1);
        #compute all solutions of lambda^a1[1] = a1[2] where a1[2] is defined as in the paper. 
        if IsList(a1) and a1[1]=1 then
            a1:=a1[2];
            Add(gens, a1*TransposedMat(m1^-1));
        fi;
    od;
    new := GModuleByMats( gens, field );
    return new; 
 end; #);

#############################################################################
##
#F  ClassicalForms_PossibleScalarsBilinear( <field>, <mat> )
# Very important: this function is meant to be called *only* for
# bilinear forms, that means, the "bar map" is the identity.
# This function returns [i0,a] such that if <mat> preserves a bilinear form B 
# modulo scalars then the scalar lambda for which mat*B*mat^T = lambda*B 
# satisfies lambda^i0 = a
##
InstallGlobalFunction( ClassicalForms_PossibleScalarsBilinear,
 function( F, M )
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

    #a:=c[1]; #this caused a mess ;-)

    #Lemma 1: if g preserves a form modulo lambda then c[d-I]\bar(c_0) = \bar(c_i) lambda^i.
    #cpol = c_0 +c_1 X + \ldots + c_d X^d
    #l will be the list of the lambda^i's. 
    l:=List([1..Length(I)-1], x ->(c[1]*c[d-I[x]+1]/c[I[x]+1])); #here is the line with identity field automorphism!

    a:=Product([1..Length(I)-1], x->l[x]^t[x]);
    # Now the scalar $\lambda$ satisfies $\lambda^{i_0}=a$

    # check: $\forall_i: c_{d-i}c_0=c_i\lambda^i
    # if any of the conditions from Lemma 1 fail, g does not preserve the form modulo this lambda. 
    if ForAny([1..Length(I)-1],x->(l[x]<>a^QuoInt(I[x],i0))) then
        Info( InfoForms, 1, "Characteristic polynomial proves that there is no scalar\n" );
      return false;
    fi;

    # compute a square root of <alpha>
    #a:=NthRoot(F,a,2*i0);
    #if a=fail then
    #    Info( InfoForms, 1,"characteristic polynomial does not reveal scalar\n" );
    #  return fail; #was false
    #fi;
    return [i0,a];
  end );


#############################################################################
##
#F  ClassicalForms_GeneratorsWithoutScalarsDual( grp )
#F  ClassicalForms_GeneratorsWithBetterScalarsBilinear( grp ) #this is the objective!
#
##
InstallGlobalFunction( ClassicalForms_GeneratorsWithoutScalarsDual,
  function( grp )
    local tries, gens, field, m1, a1, new, i, scalars, root, improvegenerator, res, newgens, champion, len;

    # the aim of this function is to replace the matrix m1 by a 
    # matrix that has as few solutions to the scalar equation 
    # lambda^a1[1] = a1[2] as possible. It checks first if a1[1] = 1, 
    # since then lambda is determined. Next we check if a1[2] is a 
    # square. And then we replace m1 by a matrix that leaves the 
    # bilinear form invariant modula the scalar 1. If none of these 
    # are possible, we try to replace m1 by a matrix that has fewer 
    # solutions to the scalar equation. 

    improvegenerator := function(m1,i,count,len)
        local a1, s, j, k, scalars;
        
        a1 := ClassicalForms_PossibleScalarsBilinear(field,m1);
        if a1 = false then
            return false; #The group does not preserve a form modulo scalars
        fi;
        
        #scalars := AsList(Group(NthRoot(field,a1[2],a1[1]))); # add all possible scalars for m1 
        
        if IsList(a1) then 
            root := NthRoot(field,a1[2],a1[1]);
            if a1[1] = 1 then # the matrix m1 has scalar a1[2]
                return [m1,[a1[2]]];
            elif IsEvenInt(LogFFE(root,PrimitiveElement(field))) then #IsEvenInt(a1[1]) then # modify m1 to a matrix that has scalar one.
                return [m1/NthRoot(field,root,2),[One(field)]];
            else
                scalars := AsList(Group(NthRoot(field,a1[2],a1[1]))); # add all possible scalars for m1 
                if count = 0 then 
                    return champion; # add all possible scalars for m1 
                elif Length(scalars) < len then
                    champion := [m1,scalars];
                    len := Length(scalars);
                fi;
                k := Random(Difference([1..Length(gens)],[i]));
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
    field := FieldOfMatrixGroup(grp);
 
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
##
InstallGlobalFunction( ClassicalForms_Signum2,
  function( field, form, quad )
    local dim, base, avoid, i, d, j, c, k, x, sgn, pol;

    # compute a new basis,  such that the symmetric form is standard
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
  ## This could be replaced by operations in "Forms"
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
##
InstallGlobalFunction( ClassicalForms_QuadraticForm2,
 function( field, form, gens, scalars )
    local dim, H, i, j, e, b, y, x, r, l;

    # raise an error if char is not two
    if Characteristic(field) <> 2  then
        Error( "characteristic must be two" );
    fi;

    # construct the upper half of the form
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
        for i  in [ 1 .. dim ]  do
            for j  in [ i+1 .. dim ]  do
                if not IsZero(r[i,j]+r[j,i]) then
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
    e := SolutionMat( TransposedMat(e), b );
    if e <> fail  then
        for i  in [ 1 .. dim ]  do
            H[i,i] := e[i];
        od;
        return ImmutableMatrix(field,H);
    else
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


#############################################################################
##
#F  ClassicalForms_InvariantFormDual( <gens_scalars> )
#   We generate the module by generators stored in gens_scalars and we know
#   that this module is the module of the original group. Then we construct
#   all possible dual modules according to the list stored in gens_scalars. 
#   Then we check which ones yield bilinear forms.
##
ClassicalForms_InvariantFormDual := function( gens_scalars )
    local   hom,  scalars,  form,  iform,  identity,  field,  root,
            q,  i,  m,  a,  quad,  sgn, gmodule, forms, biglist, x,
            dmodule, gens, scale_gens, module, output, isform, pair;

    # <dmodule> acts absolutely irreducible without scalars
    field := FieldOfMatrixGroup(Group(gens_scalars[1]));
    gmodule := GModuleByMats(gens_scalars[1],field);
    forms := [];
    biglist := Cartesian(gens_scalars[2]);
    gens := gens_scalars[1];
    for x in biglist do 
        scale_gens := List([1..Length(gens)],i->gens[i] / x[i]);
        dmodule := DualGModule(GModuleByMats(scale_gens,field)); 
        hom := MTX.Homomorphisms( gmodule, dmodule);
        #if 0 = Length(hom)  then
        #    return false;
        #elif 1 < Length(hom)  then
        #    Error( "module acts absolutely irreducibly but two forms found" );
        #fi;
        Info( InfoForms, 1, "found homomorphism between V and V^*\n" );
        Append(forms,List(hom,i->[i,x]));
    od;
    #TO DO: get rid of these checks in the final version.
    # make sure that the forms commute with the generators of <module>
    # forms is a list now. We should run the following checks for each of the elements in forms.
    
    output := [];
    root    :=  PrimitiveRoot(field);
    q        := Size(field);

    scalars  := gens_scalars[2];

    #

    for pair in forms do
        form := pair[1];
        scalars  := pair[2];
        iform    := form^-1;
        identity := One(form);
        isform := true;
        for i in [1..Length(gens)] do
            m := gens[i] * form * TransposedMat(gens[i]) * iform;
            a := m[1,1]; 
#            Print(a=scalars[i],"\n");
            if m <> a*identity  then
                Info(InfoForms, 1,
                    "form is not invariant under all generators\n" );
                isform := false;
                continue;
            fi;
         od;

        # check the type of form, replace this afterwards with default forms constructors. Then also join the two for loops.
        
        if isform = false then
            continue;
        elif TransposedMat(form) = -form  then
            Info(InfoForms, 1, "form is symplectic\n" );
            if Characteristic(field) = 2  then
                quad := ClassicalForms_QuadraticForm2(
                field, form, gens, scalars );
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
    od;
    return output;
end;


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

DualFrobeniusGModule := function(module)
local   F,  k,  dim,  mats,  dmats,  qq,  i,  j,  l;

  if SMTX.IsZeroGens(module) then
    return GModuleByMats([],module.dimension,SMTX.Field(module));
  else
    F := MTX.Field(module);
    k := DegreeOverPrimeField(F);
    if k mod 2 = 1  then
        Error( "field <F> is not a square" );
    fi;
    dim   := MTX.Dimension(module);
    mats  := MTX.Generators(module);
    dmats := List(mats,i->List(i,ShallowCopy));
    qq    := Characteristic(F) ^ ( k / 2 );
    for i  in [ 1 .. Length(mats) ]  do
      for j  in [ 1 .. dim ]  do
        for l  in [ 1 .. dim ]  do
          dmats[i][j,l] := mats[i][l,j]^qq;
        od;
      od;
      dmats[i]:=ImmutableMatrix(F,dmats[i]);
    od;

    return GModuleByMats(List(dmats,i->i^-1),F);
  fi;
end;


ClassicalForms_InvariantFormFrobenius := function( module, fmodule )
    local   fro,  hom,  form,  q,  qq,  k,  a,  scalars,  iform,
            identity,  field,  root,  i,  m,  j;

    # <fmodule> acts absolutely irreducible without scalars
    fro := DualFrobeniusGModule(fmodule);
    hom := MTX.Homomorphisms(fmodule, fro);
    if 0 = Length(hom)  then
        return false;
    elif 1 < Length(hom)  then
        Error( "module acts absolutely irreducibly but two forms found" );
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

# forms is  a record  which stores information  about which  forms the
# matrix  group grp  leaves  invariant.   It has  to  have the  record
# components .maybeDual and  .maybeFrobenius.  It sets .maybeDual
# to false if the cpoly c of  an element of the group grp excludes the
# possibility  that grp  leaves  invariant a  bilinear  form and  sets
# .maybeFrobenius to false  if the cpoly c of an  element of the group
# grp   excludes  the   possibility  that   grp  leaves   invariant  a
# sesquilinear form

InstallGlobalFunction(PossibleClassicalForms,
                       function(arg)

    local  I, d, z, f, t, i0, a, l, g, q, qq, forms, grp,  c, tM, tMi, Minv;

    if not Length(arg) in [2,3] then
        Error("Usage: PossibleClassicalForms( grp, g [,forms] )");
    fi;
    grp   := arg[1];
    g     := arg[2];
    if Length(arg) = 3 then
        forms := arg[3];
    else forms := rec();
         forms.maybeDual := true;
         f := DefaultFieldOfMatrixGroup(grp);
         forms.field := f;
         forms.maybeFrobenius := DegreeOverPrimeField(f) mod 2 = 0;
    fi;

    c := CharacteristicPolynomial(g);
    c := CoefficientsOfUnivariatePolynomial(c);

    d := DimensionOfMatrixGroup(grp);
    f := forms.field;
    z := Zero(f);
    q := Size(f);
    Minv := g^-1;
    tM := Trace(g); tMi := Trace(Minv);

    if forms.maybeFrobenius then
        qq := Characteristic(f)^(DegreeOverPrimeField(f)/2);
        tMi := tMi^qq;               # Frobenius of trace
    fi;

    if (IsZero(tM) and not IsZero(tMi)) or
       (not IsZero(tM) and IsZero(tMi)) then
        forms.maybeFrobenius := false;
        forms.maybeDual := false;
        return false;
    fi;

    I := Filtered( [0 .. d],  x -> (not IsZero(c[x+1]) ));
    # make sure that <d> and <d>-i both occur
    if ForAny( I, x ->  IsZero(c[d-x+1]) )  then
        forms.maybeFrobenius := false;
        forms.maybeDual := false;
        return false;
    fi;
    Add(I,q-1);

    # we need gcd one in order to get alpha exactly (ignoring +-)
    t := GcdRepresentation(I);
    i0:=I*t;

    if not IsZero(tM) and not IsZero(tMi) then
        return [i0, tM/tMi];
    fi;


    if forms.maybeDual  then
        a  := c[1];
        l  := List( [1..Length(I)-1],
                    x ->(a*c[d-I[x]+1]/c[I[x]+1]) );
        g  := Product( [1..Length(I)-1],x -> l[x]^t[x] );
        if ForAny( [1..Length(I)-1], x -> l[x]<>g^(I[x]/i0) )  then
            forms.maybeDual := false;
        fi;
    fi;
    if forms.maybeFrobenius  then
        a  := c[1];
        l  := List([1..Length(I)-1], x ->
                   (a*c[d-I[x]+1]^qq/c[I[x]+1]));
        g  := Product( [1..Length(I)-1],x -> l[x]^t[x] );
        if ForAny( [1..Length(I)-1], x -> l[x]<>g^(I[x]/i0) )  then
            forms.maybeFrobenius := false;
        fi;
    fi;
    if forms.maybeDual = false and forms.maybeFrobenius = false then
        return false;
    else
        return true;
    fi;
end);

#############################################################################
##
#F  PreservedSesquilinearForms( <grp> )
##
##    returns a list of forms
##
#InstallMethod( PreservedSesquilinearForms, [ IsMatrixGroup ],
#  function( grp )
#    local   forms, field, i, g, qq, c, module,  invariantforms,
#            dmodule, fmodule, form, y, newform, newforms;
#    field := DefaultFieldOfMatrixGroup(grp);

#    forms := rec();
#    forms.field := field;
#    forms.invariantforms := [];
#    forms.maybeDual      := false;
#    forms.maybeFrobenius := false;

#    # set up the module and other information
#    module := GModuleByMats(GeneratorsOfGroup(grp), field);

    # set the possibilities
#    forms.maybeDual      := true;
#    forms.maybeFrobenius := DegreeOverPrimeField(field) mod 2 = 0;

#    if forms.maybeFrobenius  then
#        qq := Characteristic(field)^(DegreeOverPrimeField(field)/2);
#    fi;

    # We first perform some inexpensive tests with a few random elements
#    for i in [1 .. 8]  do
#        g := PseudoRandom(grp);
#        if forms.maybeDual or forms.maybeFrobenius  then
#           PossibleClassicalForms( grp, g, forms );
#        fi;
#    od;

    # if all forms are excluded then we are finished
#    if not forms.maybeDual and not forms.maybeFrobenius  then
#       i := NrRows(One(g));
#        return [ BilinearFormByMatrix( NullMat(i,i,field), field ) ];
#    fi;

    # <grp> must act absolutely irreducibly
#    if not MTX.IsAbsolutelyIrreducible(module)  then
#        Info( InfoForms, 1,  "grp not absolutely irreducible\n" );
#        return [];
#    fi;

    # try to find generators without scalars
#    if forms.maybeDual  then
#        dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(grp);
#        if dmodule = false  then
#            forms.maybeDual := false;
#        fi;
#    fi;
#    if forms.maybeFrobenius  then
#        fmodule := ClassicalForms_GeneratorsWithoutScalarsFrobenius(grp);
#        if fmodule = false  then
#            forms.maybeFrobenius := false;
#        fi;
#    fi;

    # now try to find an invariant form
#    if forms.maybeDual  then
#        form := ClassicalForms_InvariantFormDual(module,dmodule);
#        if form <> false  then
#            Add( forms.invariantforms, form );
#        else
#            forms.maybeDual := false;
#        fi;
#    fi;

#    if forms.maybeFrobenius  then
#        form := ClassicalForms_InvariantFormFrobenius(module,fmodule);
#        if form <> false  then
#            Add( forms.invariantforms, form );
#        else
#            forms.maybeFrobenius := false;
#        fi;
#    fi;
    # if all forms are excluded then we are finished
#    if not forms.maybeDual and not forms.maybeFrobenius  then
#            Add( forms.invariantforms, [ "linear" ] );
#    fi;

    ## We can convert the information Frank Celler wanted
    ## to output we want...

#    newforms := [];

#    for y in forms!.invariantforms do
#       if y[1] in ["symplectic", "orthogonalplus",
#                   "orthogonalminus", "orthogonalcircle"] then
#          newform := BilinearFormByMatrix(y[2], field);
#          Add( newforms, newform );
#       elif y[1] = "unitary" then
#          newform := HermitianFormByMatrix(y[2], field);
#          Add( newforms, newform );
#       elif y[1] = "linear" then
#          i := NrRows(One(g));
#          newform := BilinearFormByMatrix( NullMat(i,i,field), field );
#          Add( newforms, newform );
#       fi;
#    od;
#
#    return newforms;
#end );

#############################################################################
##
#O  ScalarOfSimilarity( <grp> )
##
##    returns a scalar
##
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
            gens_scalars, fmodule, form, y, newform, newforms, bilinearforms;
    
    field := DefaultFieldOfMatrixGroup(grp);

    # remember to rename Dual to Bilinaer and Frobenius to Hermitian.

    forms := rec();
    forms.field := field;
    forms.invariantforms := [];
    forms.maybeDual      := false;
    forms.maybeFrobenius := false;

    # set up the module and other information
    module := GModuleByMats(GeneratorsOfGroup(grp), field);

    # set the possibilities
    forms.maybeDual      := true;
    forms.maybeFrobenius := IsEvenInt(DegreeOverPrimeField(field));

    if forms.maybeFrobenius  then
        qq := Characteristic(field)^(DegreeOverPrimeField(field)/2);
    fi;

    # We first perform some inexpensive tests with a few random elements
    #for i in [1 .. 8]  do
    #    g := PseudoRandom(grp);
    #    if forms.maybeDual or forms.maybeFrobenius  then
    #       PossibleClassicalForms( grp, g, forms );
    #    fi;
    #od;

    # if all forms are excluded then we are finished
    #if not forms.maybeDual and not forms.maybeFrobenius  then
    #    i := NrRows(One(g));
    #    return [ BilinearFormByMatrix( NullMat(i,i,field), field ) ];
    #fi;

    # <grp> must act absolutely irreducibly
    if not MTX.IsAbsolutelyIrreducible(module)  then
        Error("Currently the use of MeatAxe requires the module to be absolutely irreducible");
        #Info( InfoForms, 1,  "grp not absolutely irreducible\n" );
        #return [];
    fi;

    # try to find generators without scalars
    if forms.maybeDual then
        gens_scalars := ClassicalForms_GeneratorsWithoutScalarsDual(grp);
        if gens_scalars = false then #happens if the group does not preserve a bilinear form modulo scalars
            forms.maybeDual := false;
        fi;
    fi;
    if forms.maybeFrobenius  then
        fmodule := ClassicalForms_GeneratorsWithoutScalarsFrobenius(grp);
        if fmodule = false  then
            forms.maybeFrobenius := false;
        fi;
    fi;

    # now try to find an invariant form
    if forms.maybeDual  then
        bilinearforms := ClassicalForms_InvariantFormDual(gens_scalars);
        if bilinearforms <> false  then
            Append( forms.invariantforms, bilinearforms );
        else
            forms.maybeDual := false;
        fi;
    fi;

    if forms.maybeFrobenius  then
        form := ClassicalForms_InvariantFormFrobenius(module,fmodule);
        if form <> false  then
            Add( forms.invariantforms, form );
        else
            forms.maybeFrobenius := false;
        fi;
    fi;
    # if all forms are excluded then we are finished
    if not forms.maybeDual and not forms.maybeFrobenius  then
            Add( forms.invariantforms, [ "linear" ] );
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
##
InstallMethod( PreservedSesquilinearForms,
    "for a matrix group over a finite field",
    [ IsMatrixGroup ],
    function( grp )
    local newforms, forms, y, newform, i, field;
    newforms := [];
    field := DefaultFieldOfMatrixGroup(grp);
    forms := PreservedFormsOp(grp);
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
