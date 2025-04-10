#############################################################################
##
##  forms.gi              'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##
##  Copyright 2024, Vrije Universiteit Brussel
##  Copyright 2024, The University of Western Austalia
##
##  Implementation of quadratic and sesquilinear forms
##
#############################################################################

#############################################################################
# Fundamental methods:
#############################################################################

InstallMethod( \=, "for two trivial forms",
  [IsTrivialForm and IsFormRep, IsTrivialForm and IsFormRep],
  function( a, b )
    return a!.basefield = b!.basefield;
  end );

InstallMethod( \=, "for two forms",
  [IsForm and IsFormRep, IsForm and IsFormRep],
  function( a, b )
    return (a!.basefield = b!.basefield) and
           (a!.type = b!.type) and
           (a!.matrix = b!.matrix);
  end );

#############################################################################
# Constructor methods
# Form-by-matrix methods
#############################################################################
#############################################################################
#O  FormByMatrix( <mat>, <field>, <string> )
#   A general constructor not for users, not documented.
##
InstallMethod( FormByMatrix, "for a ffe matrix, a field and a string",
  [IsMatrix and IsFFECollColl, IsField and IsFinite, IsString],
  function( m, f, string )
    local el;
    m := ImmutableMatrix(f, m);
    el := rec( matrix := m, basefield := f, type := string );

   ## We follow a certain convention, which is outlined in the manual,
   ## in order to determine from a Gram matrix the type of the constructed
   ## form.

    if string = "hermitian" then
      if not IsInt(Sqrt(Size(f))) then
        Error("No hermitian form exist when the order of <f> is not a square" );
      fi;
      if FORMS_IsHermitianMatrix(m,f) then
        Objectify(NewType( HermitianFormFamily ,  IsFormRep),  el);
        return el;
      else
        Error("Given matrix does not define a hermitian form" );
      fi;
    elif string = "symplectic" then
      if FORMS_IsSymplecticMatrix(m,f) then
        Objectify(NewType( BilinearFormFamily ,  IsFormRep),  el);
        return el;
      else
        Error("Given matrix does not define a symplectic form" );
      fi;
    elif string = "orthogonal" then
      if Characteristic(f) = 2 then
        Error("No orthogonal forms exist in even characteristic" );
      fi;
      if FORMS_IsSymmetricMatrix(m) then
        Objectify(NewType( BilinearFormFamily ,  IsFormRep),  el);
        return el;
      else
        Error("Given matrix does not define an orthogonal form" );
      fi;
    elif string = "pseudo" then
      if IsOddInt(Size(f)) then
        Error("No pseudo forms exist in even characteristic" );
      fi;
      if (FORMS_IsSymmetricMatrix(m) and (not FORMS_IsSymplecticMatrix(m))) then
        Objectify(NewType( BilinearFormFamily ,  IsFormRep),  el);
        return el;
      else
        Error("Given matrix does not define a pseudo-form" );
      fi;
    elif string = "quadratic" then
      el.matrix := ImmutableMatrix(f, Forms_RESET(m));
      Objectify(NewType( QuadraticFormFamily ,  IsFormRep),  el);
      return el;
    else
      Error("Please specify a form properly");
    fi;
  end );

#############################################################################
#O  BilinearFormByMatrixOp( <m>, <f> )
#   A less general constructor not for users, not documented.
# <f> finite field, <m> orthogonal or symplectic matrix over <f>.
# zero matrix is allowed, then the trivial form is returned.
# if <m> is not zero, it is checked that <m> determines a bilinear form which is
# symplectic, orthogonal, pseudo or hermitian
##
InstallMethod( BilinearFormByMatrixOp, "for a ffe matrix and a field",
  [IsMatrix and IsFFECollColl, IsField and IsFinite],
  function( m, f )
    local el, n;
    n := NrRows(m);
    m := ImmutableMatrix(f, m);
    if IsZero(m) then
       el := rec( matrix := m, basefield := f, type := "trivial", vectorspace := FullRowSpace(f,n) );
       Objectify(NewType( TrivialFormFamily ,  IsFormRep),  el);
       return el;
    elif FORMS_IsSymplecticMatrix(m,f) then
       el := rec( matrix := m, basefield := f, type := "symplectic", vectorspace := FullRowSpace(f,n) );
       Objectify(NewType( BilinearFormFamily ,  IsFormRep),  el);
       return el;
    elif FORMS_IsSymmetricMatrix(m) and Characteristic(f) <> 2 then
       el := rec( matrix := m, basefield := f, type := "orthogonal", vectorspace := FullRowSpace(f,n) );
       Objectify(NewType( BilinearFormFamily ,  IsFormRep),  el);
       return el;
    elif FORMS_IsSymmetricMatrix(m) and Characteristic(f) = 2 then
       el := rec( matrix := m, basefield := f, type := "pseudo", vectorspace := FullRowSpace(f,n) );
       Objectify(NewType( BilinearFormFamily ,  IsFormRep),  el);
       return el;
    else
       Error("Invalid Gram matrix");
    fi;
  end );

# 20/01/2016: small change here: added MutableCopyMat, since it is not logical
# that changing the matrix that was used for construction afterwards changes
# the form. See test_forms12.g
#############################################################################
#O BilinearFormByMatrix( <m>, <f>): constructor for users.
# <f> finite field, <m> orthogonal or symplectic matrix over <f>.
# checks whether <m> is over <f>.
# zero matrix is allowed, then the trivial form is returned.
##
InstallMethod( BilinearFormByMatrix, "for a ffe matrix and a field",
  [IsMatrix and IsFFECollColl, IsField and IsFinite],
  function( m, f )
    local gf;
    gf := DefaultFieldOfMatrix(m);
    if not PrimitiveElement(gf) in f then
      Error("<m> is not a matrix over <f>");
    fi;
    return BilinearFormByMatrixOp( m, f);
end );

#############################################################################
#O BilinearFormByMatrix( <m> ): constructor for users. uses above constructor
# finite field <f> is determined from <m>
##
InstallMethod( BilinearFormByMatrix, "for a ffe matrix ",
  [IsMatrix and IsFFECollColl ],
  function( m )
  local f;
  f := DefaultFieldOfMatrix(m);
  return BilinearFormByMatrixOp( m, f);
end );

#############################################################################
#O QuadraticFormByMatrixOp( <m>, <f> ): constructor not for users.
#  Analogous to BilinearFormByMatrixOp
# <f> finite field, <m> matrix over <f>.
# <m> is always "Forms_RESET" to an upper triangle matrix.
# zero matrix is allowed, then the trivial form is returned.
##
InstallMethod( QuadraticFormByMatrixOp, "for a ffe matrix and a field",
  [IsMatrix and IsFFECollColl, IsField and IsFinite],
  function( m, f )
    local el, n;
    n := NrRows(m);
    m := ImmutableMatrix(f, m);
    el := rec( matrix := m, basefield := f, type := "quadratic", vectorspace := FullRowSpace(f,n) );
    if IsZero(m) then
       el.type := "trivial";
       Objectify(NewType( TrivialFormFamily ,  IsFormRep),  el);
       return el;
    else
       el.matrix := ImmutableMatrix(f, Forms_RESET(m));
       Objectify(NewType( QuadraticFormFamily ,  IsFormRep),  el);
       return el;
    fi;
  end );

# 20/01/2016: small change here: added MutableCopyMat, since it is not logical
# that changing the matrix that was used for construction afterwards changes
# the form. See test_forms12.g
#############################################################################
#O QuadraticFormByMatrix( <m>, <f> ): constructor for users.
# <f> finite field, <m> matrix over <f>.
# checks whether <m> is over <f>.
# zero matrix is allowed.
##
InstallMethod( QuadraticFormByMatrix, "for a ffe matrix and a field",
  [IsMatrix and IsFFECollColl, IsField and IsFinite],
  function( m, f )
    local gf;
    gf := DefaultFieldOfMatrix(m);
    if not PrimitiveElement(gf) in f then
      Error("<m> is not a matrix over <f>");
    fi;
    return QuadraticFormByMatrixOp( m, f);
end );

#############################################################################
#O QuadraticFormByMatrix( <m> ): constructor for users. uses above constructor
# finite field <f> is determined from <m>
##
InstallMethod( QuadraticFormByMatrix, "for a ffe matrix ",
  [IsMatrix and IsFFECollColl],
  function( m )
  local f;
  f := DefaultFieldOfMatrix(m);
  return QuadraticFormByMatrixOp( m, f);
end );

# 20/01/2016: small change here: added MutableCopyMat, since it is not logical
# that changing the matrix that was used for construction afterwards changes
# the form. See test_forms12.g
#############################################################################
#O HermitianFormByMatrix( <m>, <f>): constructor  for users.
# <f> finite field of square order, <m> hermitian matrix over <f>.
# zero matrix is allowed, then the trivial form is returned.
##
InstallMethod( HermitianFormByMatrix, "for a ffe matrix and a field",
  [IsMatrix and IsFFECollColl, IsField and IsFinite],
  function( m, f )
    local el,gf,n;
    n := NrRows(m);
    gf := DefaultFieldOfMatrix(m);
    if not PrimitiveElement(gf) in f then
      Error("<m> is not a matrix over <f>");
    fi;
    if not IsInt(Sqrt(Size(f))) then
        Error("No hermitian form exists when the order of <f> is not a square" );
    fi;
    if FORMS_IsHermitianMatrix(m,f) then
       m := ImmutableMatrix(f, m);
       el := rec( matrix := m, basefield := f, type := "hermitian", vectorspace := FullRowSpace(f,n) );
       Objectify(NewType( HermitianFormFamily ,  IsFormRep),  el);
       return el;
    else
       Error("Given matrix does not define a hermitian form" );
    fi;
  end );

#############################################################################
#O UpperTriangleMatrixByPolynomialForForm( <pol>,<ff>,<int>,<list>),
# <ff>: fin field, <pol>: polynomial over <ff>, <n>: size of the matrix to
# construct; <list>: used variables in <pol>.
# not for users, creates a Gram matrix from a given polynomial to construct a
# quadratic form, and also an orthogonal bilinear form in odd characteristic.
##
InstallMethod( UpperTriangleMatrixByPolynomialForForm,
                      [IsPolynomial, IsField, IsInt, IsList],
  function(poly, gf, n, varlist)
    local vars, mat, i, j, vals, der;

   ## UpperTriangleMatrixByPolynomialForForm returns an upper triangular
   ## matrix in all cases. It can be used for quadratic forms (all chars)
   ## and symmetric bilinear forms (odd char). This is not a restriction,
   ## since in even char, we have no orthogonal bilinear forms, as by
   ## convention, symmetric bilinear forms in odd char, hermitian forms (all char)
   ## and quadratic forms (all char) can be constructed using polynomials.
   ##  Extermely important: the polynomial of a symmetric bilinear form is
   ##  NOT the polynomial of the associated quadratic form!

    mat := NullMat(n, n, gf);
    vars := ShallowCopy(varlist);
    if IsZero(poly) then
       return mat;
    fi;
    n := Length(varlist);
    vals := List([1..n], i -> 0);
    for i in [1..n] do
      vals[i] := 1;
      mat[i,i] := Value(poly, vars, vals);
      vals[i] := 0;
    od;
    for i in [1..n-1] do
      der := Derivative(poly, vars[i]);
      for j in [i+1..n] do
        vals[j] := 1;
        mat[i,j] := Value(der,vars,vals);
        vals[j] := 0;
      od;
    od;
    vars := ShallowCopy(varlist);
    if vars*mat*vars <> poly then
      Error( "<poly> should be a homogeneous polynomial over GF(q)");
    fi;
    return mat;
end);

#############################################################################
#O GramMatrixByPolynomialForHermitianForm( <pol>,<ff>,<int>,<list>),
# <ff>: fin field, <pol>: polynomial over <ff>, <n>: size of the matrix to
# construct; <list>: used variables in <pol>.
# not for users, creates a Gram matrix from a given polynomial to construct a
# hermitian form.
##
InstallMethod( GramMatrixByPolynomialForHermitianForm,
                     [IsPolynomial, IsField, IsInt, IsList],
  function(poly, gf, n, varlist)
    local vars, varst, q, t, a, polarity, i, j, vals, der;
    q := Size(gf);
    t := Sqrt(q);
    vars := ShallowCopy(varlist);
    varst := List(vars, x -> x^t);
    polarity := NullMat(n,n,gf);
    if IsZero(poly) then
      return polarity;
    fi;
    n := Length(varlist);
    vals := [];
    vals := List([1..n], i -> 0);
    for i in [1..n] do
      vals[i] := 1;
      a := Value(poly,vars,vals);
      if a <> a^t then
        Error( "<poly> does not generate a Hermitian matrix" );
      else
        polarity[i,i] := a;
      fi;
      vals[i] := 0;
    od;
    for i in [1..n-1] do
      der := Derivative(poly,vars[i]);
      for j in [i+1..n] do
          vals[j] := 1;
          polarity[i,j] := Value(der,vars,vals);
          polarity[j,i] := polarity[i,j]^t;
          vals[j] := 0;
      od;
    od;
    vars := ShallowCopy(varlist);
      #Check whether the polynomial has the right form.
    if vars*polarity*varst <> poly then
      Error( "<poly> does not generate a Hermitian matrix" );
    fi;
    return polarity;
  end );

#############################################################################
#O BilinearFormByPolynomial( <pol>, <ring>, <int>): constructor for users.
#  <ring>: polynomial ring over finite field, <po>: polyniomial in <ring>, <int>:
#  desired dimension of the vectorspace on which the form acts.
#  an error is returned by UpperTriangleMatrixByPolynomialForForm if the given
#  <pol> is not suitable to define a bilinear form
#  zero <pol > is allowed, then the trivial form is returned.
##
InstallMethod( BilinearFormByPolynomial, "for a polynomial over a field, and a dimension",
  [IsPolynomial, IsFiniteFieldPolynomialRing, IsInt],
  function( pol, pring, n )
    local mat, form, gf, vars, polarity,el;
    gf := CoefficientsRing( pring );
    vars := IndeterminatesOfPolynomialRing( pring );
    if IsZero(pol) then
      mat := NullMat(n, n, gf);
      mat := ImmutableMatrix(gf, mat);
      el := rec( matrix := mat, basefield := gf, type := "trivial" );
      Objectify(NewType( TrivialFormFamily ,  IsFormRep),  el);
      SetPolynomialOfForm(el, pol);
      return el;
    fi;
    if Characteristic(gf) = 2 then
       Error("No orthogonal form can be associated with a quadratic polynomial in even characteristic");
    else
       mat := UpperTriangleMatrixByPolynomialForForm(pol,gf,n,vars);
       polarity := (mat + TransposedMat(mat))/(One(gf)*2);
       form := FormByMatrix(polarity,gf,"orthogonal");
       SetPolynomialOfForm(form, pol);
       return form;
    fi;
  end );

#############################################################################
#O BilinearFormByPolynomial( <pol>, <ring> ): constructor for users.
#  <ring>: polynomial ring over finite field, <po>: polyniomial in <ring>
#  dimension is the length of variables in <ring>.
##
InstallMethod( BilinearFormByPolynomial,  "no dimension",
  [IsPolynomial, IsFiniteFieldPolynomialRing],
  function( pol, pring )
    local n;
    n := Length(IndeterminatesOfPolynomialRing( pring ));
    return BilinearFormByPolynomial( pol, pring, n);
  end );

#############################################################################
#O HermitianFormByPolynomial( <pol>, <ring>, <int>): constructor for users.
#  <ring>: polynomial ring over finite field, <po>: polyniomial in <ring>, <int>:
#  desired dimension of the vectorspace on which the form acts.
#  an error is returned by GramMatrixByPolynomialForHermitianForm if the given
#  <pol> is not suitable to define a hermitian form
#  zero <pol > is allowed, then the trivial form is returned.
##
InstallMethod( HermitianFormByPolynomial, "for a polynomial over a field, and a dimension",
  [IsPolynomial, IsFiniteFieldPolynomialRing, IsInt],
  function(pol, pring, n)
    local mat, form, gf, vars, el;
    gf := CoefficientsRing( pring );
    if not IsInt(Sqrt(Size(gf))) then
        Error("the order of the underlying field of <r> must be a square" );
    fi;
    if IsZero(pol) then
      mat := NullMat(n, n, gf);
      mat := ImmutableMatrix(gf, mat);
      el := rec( matrix := mat, basefield := gf, type := "trivial" );
      Objectify(NewType( TrivialFormFamily ,  IsFormRep),  el);
      SetPolynomialOfForm(el, pol);
      return el;
    fi;
    vars := IndeterminatesOfPolynomialRing( pring );
    mat := GramMatrixByPolynomialForHermitianForm(pol,gf,n,vars);
    form := FormByMatrix(mat,gf,"hermitian");
    SetPolynomialOfForm(form, pol);
    return form;
  end );

#############################################################################
#O HermitianFormByPolynomial( <pol>, <ring> ): constructor for users.
#  <ring>: polynomial ring over finite field, <pol>: polyniomial in <ring>
#  dimension is the length of variables in <ring>.
##
InstallMethod( HermitianFormByPolynomial,  "no dimension",
  [IsPolynomial, IsFiniteFieldPolynomialRing],
  function( pol, pring )
    local n;
    n := Length(IndeterminatesOfPolynomialRing( pring ));
    return HermitianFormByPolynomial( pol, pring, n);
  end );

#############################################################################
#O QuadraticFormByPolynomial( <pol>, <ring>, <int>): constructor for users.
#  <ring>: polynomial ring over finite field, <po>: polyniomial in <ring>, <int>:
#  desired dimension of the vectorspace on which the form acts.
#  an error is returned by UpperTriangleMatrixByPolynomialForForm if the given
#  <pol> is not suitable to define a quadratic form
#  zero <pol > is allowed, then the trivial form is returned.
##
InstallMethod( QuadraticFormByPolynomial,
  "for a polynomial over a field, and a dimension",
  [IsPolynomial, IsFiniteFieldPolynomialRing, IsInt],
  function(pol, pring, n)
    local mat, form, gf, vars;
    gf := CoefficientsRing( pring );
    vars := IndeterminatesOfPolynomialRing( pring );
    mat := UpperTriangleMatrixByPolynomialForForm(pol,gf,n,vars);
    form := FormByMatrix(mat,gf,"quadratic");
    SetPolynomialOfForm(form, pol);
    return form;
  end );

#############################################################################
#O QuadraticFormByPolynomial( <pol>, <ring> ): constructor for users.
#  <ring>: polynomial ring over finite field, <pol>: polyniomial in <ring>
#  dimension is the length of variables in <ring>.
##
InstallMethod( QuadraticFormByPolynomial,  "no dimension",
  [IsPolynomial, IsFiniteFieldPolynomialRing],
  function( pol, pring )
    local n;
    n := Length(IndeterminatesOfPolynomialRing( pring ));
    return QuadraticFormByPolynomial( pol, pring, n);
  end );

#############################################################################
# Form-by-form methods (for users):
# see Package documentation for the interpretation of these functions.
#############################################################################

#############################################################################
#O BilinearFormByQuadraticForm( <form> ): constructor for users.
#  <form>: quadratic form.
##
InstallMethod( BilinearFormByQuadraticForm, [IsQuadraticForm],
  function(f)
    ## This method constructs a bilinear form from a quadratic form.
    ## very important: we use the relation 2Q(v) = f(v,v)
    local m, gf;
    m := f!.matrix;
    gf := f!.basefield;
    if Characteristic(gf) = 2 then
      Error( "No orthogonal bilinear form can be associated with a quadratic form in even characteristic");
    else
      return BilinearFormByMatrix((m+TransposedMat(m))/(One(gf)*2),gf);
    fi;
  end );

#############################################################################
#O BilinearFormByQuadraticForm( <form> ): constructor for users.
#  <form>: bilinear form.
##
InstallMethod( QuadraticFormByBilinearForm, [IsBilinearForm],
  function(f)
    ## This method constructs a quadratic form from a bilinear form.
    ## very important: we use the relation 2Q(v) = f(v,v)
    local m, gf;
    m := f!.matrix;
    gf := f!.basefield;
    if Characteristic(gf) = 2 then
       Error( "No quadratic form can be associated to a symmetric form in even characteristic");
    elif IsAlternatingForm(f) then
      Error( "No quadratic form can be associated properly to an alternating form" );
    else
      return QuadraticFormByMatrix(m,gf);
    fi;
  end );

#############################################################################
#  Attributes ... (for users).
#  see package documentation for more information.
##
InstallMethod( PolynomialOfForm, "for a trivial form",
  [IsTrivialForm],
  function(f)
    ## returns zero polynomial of the corresponding polynomialring.
    local gf, d, r;
    gf := f!.basefield;
    d := NrRows(f!.matrix);
    r := PolynomialRing(gf,d);
    return Zero(r);
  end );

InstallMethod( PolynomialOfForm, "for a quadratic form",
  [IsQuadraticForm],
  function(f)
    ## This method finds the polynomial associated to
    ## the Gram matrix of f.
    local m, gf, d, r, indets, poly;
    m := f!.matrix;
    gf := f!.basefield;
    d := NrRows(m);
    r := PolynomialRing(gf,d);
    indets := IndeterminatesOfPolynomialRing(r);
    poly := indets * m * indets;
    return poly;
  end );

InstallMethod( PolynomialOfForm, "for a hermitian form",
  [IsHermitianForm],
  function(f)
    ## This method finds the polynomial associated to
    ## the Gram matrix of f.
    local m, gf, d, r, indets, poly, q;
    gf := f!.basefield;
    q := Sqrt(Size(gf));
    m := f!.matrix;
    d := NrRows(m);
    r := PolynomialRing(gf,d);
    indets := IndeterminatesOfPolynomialRing(r);
    poly := indets * m * List(indets,t->t^q);
    return poly;
  end );

InstallMethod( PolynomialOfForm, "for a bilinear form",
  [IsBilinearForm],
  function(f)
    ## This method finds the polynomial associated to
    ## the Gram matrix of f. For a symplectic form, we
    ## get the 0 polynomial.
    ## Very important: for an orthogonal form, the polynomial is
    ## NOT the polynomial of the associated quadratic form.
    local m, gf, d, r, indets, poly;
    m := f!.matrix;
    gf := f!.basefield;
    d := NrRows(m);
    if Characteristic(gf) = 2 then
       Error( "No polynomial can be (naturally) associated to a bilinear form in even characteristic");
    fi;
    r := PolynomialRing(gf,d);
    indets := IndeterminatesOfPolynomialRing(r);
    poly := indets * m * indets;
    return poly;
  end );

InstallOtherMethod( BaseField, "for a form", [IsForm],
  function( f )
    return f!.basefield;
  end );

InstallMethod( GramMatrix, "for a form", [IsForm],
  function( f )
    return f!.matrix;
  end );

InstallMethod( CompanionAutomorphism, [ IsSesquilinearForm ],
  function( form )
    local field, aut;
    field := BaseField( form );
    if IsHermitianForm( form ) then
       aut := FrobeniusAutomorphism( field );
       aut := aut ^ (Order(aut)/2);
    else
       aut := IdentityMapping( field );
    fi;
    return aut;
  end );

InstallMethod( AssociatedBilinearForm, "for a quadratic form",
  [ IsQuadraticForm ],
  function( form )
    local gram, f;
    gram := form!.matrix;
    f := form!.basefield;
    return BilinearFormByMatrix( gram + TransposedMat(gram), f);
  end );

##next two operations are not documented. Could be in next versions.

InstallMethod( RadicalOfFormBaseMat, [IsSesquilinearForm],
  function( f )
    local m, gf, d;
    if not IsReflexiveForm( f ) then
       Error( "Form must be reflexive");
    fi;
    m := f!.matrix;
    gf := f!.basefield;
    d := NrRows(m);
    return NullspaceMat( m );
  end );

InstallMethod( RadicalOfFormBaseMat, [IsQuadraticForm],
  function( f )
    local m, null, gf, d;
    m := f!.matrix;
    m := m + TransposedMat(m);
    gf := f!.basefield;
    d := NrRows(m);
    null := NullspaceMat( m );
    if Characteristic(gf) = 2 then
      null := Filtered(SubspaceNC(gf^d,null), x -> IsZero(x^f)); #find vectors vanishing under f
    fi;
    null := Filtered(null,x-> not IsZero(x));
    return null;
  end );

InstallMethod( RadicalOfForm, "for a sesquilinear form",
  [IsSesquilinearForm],
  function( f )
    local m, null, gf, d;
    if not IsReflexiveForm( f ) then
       Error( "Form must be reflexive");
    fi;
    m := f!.matrix;
    gf := f!.basefield;
    d := NrRows(m);
    null := NullspaceMat( m );
    return Subspace( gf^d, null, "basis" );
  end );

InstallMethod( RadicalOfForm, "for a quadratic form",
  [IsQuadraticForm],
  function( f )
    local m, null, gf, d;
    m := f!.matrix;
    m := m + TransposedMat(m);
    gf := f!.basefield;
    d := NrRows(m);
    null := NullspaceMat( m );
    if Characteristic(gf) = 2 then
      null := Filtered(SubspaceNC(gf^d,null), x -> IsZero(x^f)); #find vectors vanishing under f
    fi;
    null := Filtered(null,x-> not IsZero(x));
    return Subspace( gf^d, null );
  end );

InstallMethod( RadicalOfForm, "for a trivial form",
  [IsTrivialForm],
  function( f )
    local m, null, gf, d;
    m := f!.matrix;
    gf := f!.basefield;
    d := Size(m);
    null := NullspaceMat( m );
    return Subspace( gf^d, null, "basis" );
  end );

InstallMethod( DiscriminantOfForm, [ IsQuadraticForm ],
 function( f )
   local m, gf, d, det, squares, primroot;
   m := f!.matrix;
   gf := f!.basefield;
   d := Size(m);
   if IsOddInt(d) then
      Error( "Quadratic form must be defined by a matrix of even dimension");
   fi;
   det := Determinant( m + TransposedMat(m) );
   if IsZero( det ) then
      Error( "Form must be nondegenerate" );
   fi;
   if IsOddInt(Size(gf)) then
      primroot := PrimitiveRoot( gf );
      squares := AsList( Group( primroot^2 ) );
      if det in squares then
         return "square";
      else
         return "nonsquare";
      fi;
   else
      return "square";
   fi;
 end );

InstallMethod( DiscriminantOfForm, [ IsSesquilinearForm ],
 function( f )
   local m, gf, d, det, squares, primroot;
   m := f!.matrix;
   gf := f!.basefield;
   d := Size(m);
   if IsOddInt(d) then
      Error( "Sesquilinear form must be defined by a matrix of even dimension");
   fi;
   det := Determinant( m );
   if IsZero( det ) then
      Error( "Form must be nondegenerate" );
   fi;
   if IsOddInt(Size(gf)) then
      primroot := PrimitiveRoot( gf );
      squares := AsList( Group( primroot^2 ) );
      if det in squares then
         return "square";
      else
         return "nonsquare";
      fi;
   else
      return "square";
   fi;
 end );

InstallMethod( DiscriminantOfForm, [ IsHermitianForm ],
 function( f )
   Error( "The discriminant of a hermitian form is not defined");
 end );

InstallMethod( DiscriminantOfForm, [ IsTrivialForm ],
 function( f )
   Error( "<form> must be sesquilinear or quadratic");
 end );

####
#  ... and Properties (for users).
#  see package documentation for more information.
####

InstallMethod( IsDegenerateForm, [IsSesquilinearForm],
  function( f )
    return not IsTrivial( RadicalOfForm( f ) );
  end );

InstallMethod( IsSingularForm, [IsQuadraticForm],
  function( f )
    return not IsTrivial( RadicalOfForm( f ) );
  end );

InstallMethod( IsSingularForm, [IsTrivialForm], #new in 1.2.1
  function( f )
  return true;
  end );

InstallMethod( IsDegenerateForm, [IsQuadraticForm],
  function( f )
    Info(InfoWarning,1,"Testing degeneracy of the *associated bilinear form*");
    return not IsTrivial( RadicalOfForm( AssociatedBilinearForm( f ) ) );
  end );

InstallMethod( IsDegenerateForm, [IsTrivialForm],
  function( f )
    return true;
  end );

InstallMethod( IsReflexiveForm, [IsBilinearForm],
  function( f )
    ## simply check if form is symmetric or alternating
    local m;
    m := f!.matrix;
    return m = TransposedMat(m) or m = -TransposedMat(m);
  end );

InstallMethod( IsReflexiveForm, [IsHermitianForm],
  function( f )
    return FORMS_IsHermitianMatrix( f!.matrix, f!.basefield);
  end );

InstallMethod( IsReflexiveForm, [IsTrivialForm], #new in 1.2.1
 function( f )
   return true;
end );

InstallMethod( IsReflexiveForm, [IsQuadraticForm],
  function( f )
    Error( "<form> must be sesquilinear" );
end );

InstallMethod( IsAlternatingForm, [IsBilinearForm],
  function( f )
    return FORMS_IsSymplecticMatrix( f!.matrix, f!.basefield );
 end );

InstallMethod( IsAlternatingForm, [IsHermitianForm],
  function( f )
    return false;
 end );

InstallMethod( IsAlternatingForm, [IsTrivialForm], #new in 1.2.1
 function( f )
   return true;
end );

InstallMethod( IsAlternatingForm, [IsQuadraticForm],
  function( f )
    Error( "<form> must be sesquilinear" );
end );

InstallMethod( IsSymmetricForm, [IsBilinearForm],
  function( f )
    return FORMS_IsSymmetricMatrix( f!.matrix );
 end );

InstallMethod( IsSymmetricForm, [IsHermitianForm],
  function( f )
    return false;
 end );

InstallMethod( IsSymmetricForm, [IsTrivialForm], #new in 1.2.1
 function( f )
   return true;
end );

InstallMethod( IsSymmetricForm, [IsQuadraticForm],
  function( f )
    Error( "<form> must be sesquilinear" );
end );

InstallMethod( IsSymplecticForm, "for sesquilinear forms",
  [IsSesquilinearForm and IsFormRep],
  function( f )
    local string;
    string := f!.type;
    if string = "symplectic" then
       if IsAlternatingForm( f ) then
          return true;
       else
          Error( "<form> was incorrectly specified" );
       fi;
    else
       return false;
    fi;
  end );

InstallMethod( IsSymplecticForm, [IsTrivialForm], #new in 1.2.1
 function( f )
   return true;
end );


InstallMethod( IsOrthogonalForm, "for sesquilinear forms",
  [IsSesquilinearForm and IsFormRep],
  function( f )
    local string;
    string := f!.type;
    if string = "orthogonal" then
       if IsSymmetricForm( f ) then
          return true;
       else
          Error( "Form was incorrectly specified" );
       fi;
    else
       return false;
    fi;
  end );

InstallMethod( IsOrthogonalForm, [IsTrivialForm], #new in 1.2.1
 function( f )
   return false;
end );

InstallMethod( IsOrthogonalForm, [IsQuadraticForm],
  function( f )
    Error( "<form> must be sesquilinear" );
end );

InstallMethod( IsPseudoForm, "for sesquilinear forms",
  [IsSesquilinearForm and IsFormRep],
  function( f )
    local string;
    string := f!.type;
    if string = "pseudo" then
       if (IsSymmetricForm( f ) and (not IsAlternatingForm(f))) then
          return true;
       else
          Error( "Form was incorrectly specified" );
       fi;
    else
       return false;
    fi;
  end );

InstallMethod( IsPseudoForm, [IsTrivialForm], #new in 1.2.1
 function( f )
   return false;
end );

InstallMethod( IsPseudoForm, [IsQuadraticForm],
  function( f )
    Error( "<form> must be sesquilinear" );
end );

##
#############################################################################

#############################################################################
# Base change methods (for users):
#############################################################################

#############################################################################
#O BaseChangeToCanonical( <form> ) <form>: trivial form.
#  returns warning and identity matrix.
##
InstallMethod( BaseChangeToCanonical, "for a trivial form",
  [IsTrivialForm],
  function(f)
    local b,m,n;
    m := f!.matrix;
    n := NrRows(m);
    b := IdentityMat(n,m);
    Info(InfoWarning,1,"<form> is trivial, trivial base change is returned");
    return b;
end );

#############################################################################
#O BaseChangeToCanonical( <form> ) <form>: sesquilinear form.
#  this function checks the type of the given form and calls the appropriate
#  main base change operation.
##
InstallMethod( BaseChangeToCanonical, "for a sesquilinear form",
  [IsSesquilinearForm and IsFormRep],
  function(f)
    local string,m,gf,b;
    string := f!.type;
    m := f!.matrix;
    gf := f!.basefield;
    if IsOrthogonalForm(f) then
      b := BaseChangeOrthogonalBilinear(m, gf);
      SetWittIndex(f, (b[2]+b[3]-1) / 2);
      SetIsEllipticForm(f,b[3]=0);
      SetIsParabolicForm(f,b[3]=1);
      SetIsHyperbolicForm(f,b[3]=2);
      return b[1];
    elif IsSymplecticForm(f) then
      b := BaseChangeSymplectic(m, gf);
      SetWittIndex(f, b[2]/2 );
      return b[1];
    elif string = "hermitian" then
      b := BaseChangeHermitian(m, gf);
      SetWittIndex(f, Int( (b[2]+1)/2 ));
      return b[1];
    elif string = "pseudo" then
      Error("BaseChangeToCanonical not yet implemented for pseudo forms");
    fi;
end );

#############################################################################
#O BaseChangeToCanonical( <form> ) <form>: quadratic form.
#  this function checks the type of the given form and calls the appropriate
#  main base change operation.
##
InstallMethod( BaseChangeToCanonical, "for a quadratic form",
  [IsQuadraticForm and IsFormRep],
  function(f)
    local m, gf, b, polarity;
    m := f!.matrix;
    gf := f!.basefield;
    if IsOddInt(Size(gf)) then
      polarity := (m+TransposedMat(m))/2;
      b := BaseChangeOrthogonalBilinear(polarity, gf);
    else
      b := BaseChangeOrthogonalQuadratic(m, gf);
    fi;
    SetWittIndex(f, (b[2]+b[3]-1) / 2);
    SetIsEllipticForm(f,b[3]=0);
    SetIsParabolicForm(f,b[3]=1);
    SetIsHyperbolicForm(f,b[3]=2);
    return b[1];
  end );

#############################################################################
# Overloading: Frobenius Automorphisms
#############################################################################

InstallOtherMethod( \^, "for a FFE vector and a Frobenius automorphism",
  [ IsVector and IsFFECollection, IsFrobeniusAutomorphism ],
  function( v, f )
    return List(v,x->x^f);
  end );

InstallOtherMethod( \^, "for a FFE vector and a trivial Frobenius automorphism",
  [ IsVector and IsFFECollection, IsMapping and IsOne ],
  function( v, f )
    return v;
  end );

InstallOtherMethod( \^,
  "for a compressed GF2 vector and a Frobenius automorphism",
  [ IsVector and IsFFECollection and IsGF2VectorRep, IsFrobeniusAutomorphism ],
  function( v, f )
    local w;
    w := List(v,x->x^f);
    ConvertToVectorRepNC(w,2);
    return w;
  end );

InstallOtherMethod( \^,
  "for a compressed GF2 vector and a trivial Frobenius automorphism",
  [ IsVector and IsFFECollection and IsGF2VectorRep, IsMapping and IsOne ],
  function( v, f )
    return v;
  end );

InstallOtherMethod( \^,
  "for a compressed 8bit vector and a Frobenius automorphism",
  [ IsVector and IsFFECollection and Is8BitVectorRep, IsFrobeniusAutomorphism ],
  function( v, f )
    local w;
    w := List(v,x->x^f);
    ConvertToVectorRepNC(w,Q_VEC8BIT(v));
    return w;
  end );

InstallOtherMethod( \^,
  "for a compressed 8bit vector and a trivial Frobenius automorphism",
  [ IsVector and IsFFECollection and Is8BitVectorRep, IsMapping and IsOne ],
  function( v, f )
    return v;
  end );

InstallOtherMethod( \^, "for a FFE matrix and a Frobenius automorphism",
  [ IsMatrix and IsFFECollColl, IsFrobeniusAutomorphism ],
  function( m, f )
    return List(m,v->List(v,x->x^f));
  end );

InstallOtherMethod( \^, "for a FFE matrix and a trivial Frobenius automorphism",
  [ IsMatrix and IsFFECollColl, IsMapping and IsOne ],
  function( m, f )
    return m;
  end );

InstallOtherMethod( \^,
  "for a compressed GF2 matrix and a Frobenius automorphism",
  [ IsMatrix and IsFFECollColl and IsGF2MatrixRep, IsFrobeniusAutomorphism ],
  function( m, f )
    local w,l,i;
    l := [];
    for i in [1..NrRows(m)] do
        w := List(m[i],x->x^f);
        ConvertToVectorRepNC(w,2);
        Add(l,w);
    od;
    ConvertToMatrixRepNC(l,2);
    return l;
  end );

InstallOtherMethod( \^,
  "for a compressed GF2 matrix and a trivial Frobenius automorphism",
  [ IsMatrix and IsFFECollColl and IsGF2MatrixRep, IsMapping and IsOne ],
  function( m, f )
    return m;
  end );

InstallOtherMethod( \^,
  "for a compressed 8bit matrix and a Frobenius automorphism",
  [ IsMatrix and IsFFECollColl and Is8BitMatrixRep, IsFrobeniusAutomorphism ],
  function( m, f )
    local w,l,i,q;
    l := [];
    q := Q_VEC8BIT(m[1]);
    for i in [1..NrRows(m)] do
        w := List(m[i],x->x^f);
        ConvertToVectorRepNC(w,q);
        Add(l,w);
    od;
    ConvertToMatrixRepNC(l,q);
    return l;
  end );

InstallOtherMethod( \^,
  "for a compressed 8bit matrix and a trivial Frobenius automorphism",
  [ IsMatrix and IsFFECollColl and Is8BitMatrixRep, IsMapping and IsOne ],
  function( m, f )
    return m;
  end );

#############################################################################
# Overloading: Forms
#############################################################################

InstallOtherMethod( \^, "for a pair of FFE vectors and a sesquilinear form",
  [ IsVectorList and IsFFECollColl, IsBilinearForm ],
  function( pair, f )
    if Size(pair) <> 2 then
       Error("The first argument must be a pair of vectors");
    fi;
    return pair[1] * f!.matrix * pair[2];
  end );

InstallOtherMethod( \^, "for a pair of FFE matrices and a sesquilinear form",
  [ IsFFECollCollColl, IsBilinearForm ],
  function( pair, f )
    if Size(pair) <> 2 then
       Error("The first argument must be a pair of vectors");
    fi;
    return pair[1] * f!.matrix * TransposedMat(pair[2]);
  end );

InstallOtherMethod( \^, "for a pair of FFE vectors and an hermitian form",
  [ IsVectorList and IsFFECollColl, IsHermitianForm ],
  function( pair, f )
    local frob,hh,bf;
    if Size(pair) <> 2 then
       Error("The first argument must be a pair of vectors");
    fi;
    bf := f!.basefield;
    hh := DegreeOverPrimeField(bf) / 2;
    frob := FrobeniusAutomorphism(bf)^hh;
    return pair[1] * f!.matrix * (pair[2]^frob);
  end );

InstallOtherMethod( \^, "for a pair of FFE matrices and an hermitian form",
  [ IsFFECollCollColl, IsHermitianForm ],
  function( pair, f )
    local frob,hh,bf;
    if Size(pair) <> 2 then
       Error("The first argument must be a pair of vectors");
    fi;
    bf := f!.basefield;
    hh := DegreeOverPrimeField(bf) / 2;
    frob := FrobeniusAutomorphism(bf)^hh;
    return pair[1] * f!.matrix * (TransposedMat(pair[2])^frob);
  end );

InstallOtherMethod( \^, "for a pair of FFE matrices and a trivial form", #new in 1.2.1
  [ IsVectorList and IsFFECollColl, IsTrivialForm ],
  function( pair, f )
    if Size(pair) <> 2 then
       Error("The first argument must be a pair of vectors or a vector");
    fi;
    return Zero(BaseField(f));
  end );

InstallOtherMethod( \^, "for a FFE vector and a quadratic form",
  [ IsVector and IsFFECollection, IsQuadraticForm ],
  function( v, f )
    return v * f!.matrix * v;
  end );

InstallOtherMethod( \^, "for a FFE matrix and a quadratic form",
  [ IsMatrix and IsFFECollColl, IsQuadraticForm ],
  function( m, f )
    return m * f!.matrix * TransposedMat(m);
  end );

InstallOtherMethod( \^, "for a FFE vector and a quadratic form", #new in 1.2.1
  [ IsVector and IsFFECollection, IsTrivialForm ],
  function( m, f )
    return Zero(BaseField(f));
  end );

#############################################################################
# Viewing methods:
#############################################################################
#
InstallMethod( ViewObj, [ IsTrivialForm ],
  function( f )
    Print("< trivial form >");
  end );

InstallMethod( PrintObj, [ IsTrivialForm ],
  function( f )
    Print("Trivial form\n");
    Print("Gram Matrix:\n",f!.matrix,"\n");
  end );

InstallMethod( Display, [ IsTrivialForm ],
  function( f )
    Print("Trivial form\n");
    Print("Gram Matrix:\n");
    Display(f!.matrix);
  end);


InstallMethod( ViewObj, [ IsHermitianForm ],
  function( f )
    if HasIsDegenerateForm(f) then
      if IsDegenerateForm(f) then
        Print(" < degenerate hermitian form >");
      else
        Print(" < non-degenerate hermitian form >");
      fi;
    else
      Print("< hermitian form >");
    fi;
  end );

InstallMethod( PrintObj, [ IsHermitianForm ],
  function( f )
    Print("Hermitian form\n");
    Print("Gram Matrix:\n",f!.matrix,"\n");
    if HasPolynomialOfForm( f ) then
       Print("Polynomial: ", PolynomialOfForm, "\n");
    fi;
    if HasWittIndex( f ) then
       Print("Witt Index: ", WittIndex(f), "\n");
    fi;
  end );

InstallMethod( Display, [ IsHermitianForm ],
  function( f )
    Print("Hermitian form\n");
    Print("Gram Matrix:\n");
    Display(f!.matrix);
    if HasPolynomialOfForm( f ) then
       Print("Polynomial: ");
       Display(PolynomialOfForm(f));
       Print("\n");
    fi;
    if HasWittIndex( f ) then
       Print("Witt Index: ", WittIndex(f), "\n");
    fi;
  end);

InstallMethod( ViewObj, [ IsQuadraticForm ],
  function( f )
    local string;
    string := ["< "];
    if HasIsSingularForm(f) then
      if IsSingularForm(f) then
        Add(string,"singular ");
      else
        Add(string,"non-singular ");
      fi;
    fi;
    if HasIsEllipticForm( f ) or HasIsHyperbolicForm( f ) or
      HasIsParabolicForm( f ) then
      if IsEllipticForm( f ) then
          Add(string,"elliptic ");
      elif IsHyperbolicForm( f ) then
          Add(string,"hyperbolic ");
      elif IsParabolicForm( f)  then
          Add(string,"parabolic ");
      fi;
    fi;
    Add(string,"quadratic form >");
    string := Concatenation(string);
    Print(string);
  end );

InstallMethod( PrintObj, [ IsQuadraticForm ],
  function( f )
    local string;
    string := [];
    if HasIsSingularForm(f) then
       if IsSingularForm(f) then
          Add(string,"Singular ");
       else
          Add(string,"Non-singular ");
       fi;
    fi;
    if HasIsEllipticForm( f ) or
       HasIsHyperbolicForm( f ) or
       HasIsParabolicForm( f ) then

       if IsEllipticForm( f ) then
          Add(string,"Elliptic ");
       elif IsHyperbolicForm( f ) then
          Add(string,"Hyperbolic ");
       elif IsParabolicForm( f)  then
          Add(string,"Parabolic ");
       fi;
     fi;
     Add(string,"Quadratic form\n");
     string := Concatenation(string[1],LowercaseString(Concatenation(string{[2..Length(string)]})));
     Print(string);
     Print("Gram Matrix:\n",f!.matrix,"\n");
     if HasPolynomialOfForm( f ) then
        Print("Polynomial: ", PolynomialOfForm(f), "\n");
     fi;
     if HasWittIndex( f ) then
        Print("Witt Index: ", WittIndex(f), "\n");
     fi;
  end );

InstallMethod( Display,  [ IsQuadraticForm ],
  function( f )
    local string;
    string := [];
    if HasIsSingularForm(f) then
       if IsSingularForm(f) then
          Add(string,"Singular ");
       else
          Add(string,"Non-singular ");
       fi;
    fi;
    if HasIsEllipticForm( f ) or
       HasIsHyperbolicForm( f ) or
       HasIsParabolicForm( f ) then

       if IsEllipticForm( f ) then
          Add(string,"Elliptic ");
       elif IsHyperbolicForm( f ) then
          Add(string,"Hyperbolic ");
       elif IsParabolicForm( f)  then
          Add(string,"Parabolic ");
       fi;
    fi;
    Add(string,"Quadratic form\n");
    string := Concatenation(string[1],LowercaseString(Concatenation(string{[2..Length(string)]})));
    Print(string);
    Print("Gram Matrix:\n");
    Display(f!.matrix);
    if HasPolynomialOfForm( f ) then
       Print("Polynomial: ");
       Display(PolynomialOfForm(f));
       Print("\n");
    fi;
    if HasWittIndex( f ) then
       Print("Witt Index: ", WittIndex(f), "\n");
    fi;
  end );

InstallMethod( ViewObj, [ IsBilinearForm ],
  function( f )
    local string;
    string := ["< "];
    if HasIsDegenerateForm(f) then
      if IsDegenerateForm(f) then
        Add(string,"degenerate ");
      else
        Add(string,"non-degenerate ");
      fi;
    fi;
    if HasIsOrthogonalForm(f) then
      if HasIsEllipticForm( f ) or
         HasIsHyperbolicForm( f ) or
         HasIsParabolicForm( f ) then

         if IsEllipticForm( f ) then
            Add(string,"elliptic bilinear ");
         elif IsHyperbolicForm( f ) then
            Add(string,"hyperbolic bilinear ");
         elif IsParabolicForm( f)  then
            Add(string,"parabolic bilinear ");
         fi;
      elif IsOrthogonalForm(f) then
         Add(string,"orthogonal ");
      fi;
    elif HasIsSymplecticForm(f) then
      if IsSymplecticForm(f) then
         Add(string,"symplectic ");
      fi;
    elif HasIsPseudoForm(f) then
      if IsPseudoForm(f) then
         Add(string,"pseudo ");
      fi;
    else
      Add(string,"bilinear ");
    fi;
    Add(string,"form >");
    string := Concatenation(string);
    Print(string);
  end );

InstallMethod( PrintObj, [ IsBilinearForm ],
  function( f )
    local string;
    string := [];
    if HasIsDegenerateForm(f) then
      if IsDegenerateForm(f) then
         Add(string,"Degenerate ");
      else
         Add(string,"Non-degenerate ");
      fi;
    fi;
    if HasIsOrthogonalForm(f) then
       if HasIsEllipticForm( f ) or
          HasIsHyperbolicForm( f ) or
          HasIsParabolicForm( f ) then
        if IsEllipticForm( f ) then
           Add(string,"Elliptic bilinear ");
        elif IsHyperbolicForm( f ) then
           Add(string,"Hyperbolic bilinear ");
        elif IsParabolicForm( f)  then
           Add(string,"Parabolic bilinear ");
        fi;
      elif IsOrthogonalForm(f) then
           Add(string,"Orthogonal ");
      fi;
    elif HasIsSymplecticForm(f) then
      if IsSymplecticForm(f) then
         Add(string,"Symplectic ");
      fi;
    elif HasIsPseudoForm(f) then
      if IsPseudoForm(f) then
         Add(string,"Pseudo ");
      fi;
    else
         Add(string,"Bilinear ");
    fi;
    Add(string,"form\n");
    string := Concatenation(string[1],LowercaseString(Concatenation(string{[2..Length(string)]})));
    Print(string);
    Print("Gram Matrix:\n",f!.matrix,"\n");
    if HasPolynomialOfForm( f ) then
       Print("Polynomial: ", PolynomialOfForm(f), "\n");
    fi;
    if HasWittIndex( f ) then
       Print("Witt Index: ", WittIndex(f), "\n");
    fi;
  end );

InstallMethod( Display, [ IsBilinearForm ],
  function( f )
    local string;
    string := [];
    if HasIsDegenerateForm(f) then
       if IsDegenerateForm(f) then
          Add(string,"Degenerate ");
       else
          Add(string,"Non-degenerate ");
       fi;
    fi;
    if HasIsOrthogonalForm(f) then
       if HasIsEllipticForm( f ) or
          HasIsHyperbolicForm( f ) or
          HasIsParabolicForm( f ) then

          if IsEllipticForm( f ) then
             Add(string,"Elliptic bilinear ");
          elif IsHyperbolicForm( f ) then
             Add(string,"Hyperbolic bilinear ");
          elif IsParabolicForm( f)  then
             Add(string,"Parabolic bilinear ");
          fi;
       elif IsOrthogonalForm(f) then
            Add(string,"Orthogonal ");
       fi;
    elif HasIsSymplecticForm(f) then
       if IsSymplecticForm(f) then
          Add(string,"Symplectic ");
        fi;
    elif HasIsPseudoForm(f) then
       if IsPseudoForm(f) then
          Add(string,"Pseudo ");
       fi;
    else
       Add(string,"Bilinear ");
    fi;
    Add(string,"form\n");
    string := Concatenation(string[1],LowercaseString(Concatenation(string{[2..Length(string)]})));
    Print(string);
    Print("Gram Matrix:\n");
    Display(f!.matrix);
    if HasPolynomialOfForm( f ) then
       Print("Polynomial: ");
       Display(PolynomialOfForm(f));
       Print("\n");
    fi;
    if HasWittIndex( f ) then
       Print("Witt Index: ", WittIndex(f), "\n");
    fi;
  end );

# end of Viewing methods.
#############################################################################

#############################################################################
# Functions to support Base Change operations (not for the user):
##
if IsBound(SwapMatrixColumns) and IsBound(SwapMatrixRows) then
  # For GAP >= 4.12
  BindGlobal("Forms_SwapRows", SwapMatrixRows);
  BindGlobal("Forms_SwapCols", SwapMatrixColumns);
  BindGlobal("Forms_AddRows", AddMatrixRows);
  BindGlobal("Forms_AddCols", AddMatrixColumns);
  BindGlobal("Forms_MultRow", MultMatrixRow);
  BindGlobal("Forms_MultCol", MultMatrixColumn);
else
  # For GAP <= 4.11
  BindGlobal("Forms_SwapRows", function(mat, i, j)
    mat{[i,j]} := mat{[j,i]};
  end);

  BindGlobal("Forms_SwapCols", function(mat, i, j)
    local row;
    for row in mat do
      row{[i,j]} := row{[j,i]};
    od;
  end);

  BindGlobal("Forms_AddRows", function(mat, i, j, scalar)
    mat[i] := mat[i] + mat[j] * scalar;
  end);

  BindGlobal("Forms_AddCols", function(mat, i, j, scalar)
    local row;
    for row in mat do
      row[i] := row[i] + row[j] * scalar;
    od;
  end);

  BindGlobal("Forms_MultRow", function(mat, i, scalar)
    mat[i] := mat[i] * scalar;
  end);

  BindGlobal("Forms_MultCol", function(mat, i, scalar)
    local row;
    for row in mat do
      row[i] := row[i] * scalar;
    od;
  end);

fi;

# express v as sum of two squares of elements in gf
InstallGlobalFunction(Forms_SUM_OF_SQUARES,
  function(v,gf)
    local dummy,i,v1,v2, primroot;
    primroot := PrimitiveRoot(gf);
    i := 0;
    repeat
      dummy := LogFFE(v - primroot^(2*i), primroot);
      if dummy mod 2 = 0 then
        v1 := primroot^i;
        v2 := primroot^(dummy/2);
        break;
      else
        i := i + 1;
      fi;
    until false;
    return [v1,v2];
  end );

# Apply a 2x2 transformation matrix to rows 'p1' and 'p2' of the matrix 'D'.
# The matrix 'D' may be changed in-place by this. For efficiency the matrix
# entries are passed in as arguments 'a11', 'a12', 'a21', 'a22'.
BindGlobal("Forms_TRANSFORM_2_BY_2",
  function(D,p1,p2,a11,a12,a21,a22)
    local r1,r2;
    r1 := D[p1];
    r2 := D[p2];
    D[p1] := a11 * r1 + a12 * r2;
    D[p2] := a21 * r1 + a22 * r2;
  end );

InstallGlobalFunction(Forms_REDUCE2,
  function(D,start,stop,gf)
    local n,t,i,half,primroot;
    n := NrRows(D);
    primroot := PrimitiveRoot(gf);
    half := One(gf) / 2;
    t := primroot^(LogFFE(-One(gf),primroot)/2) / 2;
    i := start;
    while i < stop do
      Forms_TRANSFORM_2_BY_2(D,i,i+1,half,t,half,-t);
      i := i + 2;
    od;
  end );

InstallGlobalFunction(Forms_REDUCE4,
  function(D,start,stop,gf)
    local n,c,d,i,dummy;
    n := NrRows(D);
    i := start;
    dummy := Forms_SUM_OF_SQUARES(-One(gf),gf);
    c := dummy[1];
    d := dummy[2];
    while i < stop do
      Forms_TRANSFORM_2_BY_2(D,i+1,i+3,c,d,d,-c);
      i := i + 4;
    od;
  end );

InstallGlobalFunction(Forms_DIFF_2_S,
  function(D,start,stop)
    local n,i,half;
    n := NrRows(D);
    i := start;
    half := One(D[1,1]) / 2;
    while i < stop do
      Forms_TRANSFORM_2_BY_2(D,i,i+1,half,half,half,-half);
      i := i + 2;
    od;
  end );

InstallGlobalFunction(Forms_HERM_CONJ,
  function(mat,t)
    local n,i,j,dummy;
    n := NrRows(mat);
    dummy := MutableTransposedMat(mat);
    for i in  [1..n] do
      for j in [1..n] do
        dummy[i,j] := dummy[i,j]^t;
      od;
    od;
    return dummy;
  end );

#Forms_RESET: given a matrix mat, computes an upper triangular matrix A such that mat and A determine
#the same quadratic form. The convention about quadratic forms in this package
#is that their Gram matrix is always an upper triangular matrix, although the user
#is free to use any matrix to construct the form.

BindGlobal("Forms_RESET_inplace",
  function(A)
    local i,j,t;
    t := Zero(A[1,1]);
    for i in [2..NrRows(A)] do
      for j in [1..i-1] do
        A[j,i] := A[j,i] + A[i,j];
        A[i,j] := t;
      od;
    od;
  end );

# HACK: ignore extra arguments to Forms_RESET, which used to
# take three arguments; the other two were redundant, so we dropped them.
# But the FinInG function also calls Forms_RESET, so for its sake, we
# ignore the other arguments
InstallGlobalFunction(Forms_RESET,
  function(A,extra...)
    A := MutableCopyMat(A);
    Forms_RESET_inplace(A);
    return A;
  end );

InstallGlobalFunction(Forms_SQRT2,
  function(a, gf)
    local z, q;
    if IsZero( a ) then
      return a;
    fi;
    q:=Size(gf);
    if q mod 2 = 0 then
      return a^(q/2);
    fi;
    z:=PrimitiveRoot(gf);
    return z^(LogFFE(a,z)/2);
  end );

InstallGlobalFunction(Forms_PERM_VAR,
  function(D,r)
    local i;
    i := Remove(D, r);
    Add(D, i, 1);
  end );

InstallGlobalFunction(Forms_C1,
  function(gf, h)
    local i, primroot;
    if h mod 2 = 0 then
      primroot := PrimitiveRoot(gf);
      i := 1;
      while i <= h - 1 do
        if not IsZero( Trace(gf, primroot^i) ) then
          return primroot^i;
        else
          i := i + 1;
        fi;
      od;
    else
      return One(gf);
    fi;
  end );

InstallGlobalFunction(Forms_QUAD_EQ,
  function(delta, gf, h)
    local i,k,dummy,result;
    k := Forms_C1(gf,h);
    dummy := Zero(gf);
    result := Zero(gf);
    for i in [1..h-1] do
      dummy := dummy + k^(2^(i-1));
      result := result + dummy*(delta^(2^i));
    od;
    return result;
  end );
##
#end support functions base change
#############################################################################

#############################################################################
# Operations to check input (most likely not for the user):
#############################################################################

InstallMethod( FORMS_IsSymplecticMatrix, [IsFFECollColl, IsField],
  function(m,f)
    if Characteristic(f) = 2 then
       if not ForAll([1..NrRows(m)], i -> IsZero(m[i,i]) ) then
           return false;
       fi;
    fi;
    return m = -TransposedMat(m);
  end );

InstallMethod( FORMS_IsSymmetricMatrix, [IsFFECollColl],
  function(m)
    return m=TransposedMat(m);
  end );

InstallMethod( FORMS_IsHermitianMatrix, [IsFFECollColl, IsField],
  function(m,f)
    return m=Forms_HERM_CONJ(m,Sqrt(Size(f)));
  end );

#############################################################################
# Main Base-change operations
# (helping methods for BaseChangeToCanonical, not for the users):
#############################################################################
#O BaseChangeOrthogonalBilinear( <mat>, <f> )
#  <f>: finite field, <mat>: matrix over <f>
# output: [D,r,w]; D = base change matrix,
#                r = number of non zero rows in D*mat*TransposedMat(D)
#                using r, it follows immediately whether the form is degenerate
#                w = character (0=elliptic, 2=hyperbolic, 1=parabolic).
##
InstallMethod( BaseChangeOrthogonalBilinear,
    [ IsMatrix and IsFFECollColl, IsField and IsFinite ],
  function(mat, gf)
    local row,i,j,k,A,b,c,d,P,D,dummy,r,w,s,v,v1,v2,
          n,q,primroot,one;
    Assert(1, mat = TransposedMat(mat));
    n := NrRows(mat);
    q := Size(gf);
    one := One(gf);

    A := MutableCopyMat(mat);
    ConvertToMatrixRep(A, gf);
    D := IdentityMat(n, gf);
    ConvertToMatrixRep(D, gf);
    row := 0;

    # Diagonalize A

    while row < n - 1 do
      row := row + 1;

      # We look for a nonzero element on the main diagonal, starting
      # from row
      i := row;
      while i <= n and IsZero(A[i,i]) do
        i := i + 1;
      od;

      if i = row then
        # do nothing since A[row,row] <> 0
      elif i <= n then
        # swap things around to ensure A[row,row] <> 0
        Forms_SwapCols(A, row, i);
        Forms_SwapRows(A, row, i);
        Forms_SwapRows(D, row, i);
      else
        # All entries on the main diagonal are zero. We now search for a
        # nonzero element off the main diagonal.
        i := row;
        while i < n do
          k := i + 1;
          while k <= n and IsZero(A[i,k]) do
            k := k + 1;
          od;
          if k = n + 1 then
             i := i + 1;
          else
             break;
          fi;
        od;

        # if i is n, then they are all zero and we can stop.
        if i = n then
          row := row - 1;
          r := row;
          break;
        fi;

        # Otherwise: Go and fetch...
        # Put it on A[row,row+1]
        if i <> row then
          Forms_SwapCols(A, row, i);
          Forms_SwapRows(A, row, i);
          Forms_SwapRows(D, row, i);
        fi;

        Forms_SwapCols(A, row + 1, k);
        Forms_SwapRows(A, row + 1, k);
        Forms_SwapRows(D, row + 1, k);

        b := 1/A[row+1,row];
        Forms_AddCols(A, row, row+1, b);
        Forms_AddRows(A, row, row+1, b);
        Forms_AddRows(D, row, row+1, b);
      fi;   # end if i = row ... elif  i <= n ... else ... fi

      # There is no zero element on the main diagonal, make the rest zero.
      c := -A[row,row]^-1;
      for i in [row+1..n] do
        b := A[i,row] * c;
        if IsZero(b) then continue; fi;
        Forms_AddCols(A, i, row, b);
        Forms_AddRows(A, i, row, b);
        Forms_AddRows(D, i, row, b);
      od;
    od;

    # Count how many variables are used.

    if row = n - 1 then
      if not IsZero(A[n,n]) then
        r := n;
      else
        r := n - 1;
      fi;
    fi;

    # Here we distinguish the quadratic from the non-quadratic

    i := 1;
    s := 0;
    primroot := PrimitiveRoot(gf);
    while i < r do
      if IsOddInt( LogFFE(A[i,i], primroot) ) then
         j := i + 1;
         repeat
            if IsEvenInt( LogFFE(A[j,j], primroot) ) then
               dummy := A[j,j];
               A[j,j] := A[i,i];
               A[i,i] := dummy;
               Forms_SwapRows(D, i, j);
               i := i + 1;
               s := s + 1;
               break;
            else
               j := j + 1;
            fi;
         until j = r + 1;
         if j = r + 1 then
            break;
         fi;
      else
        i := i + 1;
        s := s + 1;
      fi;
    od;
    if IsEvenInt( LogFFE(A[r,r], primroot) ) then
       s := s + 1;
    fi;

    # We do the form x_0^2 + ... + x_s^2 + v(x_s+1^2 + ... + x_r^2)
    # with v not quadratic.

    v := ShallowCopy(primroot);
    for i in [1..s] do
      Forms_MultRow(D,i,(primroot^(LogFFE(A[i,i], primroot)/2))^-1);
    od;
    for i in [s+1..r] do
      Forms_MultRow(D,i,(primroot^(LogFFE(A[i,i]/primroot,primroot)/2))^-1);
    od;

    # We keep as much quadratic part as we can:

    s := s - 1;
    r := r - 1;

    # Case by case:

    if not (s = -1 or r = s )  then
      # We write first v=v1^2 + v2^2
      dummy := Forms_SUM_OF_SQUARES(v,gf);
      v1 := dummy[1];
      v2 := dummy[2];

      if (r - s) mod 2 = 0 then
        v1 := v1/v;
        v2 := v2/v;
        i := s + 2;
        repeat
          Forms_TRANSFORM_2_BY_2(D,i,i+1,v1,-v2,v2,v1);
          i := i + 2;
        until i = r + 2;
        s := r;
      else
        if r mod 2 = 0 then
          i := 1;
          repeat
            Forms_TRANSFORM_2_BY_2(D,i,i+1,v1,v2,-v2,v1);
            i := i + 2;
          until i = s + 2;
          s := -1;
        elif not (s = r - 1) then
          v1 := v1/v;
          v2 := v2/v;
          i := s + 2;
          repeat
            Forms_TRANSFORM_2_BY_2(D,i,i+1,v1,-v2,v2,v1);
            i := i + 2;
          until i = r + 1;
          s := r - 1;
        fi;
      fi;
    fi;

    # Towards standard forms (uses the standard proof of
    # the classification of quadratic forms):

    if r mod 2 <> 0 then
      if s = -1 or s = r then
        if q mod 4 = 1 then
          Forms_REDUCE2(D,1,r+1,gf);
          w := 2;
        else
          if ((r-1)/2) mod 2 <> 0 then
            Forms_REDUCE4(D,1,r+1,gf);
            Forms_DIFF_2_S(D,1,r+1);
            w := 2;
          else
            Forms_REDUCE4(D,3,r+1,gf);
            Forms_DIFF_2_S(D,3,r+1);
            w := 0;
          fi;
        fi;
      else
        if q mod 4 = 1 then
          if 1 < r then
            Forms_SwapRows(D,2,r+1);
            Forms_REDUCE2(D,3,r+1,gf);
          fi;
          w := 0;
        else
          if ((r-1)/2) mod 2 <> 0 then
            Forms_SwapRows(D,4,r+1);
            if 3 < r then
              Forms_REDUCE4(D,5,r+1,gf);
            fi;
            b := primroot^(LogFFE(-v,primroot)/2);
            b := 1/(2*b);
            Forms_TRANSFORM_2_BY_2(D,3,4,one/2,-b,one/2,b);
            if 3 < r then
              Forms_DIFF_2_S(D,5,r+1);
            fi;
            w := 0;
          else
            Forms_SwapRows(D,2,r+1);
            if 1 < r then
              Forms_REDUCE4(D,3,r+1,gf);
            fi;
            b := primroot^(LogFFE(-v,primroot)/2);
            b := 1/(2*b);
            Forms_TRANSFORM_2_BY_2(D,1,2,one/2,-b,one/2,b);
            if 1 < r then
              Forms_DIFF_2_S(D,3,r+1);
            fi;
            w := 2;
          fi;
        fi;
      fi;
    elif r <> 0 then
      w := 1;
      if q mod 4 = 1 then
        Forms_REDUCE2(D,2,r+1,gf);
      else
        if r mod 4 = 0 then
          Forms_REDUCE4(D,2,r+1,gf);
          Forms_DIFF_2_S(D,2,r+1);
        else
          if 3 < r then
            Forms_REDUCE4(D,4,r+1,gf);
          fi;
          dummy := Forms_SUM_OF_SQUARES(-1,gf);
          c := dummy[1];
          d := dummy[2];
          Forms_TRANSFORM_2_BY_2(D,1,3,c,d,d,-c);
          Forms_DIFF_2_S(D,2,r+1);
          i := 3;
          while i <= r + 1 do
            Forms_MultRow(D,i,-one);
            i := i + 2;
          od;
        fi;
      fi;
    else
      w := 1;
    fi;

    return [D,r,w];
end);

#############################################################################
#O BaseChangeOrthogonalQuadratic( <mat>, <f> )
#  <f>: finite field, <mat>: matrix over <f>
# output: [D,r,w]; D = base change matrix,
#                r = number of non zero rows in D*mat*TransposedMat(D)
#                using r, it follows immediately whether the form is degenerate
#                w = character (0=elliptic, 2=hyperbolic, 1=parabolic).
##
InstallMethod(BaseChangeOrthogonalQuadratic, [ IsMatrix and IsFFECollColl, IsField and IsFinite ],
    function(mat, gf)
    local A,r,w,row,dummy,i,j,h,D,P,t,a,b,c,d,e,s,
      zeros,posr,posk,n,zero,one;
    n := NrRows(mat);
    r := n;
    row := 1;
    zero := Zero(gf);
    one := One(gf);
    A := MutableCopyMat(mat);
    D := IdentityMat(n, gf);
    ConvertToMatrixRep(D, gf);
    zeros := [];
    for i in [1..n] do
      zeros[i] := zero;
    od;
    h := DegreeOverPrimeField(gf);
    while row + 2 <= r do
      if not IsZero( A[row,row] ) then
        i := row + 1;
        # check on the main diagonal; we look for a zero
        while i <= r and not IsZero(A[i,i]) do
          i := i + 1;
        od;

        # if there is a zero somewhere, then we go and get it.

        if i <= r then
          Forms_SwapCols(A,row,i);
          Forms_SwapRows(A,row,i);
          Forms_SwapRows(D,row,i);
          Forms_RESET_inplace(A);

        # Otherwise: look in other places.

        else
          dummy := true;
          i := row;
          while i <= r - 1 and dummy do
            j := i + 1;
            while j <= r do
              if not IsZero( A[i,j] ) then
                posr := i;
                posk := j;
                dummy := false;
                break;
              else
                j := j + 1;
              fi;
            od;
            i := i + 1;
          od;

          # If all is zero, STOP
          if dummy then
            t := Forms_SQRT2(A[row,row],gf);
            Forms_MultRow(D,row,1/t);
            for i in [row + 1..r] do
              Forms_AddRows(D,i,row,Forms_SQRT2(A[i,i],gf));
            od;
            # Permutation of the variables, it is a parabolic
            r := row;
            Forms_PERM_VAR(D,r);
            w := 1;
            r := r - 1;
            return [D,r,w];
          # Otherwise: A basischange
          else
            if IsZero( A[row+1,row+2] ) then
               if posr = row + 1 then
                  Forms_SwapCols(A,posk,row+2);
                  Forms_SwapRows(A,posk,row+2);
                  Forms_SwapRows(D,posk,row+2);
               elif posk = row + 2 then
                  Forms_SwapCols(A,posr,row+1);
                  Forms_SwapRows(A,posr,row+1);
                  Forms_SwapRows(D,posr,row+1);
               elif posr = row + 2 then
                  # TODO: Does this case ever occur? I failed to find examples
                  # that trigger it
                  P := TransposedMat(PermutationMat((posk,posr,row+1),n));
                  A := P*A*TransposedMat(P);
                  D := P*D;
               else
                  Forms_SwapCols(A,posk,row+2);
                  Forms_SwapRows(A,posk,row+2);
                  Forms_SwapRows(D,posk,row+2);

                  Forms_SwapCols(A,posr,row+1);
                  Forms_SwapRows(A,posr,row+1);
                  Forms_SwapRows(D,posr,row+1);
               fi;
               Forms_RESET_inplace(A);
            fi;
            #A[row+1,row+2] <> 0
            t := A[row+1,row+2];

            b := 1/t;
            Forms_MultCol(A,row+2,b);
            Forms_MultRow(A,row+2,b);
            Forms_MultRow(D,row+2,b);

            b := A[row,row+1];
            Forms_AddCols(A,row,row+2,b);
            Forms_AddRows(A,row,row+2,b);
            Forms_AddRows(D,row,row+2,b);

            for i in [row+3..n] do
              b := A[row+1,i];
              Forms_AddCols(A,i,row+2,b);
              Forms_AddRows(A,i,row+2,b);
              Forms_AddRows(D,i,row+2,b);
            od;

            Forms_RESET_inplace(A);

            # A has now that special form a_11*X_1^2+X_1*X_2 + G(X_0,X_2,...,X_n);
            b := A[row,row];
            t :=  Forms_SQRT2(b/A[row+1,row+1],gf);
            Forms_AddCols(A,row,row+1,t);
            Forms_AddRows(A,row,row+1,t);
            Forms_AddRows(D,row,row+1,t);

            #A [row,row] is now 0
            Forms_RESET_inplace(A);
          fi;
        fi;
      fi;
      # check for zero row
      dummy := true;
      i := row + 1;
      while i <= n do
        if not IsZero( A[row,i] ) then
           dummy := false;
           break;
        else
           i := i + 1;
        fi;
      od;
      posk := i;
      # A is a zero row then...
      if dummy then
        for i in [row+1..r] do
          Forms_SwapCols(A,i,i-1);
          Forms_SwapRows(A,i,i-1);
          Forms_SwapRows(D,i,i-1);
        od;
        r := r - 1;
      else
        if posk <> row + 1 then
          Forms_SwapCols(A,posk,row+1);
          Forms_SwapRows(A,posk,row+1);
          Forms_SwapRows(D,posk,row+1);
          Forms_RESET_inplace(A);
        fi;

        # Now A[k,k+1] <> 0
        t := A[row,row+1];
        b := 1/t;
        Forms_MultCol(A,row+1,b);
        Forms_MultRow(A,row+1,b);
        Forms_MultRow(D,row+1,b);
        for i in [row+2..n] do
          b := A[row,i];
          Forms_AddCols(A,i,row+1,b);
          Forms_AddRows(A,i,row+1,b);
          Forms_AddRows(D,i,row+1,b);
        od;
        Forms_RESET_inplace(A);

        for i in [row+1..n] do
          b := A[row+1,i];
          Forms_AddCols(A,i,row,b);
          Forms_AddRows(A,i,row,b);
          Forms_AddRows(D,i,row,b);
        od;
        Forms_RESET_inplace(A);
        row := row + 2;
      fi;
    od;
    # Now there can be at most two variables left.
    # Case by case:

    if r = row then
       if IsZero(A[row,row]) then
          r := r - 1;
          w := 2;
       else
          t := Forms_SQRT2(A[r,r],gf);
          Forms_MultRow(D,r,1/t);
          Forms_PERM_VAR(D,r);
          w := 1;
       fi;
    else
       a := A[row,row];
       b := A[row,row+1];
       c := A[row+1,row+1];
       t := zero;
       if a = t then
          if b = t then
             if c = t then
             r := r - 2;
             w := 2;
          else
             Forms_MultRow(D,r,1/Forms_SQRT2(c,gf));
             Forms_PERM_VAR(D,r);
             r := r - 1;
             w := 1;
          fi;
        else
          if c = t then
            Forms_MultRow(D,r,1/b);
          else
            Forms_MultRow(D,r-1,1/b);
            Forms_AddRows(D,r,r-1,c);
          fi;
          w := 2;
        fi;
      else #a <> t
        if b = t then
          if c = t then
            Forms_MultRow(D,r-1,1/Forms_SQRT2(a,gf));
            Forms_PERM_VAR(D,r-1);
            r := r - 1;
          else
            Forms_MultRow(D,r-1,1/Forms_SQRT2(a,gf));
            Forms_AddRows(D,r,r-1,Forms_SQRT2(c,gf));
            Forms_PERM_VAR(D,r-1);
            r := r - 1;
          fi;
          w := 1;
        else
          if c = t then
            Forms_MultRow(D,r,1/b);
            Forms_AddRows(D,r-1,r,a);
            w := 2;
          else
            d := (a*c)/(b^2);
            if Trace(gf,d) = t then
              e := Forms_SQRT2(a,gf);
              s := Forms_QUAD_EQ(d,gf,h);
              Forms_TRANSFORM_2_BY_2(D, r-1, r, (s+one)/e, e/b, s/e, e/b);
              w := 2;
            else
              c := Forms_SQRT2(c,gf);
              Forms_MultRow(D,r-1,c/b);
              Forms_MultRow(D,r,1/c);
              if r > 2 then
                Forms_SwapRows(D,1,r);
                Forms_SwapRows(D,2,r-1);
              else
                Forms_SwapRows(D,1,2);
              fi;
              e := Forms_C1(gf,h);
              if e <> d then
                 a := Forms_QUAD_EQ(d+e,gf,h);
                 Forms_AddRows(D,2,1,a);
              fi;
              w := 0;
            fi;
          fi;
        fi;
      fi;
    fi;
    r := r - 1;
    return [D,r,w];
end );

#############################################################################
#O  BaseChangeHermitian( <mat>, <field> )
# input: Gram matrix of a hermitian form, field
# output: [D,r]; D = base change matrix,
#                r = number of non zero rows in D*mat*TransposedMat(D)
#                using r, it follows immediately whether the form is degenerate
##
InstallMethod(BaseChangeHermitian, [ IsMatrix and IsFFECollColl, IsField and IsFinite ],
  function(mat,gf)
    local row,i,j,k,A,a,b,P,D,t,r,n,one,A2,D2;
    n := NrRows(mat);
    one := One(gf);
    t := Sqrt(Size(gf));

    A := MutableCopyMat(mat);
    ConvertToMatrixRep(A, gf);
    D := IdentityMat(n, gf);
    ConvertToMatrixRep(D, gf);
    row := 0;

    # Diagonalize A

    while row < n - 1 do
      row := row + 1;

      # We look for a nonzero element on the main diagonal, starting
      # from row
      i := row;
      while i <= n and IsZero(A[i,i]) do
        i := i + 1;
      od;

      if i = row then
        # do nothing since A[row,row] <> 0
      elif i <= n then
        # swap things around to ensure A[row,row] <> 0
        Forms_SwapCols(A, row, i);
        Forms_SwapRows(A, row, i);
        Forms_SwapRows(D, row, i);
      else
        # All entries on the main diagonal are zero. We now search for a
        # nonzero element off the main diagonal.
        i := row;
        while i < n do
          k := i + 1;
          while k <= n and IsZero(A[i,k]) do
            k := k + 1;
          od;
          if k = n + 1 then
             i := i + 1;
          else
             break;
          fi;
        od;

        # if i is n, then they are all zero and we can stop.
        if i = n then
          row := row - 1;
          r := row;
          break;
        fi;

        # Otherwise: Go and fetch...
        # Put it on A[row,row+1]
        if i <> row then
          Forms_SwapCols(A, row, i);
          Forms_SwapRows(A, row, i);
          Forms_SwapRows(D, row, i);
        fi;

        Forms_SwapCols(A, row + 1, k);
        Forms_SwapRows(A, row + 1, k);
        Forms_SwapRows(D, row + 1, k);

        b := PrimitiveRoot(gf)/A[row+1,row];
        Forms_AddCols(A, row, row+1, b^t);
        Forms_AddRows(A, row, row+1, b);
        Forms_AddRows(D, row, row+1, b);
      fi;

      # There is no zero element on the main diagonal, make the rest zero.
      a := -A[row,row]^-1;
      for i in [row+1..n] do
        b := A[i,row] * a;
        if IsZero(b) then continue; fi;
        Forms_AddCols(A, i, row, b^t);
        Forms_AddRows(A, i, row, b);
        Forms_AddRows(D, i, row, b);
      od;
    od;

    # Count how many variables have been used

    if row = n - 1 then
      if not IsZero(A[n,n]) then
        r := n;
      else
        r := n - 1;
      fi;
    fi;

    # Take care that the diagonal elements become 1.

    for i in [1..r] do
      a := A[i,i];
      if not IsOne(a) then
        # find an element b with norm b*b^t = b^(t+1) equal to a
        b := RootFFE(gf, a, t+1);
        Forms_MultRow(D,i,1/b);
      fi;
    od;
    return [D,r-1];
end );

#############################################################################
##
#O  BaseChangeSymplectic( <mat>, <field> )
# input: Gram matrix of a symplectic form, field
# output: [D,r]; D = base change matrix,
#                r = number of non zero rows in D*mat*TransposedMat(D)
#                using r, the Witt index of the non-degenerate part can be computed.
##
InstallMethod( BaseChangeSymplectic, [IsMatrix and IsFFECollColl, IsField and IsFinite],

## This operation returns an isometry g such that g m g^T is
## the alternating form arising from the block diagonal matrix
## with each block equal to J=[[0,1],[-1,0]].

 function(m, f)
   local d, basechange, blocknr, diagpos, pos, j, a, b, offset;
   d := NrRows(m);
   basechange := IdentityMat(d, f);
   ConvertToMatrixRep(basechange, f);
   m := MutableCopyMat(m);
   for blocknr in [1 .. (Int(d/2))] do
      ## diagpos is the position of the top left corner of the block
      ## on the diagonal of m
      diagpos := 2 * blocknr - 1;

      for offset in [0..d-diagpos] do
         ## find first nonzero entry in column diagpos below the diagonal
         pos := First([diagpos+1..d], j -> not IsZero( m[diagpos+offset,j] ) );
         if pos <> fail then
            if offset <> 0 then
               Forms_SwapRows(basechange, diagpos, diagpos+offset);
               Forms_SwapRows(m, diagpos, diagpos+offset);
               Forms_SwapCols(m, diagpos, diagpos+offset);
               pos := First([diagpos+1..d], j -> not IsZero( m[diagpos,j] ) );
               Assert(0, pos <> fail);
            fi;
            break;
         fi;
      od;

      # when pos=fail, then only degeneracy is left and the base transition is complete
      if pos=fail then
         return [basechange, 2*(blocknr-1)];
                #the second entry is the number of non-zero rows after basechange
      fi;

      # normalize the non-zero entry we just located to be 1
      a := 1/m[diagpos,pos];
      Forms_MultRow(basechange, pos, a);
      Forms_MultRow(m, pos, a);
      Forms_MultCol(m, pos, a);

      #
      if pos <> diagpos+1 then
         Forms_SwapRows(basechange, diagpos+1, pos);
         Forms_SwapRows(m, diagpos+1, pos);
         Forms_SwapCols(m, diagpos+1, pos);
      fi;

      # clean up to the right and below the current block
      for j in [2 * blocknr + 1..d] do
         a := -m[j,diagpos+1];
         b := m[j,diagpos];

         Forms_AddRows(basechange, j, diagpos, a);
         Forms_AddRows(basechange, j, diagpos+1, b);

         Forms_AddRows(m, j, diagpos, a);
         Forms_AddRows(m, j, diagpos+1, b);

         Forms_AddCols(m, j, diagpos, a);
         Forms_AddCols(m, j, diagpos+1, b);
      od;
   od;
   return [basechange,2*blocknr];
          #the second entry is the number of non-zero rows after basechange
 end );

#############################################################################
# Other Operations:
#############################################################################

InstallMethod( BaseChangeHomomorphism, [ IsMatrix and IsFFECollColl, IsField ],
  function( b, gf )
  ## This function returns an intertwiner of the general linear group
  ## induced by changing the basis of its underlying vector space

    local gl, invb, hom;
    if IsZero(Determinant(b)) then
       Error("Matrix is not invertible");
    fi;
    invb := Inverse( b );
    gl := GeneralLinearGroup(Size(b), gf);
    hom := InnerAutomorphismNC( gl, invb);
    return hom;
  end );

#############################################################################
# Attributes depending on base change possibility.
#############################################################################
#A  WittIndex( <form> )
# input: bilinear form
# output <r>. r = number of non zero rows in D*mat*TransposedMat(D) =: witt index
#             (in fact only when form is non-degenerate).
# note: calling this attribute will also set properties like
#         IsElliptic,..., and BaseChangeCanonical
##
InstallMethod( WittIndex, "for a bilinear form",
  [IsBilinearForm and IsFormRep],
  function(f)
    local string, b;
    string := f!.type;
    if IsSymplecticForm(f) then
       b := BaseChangeSymplectic(f!.matrix, f!.basefield);
       SetBaseChangeToCanonical(f, b[1]);
       return (b[2]/2);
    elif IsOrthogonalForm(f) then
       b := BaseChangeOrthogonalBilinear(f!.matrix, f!.basefield);
       SetIsEllipticForm(f,b[3]=0);
       SetIsParabolicForm(f,b[3]=1);
       SetIsHyperbolicForm(f,b[3]=2);
       SetBaseChangeToCanonical(f,b[1]);
       return (b[2]+b[3]-1)/2;
    elif string = "pseudo" then
       Error("Witt Index not yet implemented for pseudo forms");
    else
       Error("Form is incorrectly specified");
    fi;
  end );

#############################################################################
#A  WittIndex( <form> )
# input: quadratic form
# output <r>. r = number of non zero rows in D*mat*TransposedMat(D) =: witt index
#             (in fact only when form is non-degenerate).
# note: calling this attribute will also set properties like IsElliptic,...
##
InstallMethod( WittIndex, "for a quadratic form",
  [IsQuadraticForm and IsFormRep],
  function(f)
    local m, gf, b, polarity;
    m := f!.matrix;
    gf := f!.basefield;
    if IsOddInt(Size(gf)) then
      polarity := (m+TransposedMat(m))/2;
      b := BaseChangeOrthogonalBilinear(polarity, gf);
    else
      b := BaseChangeOrthogonalQuadratic(m,gf);
    fi;
    SetBaseChangeToCanonical(f,b[1]);
    SetIsEllipticForm(f,b[3]=0);
    SetIsParabolicForm(f,b[3]=1);
    SetIsHyperbolicForm(f,b[3]=2);
    return (b[2]+b[3]-1)/2;
  end );

#############################################################################
#A  WittIndex( <form> )
# input: hermitian form
# output <r>. r = number of non zero rows in D*mat*TransposedMat(D) =: witt index
#             (in fact only when form is non-degenerate).
# note: calling this attribute will also set BaseChangeToCanonical
##
InstallMethod( WittIndex, "for a hermitian form",
  [IsHermitianForm and IsFormRep],
  function(f)
    local gram, gf, b;
    gram := f!.matrix;
    gf := f!.basefield;
    b := BaseChangeHermitian(gram, gf);
      #returns [D,r] with r = non-zero rows - 1 after base change
    SetBaseChangeToCanonical(f, b[1]);
    return Int( (b[2]+1)/2 );
  end );

#############################################################################
#A  WittIndex( <form> )
# input: trivial form
# output error message
##
InstallMethod( WittIndex, "for a trivial form",
  [IsTrivialForm],
  function(f)
    Error("<form> is trivial, Witt Index not defined");
end );

#############################################################################
#  Properties again (for users).
#  see package documentation for more information.
InstallMethod( IsEllipticForm, "for sesquilinear forms",
  [IsSesquilinearForm and IsFormRep],
  function(f)
    local b,gf,m;
    gf := f!.basefield;
    m := f!.matrix;
    if IsOrthogonalForm(f) then
       b := BaseChangeOrthogonalBilinear(m, gf);
       SetWittIndex(f, (b[2]+b[3]-1) / 2);
       SetBaseChangeToCanonical(f, b[1]);
       SetIsParabolicForm(f,b[3]=1);
       SetIsHyperbolicForm(f,b[3]=2);
       return b[3] = 0;
    else
       return false;
    fi;
  end );

InstallMethod( IsEllipticForm,  "for quadratic forms",
  [IsQuadraticForm and IsFormRep],
  function(f)
    local b,gf,m,polarity;
    gf := f!.basefield;
    m := f!.matrix;
    if IsOddInt(Size(gf)) then
      polarity := (m+TransposedMat(m))/2;
      b := BaseChangeOrthogonalBilinear(polarity, gf);
    else
      b := BaseChangeOrthogonalQuadratic(m,gf);
    fi;
    SetWittIndex(f, (b[2]+b[3]-1) / 2);
    SetBaseChangeToCanonical(f, b[1]);
    SetIsParabolicForm(f,b[3]=1);
    SetIsHyperbolicForm(f,b[3]=2);
    return b[3] = 0;
 end );

InstallMethod( IsEllipticForm, [IsTrivialForm], #new in 1.2.1
 function( f )
   return false;
 end );


InstallMethod( IsParabolicForm, "for sesquilinear forms",
  [IsSesquilinearForm and IsFormRep],
  function(f)
  local b,gf,m;
    gf := f!.basefield;
    m := f!.matrix;
    if IsOrthogonalForm(f) then
      b := BaseChangeOrthogonalBilinear(m, gf);
      SetWittIndex(f, (b[2]+b[3]-1) / 2);
      SetBaseChangeToCanonical(f, b[1]);
      SetIsEllipticForm(f,b[3]=0);
      SetIsHyperbolicForm(f,b[3]=2);
      return b[3] = 1;
    else
      return false;
    fi;
  end );

InstallMethod( IsParabolicForm, "for quadratic forms",
  [IsQuadraticForm and IsFormRep],
  function(f)
  local b,gf,m,polarity;
    gf := f!.basefield;
    m := f!.matrix;
    if IsOddInt(Size(gf)) then
      polarity := (m+TransposedMat(m))/2;
      b := BaseChangeOrthogonalBilinear(polarity, gf);
    else
      b := BaseChangeOrthogonalQuadratic(m,gf);
    fi;
    SetWittIndex(f, (b[2]+b[3]-1) / 2);
    SetBaseChangeToCanonical(f, b[1]);
    SetIsEllipticForm(f,b[3]=0);
    SetIsHyperbolicForm(f,b[3]=2);
    return b[3] = 1;
  end );

InstallMethod( IsParabolicForm, [IsTrivialForm],
  function( f )
    return false;
  end );


InstallMethod( IsHyperbolicForm, "for sesquilinear forms",
  [IsSesquilinearForm and IsFormRep],
  function(f)
    local b,gf,m;
    gf := f!.basefield;
    m := f!.matrix;
    if IsOrthogonalForm(f) then
      b := BaseChangeOrthogonalBilinear(m, gf);
      SetWittIndex(f, (b[2]+b[3]-1) / 2);
      SetBaseChangeToCanonical(f, b[1]);
      SetIsEllipticForm(f,b[3]=0);
      SetIsParabolicForm(f,b[3]=1);
      return b[3] = 2;
    else
      return false;
    fi;
  end );

InstallMethod( IsHyperbolicForm, "for quadratic forms",
  [IsQuadraticForm and IsFormRep],
  function(f)
    local b,gf,m,polarity;
    gf := f!.basefield;
    m := f!.matrix;
    if IsOddInt(Size(gf)) then
       polarity := (m+TransposedMat(m))/2;
       b := BaseChangeOrthogonalBilinear(polarity, gf);
    else
       b := BaseChangeOrthogonalQuadratic(m,gf);
    fi;
    SetWittIndex(f, (b[2]+b[3]-1) / 2);
    SetBaseChangeToCanonical(f, b[1]);
    SetIsEllipticForm(f,b[3]=0);
    SetIsParabolicForm(f,b[3]=1);
    return b[3] = 2;
 end );

InstallMethod( IsHyperbolicForm, [IsTrivialForm],
  function( f )
    return false;
  end );
##
#############################################################################

#############################################################################
#  Attribute again (for users).
#############################################################################

#############################################################################
#A  IsometricCanonicalForm( <form> )
# input: trivial form
# output: trivial form and a warning.
###
InstallMethod( IsometricCanonicalForm, "for bilinear forms",
  [IsTrivialForm],
  function(f)
    Info(InfoWarning,1,"<form> is trivial");
    return f;
  end );

#############################################################################
#A  IsometricCanonicalForm( <form> )
# input: bilinear form
# output: isometric canonical form of <form>. Will also set attributes and properties of output.
###
InstallMethod( IsometricCanonicalForm, "for bilinear forms",
  [IsBilinearForm and IsFormRep],
  function(f)
    local B,gram,isom,gf,form,trivial;
    gram := GramMatrix(f);
    B := BaseChangeToCanonical(f);
    isom := B*gram*TransposedMat(B);
    gf := f!.basefield;
    form := BilinearFormByMatrix(isom,gf);
    trivial := IdentityMat(NrRows(gram),gf);
    SetBaseChangeToCanonical(form,trivial);
    SetWittIndex(form, WittIndex(f));
    if IsOrthogonalForm(f) then
       SetIsOrthogonalForm(form,true);
       SetIsEllipticForm(form,IsEllipticForm(f));
       SetIsParabolicForm(form,IsParabolicForm(f));
       SetIsHyperbolicForm(form,IsHyperbolicForm(f));
    fi;
    return form;
  end );

#############################################################################
#A  IsometricCanonicalForm( <form> )
# input: hermitian form
# output: isometric canonical form of <form>. Will also set attributes and properties of output.
##
InstallMethod( IsometricCanonicalForm, "for hermitian forms",
  [IsHermitianForm and IsFormRep],
  function(f)
    local gram,gf,B,isom,form,trivial;
    gram := f!.matrix;
    gf := f!.basefield;
    trivial := IdentityMat(NrRows(gram),gf);
    B := BaseChangeToCanonical(f);
    isom := B*gram*Forms_HERM_CONJ(B,Sqrt(Size(gf)));
    form := FormByMatrix(isom,gf,"hermitian");
    SetBaseChangeToCanonical(form,trivial);
    SetWittIndex(form, WittIndex(f));
    return form;
  end );

#############################################################################
#A  IsometricCanonicalForm( <form> )
# input: quadratic form
# output: isometric canonical form of <form>. Will also set attributes and properties of output.
##
InstallMethod( IsometricCanonicalForm, "for quadratic forms",
  [IsQuadraticForm and IsFormRep],
  function(f)
    local B,gram,isom,gf,form,trivial;
    gram := GramMatrix(f);
    gf := f!.basefield;
    trivial := IdentityMat(NrRows(gram),gf);
    B := BaseChangeToCanonical(f);
    isom := Forms_RESET(B*gram*TransposedMat(B));
    form := FormByMatrix(isom,gf,"quadratic");
    SetBaseChangeToCanonical(form,trivial);
    SetWittIndex(form, WittIndex(f));
    SetIsEllipticForm(form,IsEllipticForm(f));
    SetIsParabolicForm(form,IsParabolicForm(f));
    SetIsHyperbolicForm(form,IsHyperbolicForm(f));
    return form;
end );

#############################################################################
#A  IsometricCanonicalForm( <form> )
# input: trivial form
# output: trivial form
##
InstallMethod( IsometricCanonicalForm, "for trivial forms",
  [IsTrivialForm],
  function(f)
    return f;
end );

#############################################################################
# Evaluation Methods (for users):
#############################################################################

InstallMethod( EvaluateForm,  "for a bilinear form and a pair of vectors",
  [IsBilinearForm and IsFormRep,
        IsVector and IsFFECollection, IsVector and IsFFECollection],
  function(f,v,w)
    return v*f!.matrix*w;
end );

InstallMethod( EvaluateForm,  "for a bilinear form and a pair of matrices",
  [IsBilinearForm and IsFormRep, IsFFECollColl, IsFFECollColl],
  function(f,v,w)
    return v*f!.matrix*TransposedMat(w);
end );

InstallMethod( EvaluateForm, "for an hermitian form and a pair of vectors",
  [IsHermitianForm and IsFormRep,
        IsVector and IsFFECollection, IsVector and IsFFECollection],
  function(f,v,w)
    local gf,t,hh,frob;
    gf := f!.basefield;
    hh := DegreeOverPrimeField(gf) / 2;
    frob := FrobeniusAutomorphism(gf)^hh;
    return v*f!.matrix*(w^frob);
end );

InstallMethod( EvaluateForm,  "for an hermitian form and a pair of matrices",
  [IsHermitianForm and IsFormRep, IsFFECollColl, IsFFECollColl],
  function(f,v,w)
    local gf,t,wCONJ,m,n,i,j;
    gf := f!.basefield;
    t := Sqrt(Size(gf));
    wCONJ := MutableTransposedMat(w);
    m := NrRows(wCONJ);
    n := NrCols(wCONJ);
    for i in [1..m] do
      for j in [1..n] do
        wCONJ[i,j] := wCONJ[i,j]^t;
      od;
    od;
    return v*f!.matrix*wCONJ;
end );

InstallMethod( EvaluateForm, "for quadratic forms",
  [IsQuadraticForm and IsFormRep, IsVector and IsFFECollection],
  function(f,v)
    return v*f!.matrix*v;
end );

InstallMethod( EvaluateForm,
    "for a quadratic form and an FFE matrix",
    [ IsQuadraticForm, IsMatrix and IsFFECollColl ],
    function(f,m)
        return m * f!.matrix * TransposedMat(m);
    end );

InstallMethod( EvaluateForm,  "for trivial forms and a pair of vectors",
  [IsTrivialForm,
        IsVector and IsFFECollection, IsVector and IsFFECollection],
  function(f,v,w)
    return Zero(BaseField(f));
end );

InstallMethod( EvaluateForm,  "for trivial forms and a pair of matrices",
  [IsTrivialForm, IsFFECollColl, IsFFECollColl],
  function(f,v,w)
    return Zero(BaseField(f));
end );

InstallMethod( EvaluateForm, "for trivial forms",
  [IsTrivialForm, IsVector and IsFFECollection],
  function(f,v)
    return Zero(BaseField(f));
end );

#############################################################################
# Methods to compute orthogonal subspaces to a given subspace wrt
# sesquilinear forms (for users).
#############################################################################
#############################################################################
#O OrthogonalSubspaceMat( <form>, <v> ) <form>: bil. form.
#  <v>: vector. Returns base of subspace orthogonal to <v> wrt <form>.
##
InstallMethod(OrthogonalSubspaceMat,
  "for a form and a vector",
  [IsBilinearForm, IsVector and IsFFECollection],
  function(f,v)
  local mat;
  mat := f!.matrix;
  if Length(v) <> NrRows(mat) then
    Error("<v> has the wrong dimension");
  fi;
  return NullspaceMat(TransposedMat([v*mat]));
end );

#############################################################################
#O OrthogonalSubspaceMat( <form>, <sub> ) <form>: bil. form.
#  <sub>: base of subspace. Returns base of subspace orthogonal to <sub> wrt <form>.
##
InstallMethod(OrthogonalSubspaceMat,
  "for a form and a basis of a subspace",
  [IsBilinearForm, IsMatrix],
  function(f,sub)
  local mat,perp;
  mat := f!.matrix;
  if Length(sub[1]) <> NrRows(mat) then
    Error("<sub> contains vectors of wrong dimension");
  fi;
  perp := TransposedMat(sub*mat);
  return NullspaceMat(perp);
end );

#############################################################################
#O OrthogonalSubspaceMat( <form>, <v> ) <form>: herm. form.
#  <v>: vector. Returns base of subspace orthogonal to <v> wrt <form>.
##
InstallMethod(OrthogonalSubspaceMat,
  "for a form and a vector",
  [IsHermitianForm, IsVector and IsFFECollection],
  function(f,v)
  local mat,gf,t,vt;
  mat := f!.matrix;
  if Length(v) <> NrRows(mat) then
    Error("<v> has the wrong dimension");
  fi;
  gf := f!.basefield;
  t := Sqrt(Size(gf));
  vt := List(v,x->x^t);
  return NullspaceMat(mat*TransposedMat([vt]));
end );

#############################################################################
#O OrthogonalSubspaceMat( <form>, <sub> ) <form>: herm. form.
#  <sub>: base of subspace. Returns base of subspace orthogonal to <sub> wrt <form>.
##
InstallMethod(OrthogonalSubspaceMat,
  "for a form and a basis of a subspace",
  [IsHermitianForm, IsMatrix],
  function(f,sub)
  local mat,gf,t,subt;
  mat := f!.matrix;
  if Length(sub[1]) <> NrRows(mat) then
    Error("<sub> contains vectors of wrong dimension");
  fi;
  gf := f!.basefield;
  t := Sqrt(Size(gf));
  subt := List(sub,x->List(x,y->y^t));
  return NullspaceMat(mat*TransposedMat(subt));
end );

#############################################################################
# Methods to compute orthogonal subspaces to a given subspace wrt
# quadratic forms (for users). Remember the subtle difference between
# orthogonality and singularity for subspaces.
#############################################################################
#############################################################################
#O OrthogonalSubspaceMat( <form>, <v> ) <form>: quad. form.
#  <v>: vector. Returns base of subspace orthogonal to <v> wrt associated
# bilinear form of <form>.
##
InstallMethod(OrthogonalSubspaceMat,
  "for a form and a vector",
  [IsQuadraticForm, IsVector and IsFFECollection],
  function(f,v)
  local bilf;
  bilf := AssociatedBilinearForm(f);
  return OrthogonalSubspaceMat(bilf,v); #note that this call will perform dim
end );                  #checks

#############################################################################
#O OrthogonalSubspaceMat( <form>, <sub> ) <form>: quad. form.
#  <sub>: base of subspace. Returns base of subspace orthogonal to <sub> wrt
# associated bilinear form of <form>.
InstallMethod(OrthogonalSubspaceMat,
  "for a form and a basis of a subspace",
  [IsQuadraticForm, IsMatrix],
  function(f,sub)
  local bilf;
  bilf := AssociatedBilinearForm(f);
  return OrthogonalSubspaceMat(bilf,sub); #note that this call will perform dim
end );                  #checks

#############################################################################
# Methods for IsIsotropicVector and IsTotallyIsotropicSubspace for sesquilinear
# forms (for users)
#############################################################################
#############################################################################
#O IsIsotropicVector( <form>, <v> ) <form>: sesq. form.
#  returns true if and only if <v> is a totally isotropic vector wrt <form>.
##
InstallMethod(IsIsotropicVector,
  "for a form and a vector",
  [IsSesquilinearForm, IsVector and IsFFECollection],
  function(f,v)
  return IsZero( [v,v]^f );
end );

#############################################################################
#O IsTotallyIsotropicSubspace( <form>, <sub> ) <form>: sesq. form.
#  returns true if and only if <sub> is totally isotropic subspace wrt <form>.
##
InstallMethod(IsTotallyIsotropicSubspace,
  "for a form and a basis of a subspace",
  [IsSesquilinearForm, IsMatrix],
  function(f,sub)
    local mat;
    mat := f!.matrix;
    if f!.type = "hermitian" then
       #return IsZero( (sub^CompanionAutomorphism( f )) * mat * TransposedMat(sub) );
       #the next line replaces the previous one, repairing a bug found by John Bamberg.
       return IsZero( sub * mat * TransposedMat(sub^CompanionAutomorphism( f )) );
    else
       return IsZero( sub * mat * TransposedMat(sub) );
    fi;
end );

#############################################################################
# Methods for IsIsotropicVector and IsTotallyIsotropicSubspace for quadratic
# forms (look at the definitions!);
#############################################################################
#############################################################################
#O IsIsotropicVector( <form>, <v> ) <form>: quad. form.
# returns true if and only if <v> is a totally isotropic vector wrt to
# associated bilinear form of<form>.
##
InstallMethod(IsIsotropicVector,
  "for a quadratic form and a vector",
  [IsQuadraticForm, IsVector and IsFFECollection],
  function(f,v)
  return IsIsotropicVector(AssociatedBilinearForm(f),v); #performs dim checks
end );

#############################################################################
#O IsTotallyIsotropicSubspace( <form>, <sub> ) <form>: quad. form.
#  returns true if and only if <sub> is totally isotropic subspace wrt to
# associated bilinear form of <form>.
##
InstallMethod(IsTotallyIsotropicSubspace,
  "for a quadratic form and a basis of a subspace",
  [IsQuadraticForm, IsMatrix],
  function(f,sub)
  return IsTotallyIsotropicSubspace(AssociatedBilinearForm(f),sub); #performs dc.
end );

#############################################################################
# Methods for IsSingularVector and IsTotallySingularSubspace for quadratic
# forms (look at the definitions!);
#############################################################################
#############################################################################
#O IsSingularVector( <form>, <v> ) <form>: quad. form.
# returns true if and only if <v> is a singular vector wrt to <form>.
##
InstallMethod(IsSingularVector,
  "for a quadratic form and a vector",
  [IsQuadraticForm, IsVector and IsFFECollection],
  function(f,v)
  return v^f=Zero(BaseField(f));
end );

#############################################################################
#O IsTotallySingularSubspace( <form>, <sub> ) <form>: quad. form.
#  returns true if and only if <sub> is totally singular subspace wrt <form>.
##
InstallMethod(IsTotallySingularSubspace,
  "for a quadratic form and a basis of a subspace",
  [IsQuadraticForm, IsMatrix],
  function(f,sub)
  local fsub;
  fsub := Filtered(sub,x-> IsZero(x^f));
  if Length(fsub) <> Length(sub) then
    return false;
  else
    return IsTotallyIsotropicSubspace(AssociatedBilinearForm(f),sub);
  fi;
end );

#############################################################################
#A TypeOfForm( <form> ) <form>: sesquilinear or quadratic form
#  returns one of 0, 1, -1, 1/2, or -1/2.
# Let <f> be a form on V(n,q), with radical R, a k-dimensional subspace
# of V(n,q), 0 <= k <= n. Then <f> induces a non-degenerate/non-singular
# form <g> on V/R. When dim(R)=0, this induced form <g> is <f> itself of course.
# TypeOfForm(f) returns:
#  - 0 when g is symplectic (1) or parabolic (2)
#  - +1 when g is hyperbolic (3)
#  - -1 when g is elliptic (4)
#  - -1/2 when g is hermitian in odd dimension (5)
#  - 1/2 when g is hermitian in even dimension (6)
#  - An error when f is a pseudo form (7);
#    note that no method is installed for trivial forms.
#  One way to remember it is that the number of points of the
#  corresponding rank 1 polar space is q^(1 - sign) + 1, q the order
#  of the base field of the form.
#  The above description is valid for sesquilinear and quadratic forms,
#  of course case (1) can only occur for bilinear forms, cases (2), (3) and (4)
#  for both quadratic and bilinear forms, cases (5) and (6) only for hermitian
#  forms and case (7) only for bilinear forms.
#
#  In principle, for the non-degenerate/non-singular orthogonal forms,
#  it could be sufficient to investigate whether the determinant is square or not,
#  however, since degenerate/singular forms are allowed, in those cases, we must
#  compute the induced form to use the determinant method. Taking into account
#  that all is already available by simply testing for IsEllipticForm, IsHyperbolicForm
#  IsParabolicForm (which goes through the base change mechanisms), I opted to use
#  simply these tests. The hermitian case is straightforwardly done.
##
#############################################################################
#A TypeOfForm( <form> ) <form>: bilinear form
#
InstallMethod( TypeOfForm,
    "for a bilinear form",
    [IsBilinearForm],
    function(f)
    if IsSymplecticForm(f)
        then return 0;
    elif IsEllipticForm(f)
        then return -1;
    elif IsHyperbolicForm(f)
        then return 1;
    elif IsParabolicForm(f)
        then return 0;
    else
        Error("<f> is a pseudo form and has no defined type");
    fi;
end );

#############################################################################
#A TypeOfForm( <form> ) <form>: quadratic form
#
InstallMethod( TypeOfForm,
    "for a quadratic form",
    [IsQuadraticForm],
    function(f)
    if IsEllipticForm(f)
        then return -1;
    elif IsHyperbolicForm(f)
        then return 1;
    elif IsParabolicForm(f)
        then return 0;
    fi;
  end );

#############################################################################
#A TypeOfForm( <form> ) <form>: hermitian form
#
InstallMethod( TypeOfForm,
  "for a unitary form",
  [IsHermitianForm],
  function(f)
  local m,dim;
    m := f!.matrix;
    dim := Dimension(RadicalOfForm(f));
    if IsOddInt(NrRows(m)-dim) then
       return -1/2;
    else
       return 1/2;
    fi;
  end );


