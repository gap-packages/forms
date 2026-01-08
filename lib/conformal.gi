
# auxiliary function analogous to `DescribesInvariantBilinearForm`
BindGlobal( "DescribesInvariantBilinearFormUpToScalars",
    obj -> IsMatrixOrMatrixObj( obj ) or
           IsBilinearForm( obj ) or
           ( IsGroup( obj ) and HasInvariantBilinearForm( obj ) ) or
           ( IsGroup( obj ) and HasInvariantBilinearFormUpToScalars( obj ) ) );


#############################################################################
##
#F  ConformalSymplecticGroup( [<filt>, ]<d>, <q>[, <form>] ) conf. sympl. gp.
#F  ConformalSymplecticGroup( [<filt>, ]<d>, <R>[, <form>] ) conf. sympl. gp.
#F  ConformalSymplecticGroup( [<filt>, ]<form> )   conformal symplectic group
#F  CSp( [<filt>, ]<d>, <q>[, <form>] )            conformal symplectic group
#F  CSp( [<filt>, ]<d>, <R>[, <form>] )            conformal symplectic group
#F  CSp( [<filt>, ]<form> )                        conformal symplectic group
##
InstallGlobalFunction( ConformalSymplecticGroup, function ( arg )
  local filt, form;

  if IsFilter( First( arg ) ) then
    filt:= Remove( arg, 1 );
  else
    filt:= IsMatrixGroup;
  fi;
  if DescribesInvariantBilinearFormUpToScalars( Last( arg ) ) then
    # interpret this argument (matrix or form or group with stored form)
    # as "up to scalars"
    form:= Remove( arg );
    if Length( arg ) = 0 then
      # ( [<filt>, ]<form> )
      return ConformalSymplecticGroupCons( filt, form );
    elif Length( arg ) = 2 and IsPosInt( arg[1] )
                           and ( IsRing( arg[2] ) or IsPosInt( arg[2] ) ) then
      # ( [<filt>, ]<d>, <R>, <form> ) or ( [<filt>, ]<d>, <q>, <form> )
      return ConformalSymplecticGroupCons( filt, arg[1], arg[2], form );
    fi;
  elif Length( arg ) = 2 and IsPosInt( arg[1] )
                         and ( IsRing( arg[2] ) or IsPosInt( arg[2] ) ) then
    # ( [<filt>, ]<d>, <R> ) or ( [<filt>, ]<d>, <q> )
    return ConformalSymplecticGroupCons( filt, arg[1], arg[2] );
  fi;
  Error( "usage: ConformalSymplecticGroup( [<filt>, ]<d>, <R>[, <form>] )\n",
         "or ConformalSymplecticGroup( [<filt>, ]<d>, <q>[, <form>] )\n",
         "or ConformalSymplecticGroup( [<filt>, ]<form> )" );
end );


#############################################################################
##
#M  ConformalSymplecticGroupCons( <IsMatrixGroup>, <d>, <F> )
##
InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension and finite field",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsField and IsFinite" ],
  function( filter, d, F )
  local q, z, o, filt, g, mat1, mat2, i, mat3, size, qi, c;

  # the dimension must be even
  if d mod 2 = 1 then
    Error( "the dimension <d> must be even" );
  fi;
  q:= Size( F );
  z:= PrimitiveRoot( F );
  o:= One( F );

  # Decide about the internal representation of group generators.
  filt:= ValueOption( "ConstructingFilter" );
  if filt = fail then
    filt:= IsPlistRep;
  fi;

  # if the dimension is two it is a general linear group
  if d = 2 then
    g:= GL( 2, F );
#T TODO/Note:
#T Currently `filt` is ignored here.
#T This will be fixed automatically as soon as
#T `GL` also supports the global option.
    c:= List( GeneratorsOfGroup( g ), x -> Matrix( filt, F, x ) );
    c:= GroupWithGenerators( c );
    SetDimensionOfMatrixGroup( c, d );
    SetFieldOfMatrixGroup( c, F );
    SetName( c, Name( g ) );
    SetSize( c, Size( g ) );
    g:= c;
  else
      # CSp(4,2)
      if d = 4 and q = 2  then
        mat1:= Matrix( filt, F, [1,0,1,1,1,0,0,1,0,1,0,1,1,1,1,1] * o, 4 );
        mat2:= Matrix( filt, F, [0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0] * o, 4 );

      # CSp(d,q)
      else
        mat1 := IdentityMatrix( filt, F, d );
        mat2 := ZeroMatrix( filt, F, d, d );
        for i  in [ 2 .. d/2 ]      do mat2[i,i-1]:= o;  od;
        for i  in [ d/2+1 .. d-1 ]  do mat2[i,i+1]:= o;  od;

      if q mod 2 = 1  then
        mat1[  1,    1] := z;
        mat1[  d,    d] := z^-1;
        mat2[  1,    1] := o;
        mat2[  1,d/2+1] := o;
        mat2[d-1,  d/2] := o;
        mat2[  d,  d/2] := -o;

      elif q <> 2  then
        mat1[    1,    1] := z;
        mat1[  d/2,  d/2] := z;
        mat1[d/2+1,d/2+1] := z^-1;
        mat1[    d,    d] := z^-1;
        mat2[    1,d/2-1] := o;
        mat2[    1,  d/2] := o;
        mat2[    1,d/2+1] := o;
        mat2[d/2+1,  d/2] := o;
        mat2[    d,  d/2] := o;

      else
        mat1[    1,  d/2] := o;
        mat1[    1,    d] := o;
        mat1[d/2+1,    d] := o;
        mat2[    1,d/2+1] := o;
        mat2[    d,  d/2] := o;
      fi;
    fi;

    mat3:= IdentityMatrix( filt, F, d );
    for i in [ 1 .. d/2 ] do
      mat3[i, i]:= z;
    od;

    mat1:= ImmutableMatrix( F, mat1, true );
    mat2:= ImmutableMatrix( F, mat2, true );
    mat3:= ImmutableMatrix( F, mat3, true );

    # avoid to call 'Group' because this would check invertibility ...
    g:= GroupWithGenerators( [ mat1, mat2, mat3 ] );
    SetName( g, Concatenation( "CSp(", String(d), ",", String(q), ")" ) );
    SetDimensionOfMatrixGroup( g, d );
    SetFieldOfMatrixGroup( g, F );

    # add the size
    size := 1;
    qi   := 1;
    for i in [ 1 .. d/2 ] do
      qi   := qi * q^2;
      size := size * (qi-1);
    od;
    SetSize( g, q^((d/2)^2) * size * (q-1) );
  fi;

  # construct the form
  c:= ZeroMatrix( filt, F, d, d );
  for i in [ 1 .. d/2 ] do
    c[i,d-i+1]:= o;
    c[d/2+i,d/2-i+1]:= -o;
  od;
  SetInvariantBilinearFormUpToScalars( g,
      rec( matrix:= ImmutableMatrix( F, c, true ) ) );

  SetIsFullSubgroupGLRespectingBilinearFormUpToScalars( g, true );

  # and return
  return g;
end );


#############################################################################
##
#M  ConformalSymplecticGroupCons( <IsMatrixGroup>, <d>, <q> )
##
InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension and finite field size",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt" ],
  { filt, n, q } -> ConformalSymplecticGroupCons( filt, n, GF(q) ) );


#############################################################################
##
#M  ConformalSymplecticGroupCons( <filt>, <form> )
##
InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for matrix of form",
  [ "IsMatrixGroup and IsFinite", "IsMatrixOrMatrixObj" ],
  { filt, mat } -> ConformalSymplecticGroupCons( filt,
                     BilinearFormByMatrix( mat, BaseDomain( mat ) ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for group with form",
  [ "IsMatrixGroup and IsFinite", "IsGroup and HasInvariantBilinearForm" ],
  { filt, G } -> ConformalSymplecticGroupCons( filt,
                   BilinearFormByMatrix(
                     InvariantBilinearForm( G ).matrix,
                     FieldOfMatrixGroup( G ) ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsGroup and HasInvariantBilinearFormUpToScalars" ],
  { filt, G } -> ConformalSymplecticGroupCons( filt,
                   BilinearFormByMatrix(
                     InvariantBilinearFormUpToScalars( G ).matrix,
                     FieldOfMatrixGroup( G ) ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for form",
  [ "IsMatrixGroup and IsFinite", "IsBilinearForm" ],
  { filt, form } -> ConformalSymplecticGroupCons( filt,
                      NumberRows( form!.matrix ),
                      form!.basefield, form ) );


#############################################################################
##
#M  ConformalSymplecticGroupCons( <filt>, <d>, <q>, <form> )
##
InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field size, matrix of form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsMatrixOrMatrixObj" ],
  { filt, d, q, mat } -> ConformalSymplecticGroupCons( filt, d, GF(q),
                           BilinearFormByMatrix( mat, GF(q) ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field size, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsGroup and HasInvariantBilinearForm" ],
  { filt, d, q, G } -> ConformalSymplecticGroupCons( filt, d, GF(q),
                         BilinearFormByMatrix(
                           InvariantBilinearForm( G ).matrix, GF(q) ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field size, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsGroup and HasInvariantBilinearFormUpToScalars" ],
  { filt, d, q, G } -> ConformalSymplecticGroupCons( filt, d, GF(q),
                         BilinearFormByMatrix(
                           InvariantBilinearFormUpToScalars( G ).matrix, GF(q) ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field size, form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsBilinearForm" ],
  { filt, d, q, form } -> ConformalSymplecticGroupCons( filt, d, GF(q), form ) );


#############################################################################
##
#M  ConformalSymplecticGroupCons( <filt>, <d>, <R>, <form> )
##
InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field, matrix of form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsField and IsFinite",
    "IsMatrixOrMatrixObj" ],
  { filt, d, F, form } -> ConformalSymplecticGroupCons( filt, d, F,
                            BilinearFormByMatrix( form, F ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsField and IsFinite",
    "IsGroup and HasInvariantBilinearForm" ],
  { filt, d, F, G } -> ConformalSymplecticGroupCons( filt, d, F,
                         BilinearFormByMatrix(
                           InvariantBilinearForm( G ).matrix, F ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsField and IsFinite",
    "IsGroup and HasInvariantBilinearFormUpToScalars" ],
  { filt, d, F, G } -> ConformalSymplecticGroupCons( filt, d, F,
                         BilinearFormByMatrix(
                           InvariantBilinearFormUpToScalars( G ).matrix, F ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for dimension, finite field, form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsField and IsFinite",
    "IsBilinearForm" ],
  function( filt, d, F, form )
  local g, stored, wanted, mat1, mat2, mat, matinv, gens, gg;

  # Create the default generators and form.
  g:= ConformalSymplecticGroupCons( filt, d, F );
  stored:= InvariantBilinearFormUpToScalars( g ).matrix;

  # If the prescribed form fits then just return.
  if stored = form!.matrix then
    return g;
  fi;

  # Compute a base change matrix.
  # (Check that the canonical forms are equal.)
  wanted:= BilinearFormByMatrix( stored, F );
  mat1:= BaseChangeToCanonical( form );
  mat2:= BaseChangeToCanonical( wanted );
  if mat1 * form!.matrix * TransposedMat( mat1 ) <>
     mat2 * stored * TransposedMat( mat2 ) then
    Error( "canonical forms of <form> and <wanted> differ" );
  fi;
  mat:= mat2^-1 * mat1;
  matinv:= mat^-1;

  # Create the group w.r.t. the prescribed form.
  gens:= List( GeneratorsOfGroup( g ),
               x -> Matrix( matinv * x * mat, stored ) );
  gg:= GroupWithGenerators( gens );

  UseIsomorphismRelation( g, gg );

  if HasName( g ) then
    SetName( gg, Name( g ) );
  fi;

  SetInvariantBilinearFormUpToScalars( gg,
    rec( matrix:= Matrix( form!.matrix, stored ) ) );

  if HasIsFullSubgroupGLRespectingBilinearFormUpToScalars( g ) then
    SetIsFullSubgroupGLRespectingBilinearFormUpToScalars( gg,
        IsFullSubgroupGLRespectingBilinearFormUpToScalars( g ) );
  fi;

  return gg;
end );


#############################################################################
##
##  Support `IsPermGroup` as first argument in `ConformalSymplecticGroup`.
##
PermConstructor( ConformalSymplecticGroupCons,
  [ IsPermGroup, IsInt, IsObject ],
  IsMatrixGroup and IsFinite);


#############################################################################
##
#M  \in( <mat>, <CSp> ) . . . . . . . .  membership test method based on form
##
InstallMethod( \in, "respecting bilinear form", IsElmsColls,
  [ "IsMatrixOrMatrixObj",
  "IsFullSubgroupGLRespectingBilinearFormUpToScalars" ],
  {} -> RankFilter( IsHandledByNiceMonomorphism ), # override nice mon. method
  function( mat, G )
  local inv;

  if not IsSubset( FieldOfMatrixGroup( G ),
                   FieldOfMatrixList( [ mat ] ) ) then
    return false;
  fi;
  inv:= InvariantBilinearFormUpToScalars( G ).matrix;
  return _IsEqualModScalars( inv, mat * inv * TransposedMat( mat ) );
end );


#############################################################################
##
##  The following methods are currently needed to make the code work
##  in case one creates groups whose elements are in `IsMatrixObj`.
##  Eventually we want to get rid of them (or or to add something similar
##  to the GAP library).
##

# Strictly speaking, the following is not correct,
# according to the definition of `DegreeFFE`.
# Eventually we should fix the use of `FieldOfMatrixList`
# and `FieldOfMatrixGroup`, then `DegreeFFE` will not be important anymore.
InstallOtherMethod( DegreeFFE,
  [ "IsMatrixObj and IsFFECollColl" ],
  mat -> DegreeOverPrimeField( BaseDomain( mat ) ) );

# This is a really ugly hack.
# Without this method, multiplying a matrix object <M1> with a matrix <M2>
# is possible but yields the matrix product of <M1> with the *transposed*
# of <M2>.
# (This is because GAP regards <M1> as a scalar and computes the list of
# products of <M1> with the rows of <M2>.)
InstallOtherMethod( \*,
  [ "IsMatrixObj", "IsMatrix" ],
  { matobj, mat } -> Unpack( matobj ) * mat );


InstallOtherMethod( BilinearFormByMatrix,
  "for a ffe matrix object and a field",
  [ "IsMatrixObj and IsFFECollColl", "IsField and IsFinite" ],
  { m, F } -> BilinearFormByMatrix( Unpack( m ), F ) );

InstallOtherMethod( BilinearFormByMatrix,
  "for a ffe matrix object",
  [ "IsMatrixObj and IsFFECollColl" ],
  m -> BilinearFormByMatrix( Unpack( m ) ) );

