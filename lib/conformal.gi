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
                     Forms_FieldOfDefinition( G, InvariantBilinearForm( G ) ) ) ) );

InstallMethod( ConformalSymplecticGroupCons,
  "matrix group for group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsGroup and HasInvariantBilinearFormUpToScalars" ],
  { filt, G } -> ConformalSymplecticGroupCons( filt,
                   BilinearFormByMatrix(
                     InvariantBilinearFormUpToScalars( G ).matrix,
                     Forms_FieldOfDefinition( G, InvariantBilinearFormUpToScalars( G ) ) ) ) );

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
    rec( matrix:= Matrix( form!.matrix, stored ), baseDomain:= F ) );

  if HasIsFullSubgroupGLRespectingBilinearFormUpToScalars( g ) then
    SetIsFullSubgroupGLRespectingBilinearFormUpToScalars( gg,
        IsFullSubgroupGLRespectingBilinearFormUpToScalars( g ) );
  fi;

  return gg;
end );


#############################################################################
##
##  The following methods are currently needed to make the code work
##  in case one creates groups whose elements are in `IsMatrixObj`.
##  Eventually we must support `IsMatrixObj` matrices in form objects.
##

InstallOtherMethod( BilinearFormByMatrix,
  "for a ffe matrix object and a field",
  [ "IsMatrixObj and IsFFECollColl", "IsField and IsFinite" ],
  { m, F } -> BilinearFormByMatrix( Unpack( m ), F ) );

InstallOtherMethod( BilinearFormByMatrix,
  "for a ffe matrix object",
  [ "IsMatrixObj and IsFFECollColl" ],
  m -> BilinearFormByMatrix( Unpack( m ) ) );


# The following is apparently needed in the tests in `tst/adv/conformal.tst`.
# Strictly speaking, the following is not correct,
# according to the definition of `DegreeFFE`.
# Eventually we should fix the use of `FieldOfMatrixList`
# and `FieldOfMatrixGroup`, then `DegreeFFE` will not be important anymore.
InstallOtherMethod( DegreeFFE,
  [ "IsMatrixObj and IsFFECollColl" ],
  mat -> DegreeOverPrimeField( BaseDomain( mat ) ) );
