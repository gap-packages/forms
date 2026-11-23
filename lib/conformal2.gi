#############################################################################
##
##  conformal2.gi         'Forms' package
##
##  Provide the methods for 'ConformalOrthogonalGroupCons' and
##  'ConformalUnitaryGroupCons' that involve a prescribed form.
##


#############################################################################
##
#F  Forms_ConformalOrthogonalGroup( <g>, <form> )
##
##  Return a conjugate of the conformal orthogonal group <g> that respects
##  the quadratic form <form> up to scalars.
##
BindGlobal( "Forms_ConformalOrthogonalGroup", function( g, form )
    local stored, gf, wanted, mat1, mat2, mat, matinv, gens, gg;

    stored:= InvariantQuadraticFormUpToScalars( g );
    gf:= stored.baseDomain;
    stored:= stored.matrix;

    # Check that 'form' lives over the intended field.
    if not IsSubset( gf, form!.basefield ) then
      Error( "the defining field of <form> does not fit to <gf>" );
    fi;

    # If the prescribed form fits then just return.
    if stored = form!.matrix then
      return g;
    fi;

    # Compute a base change matrix.
    # (Check that the canonical forms are equal up to scalars.)
    wanted:= QuadraticFormByMatrix( stored, gf );
    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( wanted );
    if not _IsEqualModScalars(
               Forms_RESET( mat1 * form!.matrix * TransposedMat( mat1 ) ),
               Forms_RESET( mat2 * stored * TransposedMat( mat2 ) ) ) then
      Error( "canonical forms of <form> and <wanted> differ" );
    fi;
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    # Create the group w.r.t. the prescribed form.
    gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
    gg:= GroupWithGenerators( gens );

    UseIsomorphismRelation( g, gg );

    if HasName( g ) then
      SetName( gg, Name( g ) );
    fi;

    SetInvariantQuadraticFormUpToScalars( gg,
        rec( matrix:= form!.matrix, baseDomain:= gf ) );
    if HasIsFullSubgroupGLRespectingQuadraticFormUpToScalars( g ) then
      SetIsFullSubgroupGLRespectingQuadraticFormUpToScalars( gg,
          IsFullSubgroupGLRespectingQuadraticFormUpToScalars( g ) );
    fi;

    SetInvariantBilinearFormUpToScalars( gg,
        rec( matrix:= form!.matrix + TransposedMat( form!.matrix ),
             baseDomain:= gf ) );
    if Characteristic( gf ) <> 2 and
       HasIsFullSubgroupGLRespectingBilinearFormUpToScalars( g ) then
      SetIsFullSubgroupGLRespectingBilinearFormUpToScalars( gg,
          IsFullSubgroupGLRespectingBilinearFormUpToScalars( g ) );
    fi;

    return gg;
end );


#############################################################################
##
#M  ConformalOrthogonalGroupCons( <filt>, <form> )
##
InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for matrix of form",
  [ "IsMatrixGroup and IsFinite", "IsMatrixOrMatrixObj" ],
  { filt, mat } -> ConformalOrthogonalGroupCons( filt,
                     QuadraticFormByMatrix( mat, BaseDomain( mat ) ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for group with form",
  [ "IsMatrixGroup and IsFinite", "IsGroup and HasInvariantQuadraticForm" ],
  { filt, G } -> ConformalOrthogonalGroupCons( filt,
                   QuadraticFormByMatrix(
                     InvariantQuadraticForm( G ).matrix,
                     Forms_FieldOfDefinition( InvariantQuadraticForm( G ), G ) ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsGroup and HasInvariantQuadraticFormUpToScalars" ],
  { filt, G } -> ConformalOrthogonalGroupCons( filt,
                   QuadraticFormByMatrix(
                     InvariantQuadraticFormUpToScalars( G ).matrix,
                     InvariantQuadraticFormUpToScalars( G ).baseDomain ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for form",
  [ "IsMatrixGroup and IsFinite", "IsQuadraticForm" ],
  function( filt, form )
  local d, q, e;

  d:= NumberRows( form!.matrix );
  q:= Size( form!.basefield );
  if IsOddInt( d ) then
    e:= 0;
  elif IsEllipticForm( form ) then
    e:= -1;
  else
    e:= 1;
  fi;
  return ConformalOrthogonalGroupCons( filt, e, d, q, form );
end );


#############################################################################
##
#M  ConformalOrthogonalGroupCons( <filt>, <e>, <d>, <q>, <form> )
##
InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field size, matrix of form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsPosInt",
    "IsMatrixOrMatrixObj" ],
  { filt, e, d, q, mat } -> ConformalOrthogonalGroupCons( filt, e, d, q,
                              QuadraticFormByMatrix( mat, GF(q) ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field size, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsPosInt",
    "IsGroup and HasInvariantQuadraticForm" ],
  { filt, e, d, q, G } -> ConformalOrthogonalGroupCons( filt, e, d, q,
                            QuadraticFormByMatrix(
                              InvariantQuadraticForm( G ).matrix, GF(q) ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field size, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsPosInt",
    "IsGroup and HasInvariantQuadraticFormUpToScalars" ],
  { filt, e, d, q, G } -> ConformalOrthogonalGroupCons( filt, e, d, q,
                            QuadraticFormByMatrix(
                              InvariantQuadraticFormUpToScalars( G ).matrix,
                              GF(q) ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field size, form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsPosInt",
    "IsQuadraticForm" ],
  { filt, e, d, q, form } -> Forms_ConformalOrthogonalGroup(
                               ConformalOrthogonalGroupCons( filt, e, d, q ),
                               form ) );


#############################################################################
##
#M  ConformalOrthogonalGroupCons( <filt>, <e>, <d>, <R>, <form> )
##
InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field, matrix of form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsField and IsFinite",
    "IsMatrixOrMatrixObj" ],
  { filt, e, d, F, mat } -> ConformalOrthogonalGroupCons( filt, e, d,
                              Size( F ), QuadraticFormByMatrix( mat, F ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsField and IsFinite",
    "IsGroup and HasInvariantQuadraticForm" ],
  { filt, e, d, F, G } -> ConformalOrthogonalGroupCons( filt, e, d,
                            Size( F ),
                            QuadraticFormByMatrix(
                              InvariantQuadraticForm( G ).matrix, F ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsField and IsFinite",
    "IsGroup and HasInvariantQuadraticFormUpToScalars" ],
  { filt, e, d, F, G } -> ConformalOrthogonalGroupCons( filt, e, d,
                            Size( F ),
                            QuadraticFormByMatrix(
                              InvariantQuadraticFormUpToScalars( G ).matrix,
                              F ) ) );

InstallMethod( ConformalOrthogonalGroupCons,
  "matrix group for <e>, dimension, finite field, form",
  [ "IsMatrixGroup and IsFinite",
    "IsInt",
    "IsPosInt",
    "IsField and IsFinite",
    "IsQuadraticForm" ],
  { filt, e, d, F, form } -> ConformalOrthogonalGroupCons( filt, e, d,
                               Size( F ), form ) );


#############################################################################
##
#M  ConformalUnitaryGroupCons( <filt>, <form> )
##
InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for matrix of form",
  [ "IsMatrixGroup and IsFinite", "IsMatrixOrMatrixObj" ],
  function( filt, mat )
  local q, form;

  # Guess the intended field.
  q:= Size( BaseDomain( mat ) );
  if q = RootInt( q, 2 )^2 then
    form:= HermitianFormByMatrix( mat, GF(q) );
  else
    form:= HermitianFormByMatrix( mat, GF(q^2) );
  fi;
  return ConformalUnitaryGroupCons( filt, form );
end );

InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsGroup and HasInvariantSesquilinearForm" ],
  { filt, G } -> ConformalUnitaryGroupCons( filt,
                   HermitianFormByMatrix(
                     InvariantSesquilinearForm( G ).matrix,
                     Forms_FieldOfDefinition( InvariantSesquilinearForm( G ), G ) ) ) );

InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsGroup and HasInvariantSesquilinearFormUpToScalars" ],
  { filt, G } -> ConformalUnitaryGroupCons( filt,
                   HermitianFormByMatrix(
                     InvariantSesquilinearFormUpToScalars( G ).matrix,
                     InvariantSesquilinearFormUpToScalars( G ).baseDomain ) ) );

InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for form",
  [ "IsMatrixGroup and IsFinite", "IsHermitianForm" ],
  { filt, form } -> ConformalUnitaryGroupCons( filt,
                      NumberRows( form!.matrix ),
                      RootInt( Size( form!.basefield ), 2 ), form ) );


#############################################################################
##
#M  ConformalUnitaryGroupCons( <filt>, <d>, <q>, <form> )
##
InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for dimension, finite field size, matrix of form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsMatrixOrMatrixObj" ],
  { filt, d, q, mat } -> ConformalUnitaryGroupCons( filt, d, q,
                           HermitianFormByMatrix( mat, GF(q^2) ) ) );

InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for dimension, finite field size, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsGroup and HasInvariantSesquilinearForm" ],
  { filt, d, q, G } -> ConformalUnitaryGroupCons( filt, d, q,
                         HermitianFormByMatrix(
                           InvariantSesquilinearForm( G ).matrix, GF(q^2) ) ) );

InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for dimension, finite field size, group with form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsGroup and HasInvariantSesquilinearFormUpToScalars" ],
  { filt, d, q, G } -> ConformalUnitaryGroupCons( filt, d, q,
                         HermitianFormByMatrix(
                           InvariantSesquilinearFormUpToScalars( G ).matrix,
                           GF(q^2) ) ) );

InstallMethod( ConformalUnitaryGroupCons,
  "matrix group for dimension, finite field size, form",
  [ "IsMatrixGroup and IsFinite",
    "IsPosInt",
    "IsPosInt",
    "IsHermitianForm" ],
  function( filt, d, q, form )
  local g, stored, F, wanted, mat1, mat2, mat, matinv, gens, gg;

  # Create the default generators and form.
  g:= ConformalUnitaryGroupCons( filt, d, q );
  stored:= InvariantSesquilinearFormUpToScalars( g ).matrix;

  # Check that 'form' lives over the intended field.
  F:= GF(q^2);
  if not IsSubset( F, form!.basefield ) then
    Error( "the defining field of <form> does not fit to <q>" );
  fi;

  # If the prescribed form fits then just return.
  if stored = form!.matrix then
    return g;
  fi;

  # Compute a base change matrix.
  # (Check that the canonical forms are equal up to scalars.)
  wanted:= HermitianFormByMatrix( stored, F );
  mat1:= BaseChangeToCanonical( form );
  mat2:= BaseChangeToCanonical( wanted );
  if not _IsEqualModScalars(
             mat1 * form!.matrix * TransposedFrobeniusMat( mat1, q ),
             mat2 * stored * TransposedFrobeniusMat( mat2, q ) ) then
    Error( "canonical forms of <form> and <wanted> differ" );
  fi;
  mat:= mat2^-1 * mat1;
  matinv:= mat^-1;

  # Create the group w.r.t. the prescribed form.
  gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
  gg:= GroupWithGenerators( gens );

  UseIsomorphismRelation( g, gg );

  if HasName( g ) then
    SetName( gg, Name( g ) );
  fi;

  SetInvariantSesquilinearFormUpToScalars( gg,
      rec( matrix:= form!.matrix, baseDomain:= F ) );
  if HasIsFullSubgroupGLRespectingSesquilinearFormUpToScalars( g ) then
    SetIsFullSubgroupGLRespectingSesquilinearFormUpToScalars( gg,
        IsFullSubgroupGLRespectingSesquilinearFormUpToScalars( g ) );
  fi;

  return gg;
end );
