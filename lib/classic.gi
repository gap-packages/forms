#############################################################################
##
##  classic.gi            'Forms' package
##
##  Provide the methods for 'GO', 'SO', 'Omega', 'GU', 'SU', and 'Sp' that
##  involve a prescribed invariant form.
##


# Compatibility with GAP < 4.12
if not IsBound(IsMatrixOrMatrixObj) then
    BindGlobal("IsMatrixOrMatrixObj", IsMatrixObj);
fi;

#############################################################################
##
#O  GeneralOrthogonalGroupCons( <filter>, <form> )
#O  GeneralOrthogonalGroupCons( <filter>, <e>, <d>, <q>, <form> )
#O  GeneralOrthogonalGroupCons( <filter>, <e>, <d>, <R>, <form> )
##
##  'GeneralOrthogonalGroup' is a plain function that is defined in the GAP
##  library.
##  It calls 'GeneralOrthogonalGroupCons',
##  thus we have to declare the variants involving a quadratic form,
##  and install the corresponding methods.
##
Perform(
    [ IsMatrixOrMatrixObj, IsQuadraticForm, IsGroup and HasInvariantQuadraticForm ],
    function( obj )
      DeclareConstructor( "GeneralOrthogonalGroupCons",
        [ IsGroup, obj ] );
      DeclareConstructor( "GeneralOrthogonalGroupCons",
        [ IsGroup, IsInt, IsPosInt, IsPosInt, obj ] );
      DeclareConstructor( "GeneralOrthogonalGroupCons",
        [ IsGroup, IsInt, IsPosInt, IsRing, obj ] );
    end );


#############################################################################
##
#M  GeneralOrthogonalGroupCons( <filt>, <form> )
##
InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for matrix of form",
    [ IsMatrixGroup and IsFinite, IsMatrixOrMatrixObj ],
    { filt, mat } -> GeneralOrthogonalGroupCons( filt,
                       QuadraticFormByMatrix( mat, BaseDomain( mat ) ) ) );

InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for group with form",
    [ IsMatrixGroup and IsFinite, IsGroup and HasInvariantQuadraticForm ],
    { filt, G } -> GeneralOrthogonalGroupCons( filt,
                     QuadraticFormByMatrix(
                       InvariantQuadraticForm( G ).matrix,
                       FieldOfMatrixGroup( G ) ) ) );

InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for form",
    [ IsMatrixGroup and IsFinite, IsQuadraticForm ],
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
    return GeneralOrthogonalGroupCons( filt, e, d, q, form );
end );


#############################################################################
##
#M  GeneralOrthogonalGroupCons( <filt>, <e>, <d>, <q>, <form> )
##
InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field size, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsMatrixOrMatrixObj ],
    { filt, e, d, q, mat } -> GeneralOrthogonalGroupCons( filt, e, d, q,
                                QuadraticFormByMatrix( mat, GF(q) ) ) );

InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field size, group with form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsGroup and HasInvariantQuadraticForm ],
    { filt, e, d, q, G } -> GeneralOrthogonalGroupCons( filt, e, d, q,
                              QuadraticFormByMatrix(
                                InvariantQuadraticForm( G ).matrix, GF(q) ) ) );

InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field size, form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsQuadraticForm ],
    function( filt, e, d, q, form )
    local g, stored, wanted, mat1, mat2, mat, matinv, gens, gg;

    # Create the default generators and form.
    g:= GeneralOrthogonalGroupCons( filt, e, d, q );
    stored:= InvariantQuadraticForm( g ).matrix;

    # If the prescribed form fits then just return.
    if stored = form!.matrix then
      return g;
    fi;

    # Compute a base change matrix.
    # (Check that the canonical forms are equal.)
    wanted:= QuadraticFormByMatrix( stored, GF(q) );
    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( wanted );
    if Forms_RESET( mat1 * form!.matrix * TransposedMat( mat1 ), d, q ) <>
       Forms_RESET( mat2 * stored * TransposedMat( mat2 ), d, q ) then
      Error( "canonical forms of <form> and <wanted> differ" );
    fi;
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( QuadraticFormByMatrix( stored, GF(q) ) );
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    # Create the group w.r.t. the prescribed form.
    gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
    gg:= GroupWithGenerators( gens );

    if HasSize( g ) then
      SetSize( gg, Size( g ) );
    fi;

    if HasName( g ) then
      SetName( gg, Name( g ) );
    fi;

    SetInvariantQuadraticForm( gg, rec( matrix:= form!.matrix ) );
    if HasIsFullSubgroupGLorSLRespectingQuadraticForm( g ) then
      SetIsFullSubgroupGLorSLRespectingQuadraticForm( gg,
          IsFullSubgroupGLorSLRespectingQuadraticForm( g ) );
    fi;

    SetInvariantBilinearForm( gg, rec( matrix:= matinv * InvariantBilinearForm( g ).matrix * TransposedMat( matinv ) ) );
    if q mod 2 = 1 and
       HasIsFullSubgroupGLorSLRespectingBilinearForm( g ) then
      SetIsFullSubgroupGLorSLRespectingBilinearForm( gg,
          IsFullSubgroupGLorSLRespectingBilinearForm( g ) );
    fi;

    return gg;
end );


#############################################################################
##
#M  GeneralOrthogonalGroupCons( <filt>, <e>, <d>, <R>, <form> )
##
InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsField and IsFinite,
      IsMatrixOrMatrixObj ],
    { filt, e, d, F, form } -> GeneralOrthogonalGroupCons( filt, e, d,
                                 Size( F ),
                                 QuadraticFormByMatrix( form, F ) ) );

InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field, group with form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsField and IsFinite,
      IsGroup and HasInvariantQuadraticForm ],
    { filt, e, d, F, G } -> GeneralOrthogonalGroupCons( filt, e, d,
                              Size( F ),
                              QuadraticFormByMatrix(
                                InvariantQuadraticForm( G ).matrix, F ) ) );

InstallMethod( GeneralOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field, form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsField and IsFinite,
      IsQuadraticForm ],
    { filt, e, d, F, form } -> GeneralOrthogonalGroupCons( filt, e, d,
                                 Size( F ), form ) );


#############################################################################
##
#O  SpecialOrthogonalGroupCons( <filter>, <form> )
#O  SpecialOrthogonalGroupCons( <filter>, <e>, <d>, <q>, <form> )
#O  SpecialOrthogonalGroupCons( <filter>, <e>, <d>, <R>, <form> )
##
##  'SpecialOrthogonalGroup' is a plain function that is defined in the GAP
##  library.
##  It calls 'SpecialOrthogonalGroupCons',
##  thus we have to declare the variants involving a quadratic form,
##  and install the corresponding methods.
##
Perform(
    [ IsMatrixOrMatrixObj, IsQuadraticForm, IsGroup and HasInvariantQuadraticForm ],
    function( obj )
      DeclareConstructor( "SpecialOrthogonalGroupCons",
        [ IsGroup, obj ] );
      DeclareConstructor( "SpecialOrthogonalGroupCons",
        [ IsGroup, IsInt, IsPosInt, IsPosInt, obj ] );
      DeclareConstructor( "SpecialOrthogonalGroupCons",
        [ IsGroup, IsInt, IsPosInt, IsRing, obj ] );
    end );


#############################################################################
##
#M  SpecialOrthogonalGroupCons( <filt>, <form> )
##
InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for matrix of form",
    [ IsMatrixGroup and IsFinite, IsMatrixOrMatrixObj ],
    { filt, mat } -> SpecialOrthogonalGroupCons( filt,
                       QuadraticFormByMatrix( mat, BaseDomain( mat ) ) ) );

InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for group with form",
    [ IsMatrixGroup and IsFinite, IsGroup and HasInvariantQuadraticForm ],
    { filt, G } -> SpecialOrthogonalGroupCons( filt,
                     QuadraticFormByMatrix(
                       InvariantQuadraticForm( G ).matrix,
                       FieldOfMatrixGroup( G ) ) ) );

InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for form",
    [ IsMatrixGroup and IsFinite, IsQuadraticForm ],
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
    return SpecialOrthogonalGroupCons( filt, e, d, q, form );
end );


#############################################################################
##
#M  SpecialOrthogonalGroupCons( <filt>, <e>, <d>, <q>, <form> )
##
InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field size, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsMatrixOrMatrixObj ],
    { filt, e, d, q, mat } -> SpecialOrthogonalGroupCons( filt, e, d, q,
                                QuadraticFormByMatrix( mat, GF(q) ) ) );

InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field size, group with form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsGroup and HasInvariantQuadraticForm ],
    { filt, e, d, q, G } -> SpecialOrthogonalGroupCons( filt, e, d, q,
                              QuadraticFormByMatrix(
                                InvariantQuadraticForm( G ).matrix, GF(q) ) ) );

InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field size, form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsQuadraticForm ],
    function( filt, e, d, q, form )
    local g, stored, wanted, mat1, mat2, mat, matinv, gens, gg;

    # Create the default generators and form.
    g:= SpecialOrthogonalGroupCons( filt, e, d, q );
    stored:= InvariantQuadraticForm( g ).matrix;

    # If the prescribed form fits then just return.
    if stored = form!.matrix then
      return g;
    fi;

    # Compute a base change matrix.
    # (Check that the canonical forms are equal.)
    wanted:= QuadraticFormByMatrix( stored, GF(q) );
    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( wanted );
    if Forms_RESET( mat1 * form!.matrix * TransposedMat( mat1 ), d, q ) <>
       Forms_RESET( mat2 * stored * TransposedMat( mat2 ), d, q ) then
      Error( "canonical forms of <form> and <wanted> differ" );
    fi;
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    # Create the group w.r.t. the prescribed form.
    gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
    gg:= GroupWithGenerators( gens );

    if HasSize( g ) then
      SetSize( gg, Size( g ) );
    fi;

    if HasName( g ) then
      SetName( gg, Name( g ) );
    fi;

    SetInvariantQuadraticForm( gg, rec( matrix:= form!.matrix ) );
    if HasIsFullSubgroupGLorSLRespectingQuadraticForm( g ) then
      SetIsFullSubgroupGLorSLRespectingQuadraticForm( gg,
          IsFullSubgroupGLorSLRespectingQuadraticForm( g ) );
    fi;

    SetInvariantBilinearForm( gg, rec( matrix:= matinv * InvariantBilinearForm( g ).matrix * TransposedMat( matinv ) ) );
    if q mod 2 = 1 and
       HasIsFullSubgroupGLorSLRespectingBilinearForm( g ) then
      SetIsFullSubgroupGLorSLRespectingBilinearForm( gg,
          IsFullSubgroupGLorSLRespectingBilinearForm( g ) );
    fi;

    return gg;
end );


#############################################################################
##
#M  SpecialOrthogonalGroupCons( <filt>, <e>, <d>, <R>, <form> )
##
InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsField and IsFinite,
      IsMatrixOrMatrixObj ],
    { filt, e, d, F, form } -> SpecialOrthogonalGroupCons( filt, e, d,
                                 Size( F ),
                                 QuadraticFormByMatrix( form, F ) ) );

InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field, group with form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsField and IsFinite,
      IsGroup and HasInvariantQuadraticForm ],
    { filt, e, d, F, G } -> SpecialOrthogonalGroupCons( filt, e, d,
                              Size( F ),
                              QuadraticFormByMatrix(
                                InvariantQuadraticForm( G ).matrix, F ) ) );

InstallMethod( SpecialOrthogonalGroupCons,
    "matrix group for <e>, dimension, finite field, form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsField and IsFinite,
      IsQuadraticForm ],
    { filt, e, d, F, form } -> SpecialOrthogonalGroupCons( filt, e, d,
                                 Size( F ), form ) );


#############################################################################
##
#O  OmegaCons( <filt>, [<e>, <d>, <q>, ]<form> )
#O  Omega( [<filt>, ]<form> )
#O  Omega( [<filt>, ][<e>, ]<d>, <q>, <form> )
#O  Omega( [<filt>, ][<e>, ]<d>, <R>, <form> )
##
##  'Omega' is an operation hat is defined in the GAP library.
##  Thus we have to declare the variants involving a quadratic form,
##  and install the corresponding 'Omega' methods that call 'OmegaCons'.
##
##  Install the methods involving <form>, which may be either a matrix or
##  a quadratic form or a group with stored 'InvariantQuadraticForm'.
##  (The other methods have been installed in the GAP library.)
##
Perform(
    [ IsMatrixOrMatrixObj, IsQuadraticForm, IsGroup and HasInvariantQuadraticForm ],
    function( obj )
      DeclareConstructor( "OmegaCons", [ IsGroup, obj ] );
      DeclareConstructor( "OmegaCons", [ IsGroup, IsInt, IsPosInt, IsPosInt, obj ] );

      DeclareOperation( "Omega", [ obj ] );
      InstallMethod( Omega, [ obj ], x -> OmegaCons( IsMatrixGroup, x ) );

      DeclareOperation( "Omega", [ IsFunction, obj ] );
      InstallMethod( Omega, [ IsFunction, obj ], OmegaCons );

      DeclareOperation( "Omega", [ IsPosInt, IsPosInt, obj ] );
      InstallMethod( Omega, [ IsPosInt, IsPosInt, obj ],
        { d, q, x } -> OmegaCons( IsMatrixGroup, 0, d, q, x ) );

      DeclareOperation( "Omega", [ IsPosInt, IsRing, obj ] );
      InstallMethod( Omega, [ IsPosInt, IsField and IsFinite, obj ],
        { d, R, x } -> OmegaCons( IsMatrixGroup, 0, d, Size( R ), x ) );

      DeclareOperation( "Omega", [ IsFunction, IsPosInt, IsPosInt, obj ] );
      InstallMethod( Omega, [ IsFunction, IsPosInt, IsPosInt, obj ],
        { filt, d, q, x } -> OmegaCons( filt, 0, d, q, x ) );

      DeclareOperation( "Omega", [ IsFunction, IsPosInt, IsRing, obj ] );
      InstallMethod( Omega, [ IsFunction, IsPosInt, IsField and IsFinite, obj ],
        { filt, d, R, x } -> OmegaCons( filt, 0, d, Size( R ), x ) );

      DeclareOperation( "Omega", [ IsInt, IsPosInt, IsPosInt, obj ] );
      InstallMethod( Omega, [ IsInt, IsPosInt, IsPosInt, obj ],
        { e, d, q, x } -> OmegaCons( IsMatrixGroup, e, d, q, x ) );

      DeclareOperation( "Omega", [ IsInt, IsPosInt, IsRing, obj ] );
      InstallMethod( Omega, [ IsInt, IsPosInt, IsField and IsFinite, obj ],
        { e, d, R, x } -> OmegaCons( IsMatrixGroup, e, d, Size( R ), x ) );

      DeclareOperation( "Omega", [ IsFunction, IsInt, IsPosInt, IsPosInt, obj ] );
      InstallMethod( Omega, [ IsFunction, IsInt, IsPosInt, IsPosInt, obj ],
        OmegaCons );

      DeclareOperation( "Omega", [ IsFunction, IsInt, IsPosInt, IsRing, obj ] );
      InstallMethod( Omega, [ IsFunction, IsInt, IsPosInt, IsField and IsFinite, obj ],
        { filt, e, d, R, x } -> OmegaCons( filt, e, d, Size( R ), x ) );
    end );


#############################################################################
##
#M  OmegaCons( <filt>, <form> )
##
InstallMethod( OmegaCons,
    "matrix group for matrix of form",
    [ IsMatrixGroup and IsFinite, IsMatrixOrMatrixObj ],
    { filt, mat } -> OmegaCons( filt,
                       QuadraticFormByMatrix( mat, BaseDomain( mat ) ) ) );

InstallMethod( OmegaCons,
    "matrix group for group with form",
    [ IsMatrixGroup and IsFinite, IsGroup and HasInvariantQuadraticForm ],
    { filt, G } -> OmegaCons( filt,
                     QuadraticFormByMatrix(
                       InvariantQuadraticForm( G ).matrix,
                       FieldOfMatrixGroup( G ) ) ) );

InstallMethod( OmegaCons,
    "matrix group for form",
    [ IsMatrixGroup and IsFinite, IsQuadraticForm ],
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
    return OmegaCons( filt, e, d, q, form );
end );


#############################################################################
##
#M  OmegaCons( <filt>, <e>, <d>, <q>, <form> )
##
InstallMethod( OmegaCons,
    "matrix group for <e>, dimension, finite field size, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsMatrixOrMatrixObj ],
    { filt, e, d, q, mat } -> OmegaCons( filt, e, d, q,
                                QuadraticFormByMatrix( mat, GF(q) ) ) );

InstallMethod( OmegaCons,
    "matrix group for <e>, dimension, finite field size, group with form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsGroup and HasInvariantQuadraticForm ],
    { filt, e, d, q, G } -> OmegaCons( filt, e, d, q,
                              QuadraticFormByMatrix(
                                InvariantQuadraticForm( G ).matrix, GF(q) ) ) );

InstallMethod( OmegaCons,
    "matrix group for <e>, dimension, finite field size, form",
    [ IsMatrixGroup and IsFinite,
      IsInt,
      IsPosInt,
      IsPosInt,
      IsQuadraticForm ],
    function( filt, e, d, q, form )
    local g, stored, wanted, mat1, mat2, mat, matinv, gens, gg;

    # Create the default generators and form.
    g:= OmegaCons( filt, e, d, q );
    stored:= InvariantQuadraticForm( g ).matrix;

    # If the prescribed form fits then just return.
    if stored = form!.matrix then
      return g;
    fi;

    # Compute a base change matrix.
    # (Check that the canonical forms are equal.)
    wanted:= QuadraticFormByMatrix( stored, GF(q) );
    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( wanted );
    if Forms_RESET( mat1 * form!.matrix * TransposedMat( mat1 ), d, q ) <>
       Forms_RESET( mat2 * stored * TransposedMat( mat2 ), d, q ) then
      Error( "canonical forms of <form> and <wanted> differ" );
    fi;
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    # Create the group w.r.t. the prescribed form.
    gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
    gg:= GroupWithGenerators( gens );

    if HasSize( g ) then
      SetSize( gg, Size( g ) );
    fi;

    if HasName( g ) then
      SetName( gg, Name( g ) );
    fi;

    SetInvariantQuadraticForm( gg, rec( matrix:= form!.matrix ) );
    if HasIsFullSubgroupGLorSLRespectingQuadraticForm( g ) then
      SetIsFullSubgroupGLorSLRespectingQuadraticForm( gg,
          IsFullSubgroupGLorSLRespectingQuadraticForm( g ) );
    fi;

    SetInvariantBilinearForm( gg, rec( matrix:= matinv * InvariantBilinearForm( g ).matrix * TransposedMat( matinv ) ) );
    if q mod 2 = 1 and
       HasIsFullSubgroupGLorSLRespectingBilinearForm( g ) then
      SetIsFullSubgroupGLorSLRespectingBilinearForm( gg,
          IsFullSubgroupGLorSLRespectingBilinearForm( g ) );
    fi;

    return gg;
end );


#############################################################################
##
#O  GeneralUnitaryGroupCons( <filter>, <form> )
#O  GeneralUnitaryGroupCons( <filter>, <d>, <q>, <form> )
##
##  'GeneralUnitaryGroup' is a plain function that is defined in the GAP
##  library.
##  It calls 'GeneralUnitaryGroupCons',
##  thus we have to declare the variants involving a quadratic form,
##  and install the corresponding methods.
##
Perform(
    [ IsMatrixOrMatrixObj, IsHermitianForm, IsGroup and HasInvariantSesquilinearForm ],
    function( obj )
      DeclareConstructor( "GeneralUnitaryGroupCons", [ IsGroup, obj ] );
      DeclareConstructor( "GeneralUnitaryGroupCons",
        [ IsGroup, IsPosInt, IsPosInt, obj ] );
    end );


#############################################################################
##
#M  GeneralUnitaryGroupCons( <filt>, <form> )
##
InstallMethod( GeneralUnitaryGroupCons,
    "matrix group for matrix of form",
    [ IsMatrixGroup and IsFinite, IsMatrixOrMatrixObj ],
    function( filt, mat )
    local q, form;

    # Guess the intended field.
    q:= Size( BaseDomain( mat ) );
    if q = RootInt( q, 2 )^2 then
      form:= HermitianFormByMatrix( mat, GF(q) );
    else
      form:= HermitianFormByMatrix( mat, GF(q^2) );
    fi;
    return GeneralUnitaryGroupCons( filt, form );
    end );

InstallMethod( GeneralUnitaryGroupCons,
    "matrix group for group with form",
    [ IsMatrixGroup and IsFinite, IsGroup and HasInvariantSesquilinearForm ],
    { filt, G } -> GeneralUnitaryGroupCons( filt,
                     HermitianFormByMatrix(
                       InvariantSesquilinearForm( G ).matrix,
                       FieldOfMatrixGroup( G ) ) ) );

InstallMethod( GeneralUnitaryGroupCons,
    "matrix group for form",
    [ IsMatrixGroup and IsFinite, IsHermitianForm ],
    { filt, form } -> GeneralUnitaryGroupCons( filt,
                        NumberRows( form!.matrix ),
                        RootInt( Size( form!.basefield ), 2 ), form ) );


#############################################################################
##
#M  GeneralUnitaryGroupCons( <filt>, <d>, <q>, <form> )
##
InstallMethod( GeneralUnitaryGroupCons,
    "matrix group for dimension, finite field size, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsMatrixOrMatrixObj ],
    { filt, d, q, mat } -> GeneralUnitaryGroupCons( filt, d, q,
                             HermitianFormByMatrix( mat, GF(q^2) ) ) );

InstallMethod( GeneralUnitaryGroupCons,
    "matrix group for dimension, finite field size, group with form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsGroup and HasInvariantSesquilinearForm ],
    { filt, d, q, G } -> GeneralUnitaryGroupCons( filt, d, q,
                           HermitianFormByMatrix(
                             InvariantSesquilinearForm( G ).matrix, GF(q^2) ) ) );

InstallMethod( GeneralUnitaryGroupCons,
    "matrix group for dimension, finite field size, form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsHermitianForm ],
    function( filt, d, q, form )
    local g, stored, wanted, mat1, mat2, mat, matinv, gens, gg;

    # Create the default generators and form.
    g:= GeneralUnitaryGroupCons( filt, d, q );
    stored:= InvariantSesquilinearForm( g ).matrix;

    # If the prescribed form fits then just return.
    if stored = form!.matrix then
      return g;
    fi;

    # Compute a base change matrix.
    # (Check that the canonical forms are equal.)
    wanted:= HermitianFormByMatrix( stored, GF(q^2) );
    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( wanted );
    if mat1 * form!.matrix * TransposedFrobeniusMat( mat1, q ) <>
       mat2 * stored * TransposedFrobeniusMat( mat2, q ) then
      Error( "canonical forms of <form> and <wanted> differ" );
    fi;
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    # Create the group w.r.t. the prescribed form.
    gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
    gg:= GroupWithGenerators( gens );

    if HasSize( g ) then
      SetSize( gg, Size( g ) );
    fi;

    if HasName( g ) then
      SetName( gg, Name( g ) );
    fi;

    SetInvariantSesquilinearForm( gg, rec( matrix:= form!.matrix ) );
    if HasIsFullSubgroupGLorSLRespectingSesquilinearForm( g ) then
      SetIsFullSubgroupGLorSLRespectingSesquilinearForm( gg,
          IsFullSubgroupGLorSLRespectingSesquilinearForm( g ) );
    fi;

    return gg;
end );


#############################################################################
##
#O  SpecialUnitaryGroupCons( <filter>, <form> )
#O  SpecialUnitaryGroupCons( <filter>, <d>, <q>, <form> )
##
##  'SpecialUnitaryGroup' is a plain function that is defined in the GAP
##  library.
##  It calls 'SpecialUnitaryGroupCons',
##  thus we have to declare the variants involving a quadratic form,
##  and install the corresponding methods.
##
Perform(
    [ IsMatrixOrMatrixObj, IsHermitianForm, IsGroup and HasInvariantSesquilinearForm ],
    function( obj )
      DeclareConstructor( "SpecialUnitaryGroupCons", [ IsGroup, obj ] );
      DeclareConstructor( "SpecialUnitaryGroupCons",
        [ IsGroup, IsPosInt, IsPosInt, obj ] );
    end );


#############################################################################
##
#M  SpecialUnitaryGroupCons( <filt>, <form> )
##
InstallMethod( SpecialUnitaryGroupCons,
    "matrix group for matrix of form",
    [ IsMatrixGroup and IsFinite, IsMatrixOrMatrixObj ],
    function( filt, mat )
    local q, form;

    # Guess the intended field.
    q:= Size( BaseDomain( mat ) );
    if q = RootInt( q, 2 )^2 then
      form:= HermitianFormByMatrix( mat, GF(q) );
    else
      form:= HermitianFormByMatrix( mat, GF(q^2) );
    fi;
    return SpecialUnitaryGroupCons( filt, form );
    end );

InstallMethod( SpecialUnitaryGroupCons,
    "matrix group for group with form",
    [ IsMatrixGroup and IsFinite, IsGroup and HasInvariantSesquilinearForm ],
    { filt, G } -> SpecialUnitaryGroupCons( filt,
                     HermitianFormByMatrix(
                       InvariantSesquilinearForm( G ).matrix,
                       FieldOfMatrixGroup( G ) ) ) );

InstallMethod( SpecialUnitaryGroupCons,
    "matrix group for form",
    [ IsMatrixGroup and IsFinite, IsHermitianForm ],
    { filt, form } -> SpecialUnitaryGroupCons( filt,
                        NumberRows( form!.matrix ),
                        RootInt( Size( form!.basefield ), 2 ), form ) );


#############################################################################
##
#M  SpecialUnitaryGroupCons( <filt>, <d>, <q>, <form> )
##
InstallMethod( SpecialUnitaryGroupCons,
    "matrix group for dimension, finite field size, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsMatrixOrMatrixObj ],
    { filt, d, q, mat } -> SpecialUnitaryGroupCons( filt, d, q,
                             HermitianFormByMatrix( mat, GF(q^2) ) ) );

InstallMethod( SpecialUnitaryGroupCons,
    "matrix group for dimension, finite field size, group with form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsGroup and HasInvariantSesquilinearForm ],
    { filt, d, q, G } -> SpecialUnitaryGroupCons( filt, d, q,
                           HermitianFormByMatrix(
                             InvariantSesquilinearForm( G ).matrix, GF(q^2) ) ) );

InstallMethod( SpecialUnitaryGroupCons,
    "matrix group for dimension, finite field size, form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsHermitianForm ],
    function( filt, d, q, form )
    local g, stored, wanted, mat1, mat2, mat, matinv, gens, gg;

    # Create the default generators and form.
    g:= SpecialUnitaryGroupCons( filt, d, q );
    stored:= InvariantSesquilinearForm( g ).matrix;

    # If the prescribed form fits then just return.
    if stored = form!.matrix then
      return g;
    fi;

    # Compute a base change matrix.
    # (Check that the canonical forms are equal.)
    wanted:= HermitianFormByMatrix( stored, GF(q^2) );
    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( wanted );
    if mat1 * form!.matrix * TransposedFrobeniusMat( mat1, q ) <>
       mat2 * stored * TransposedFrobeniusMat( mat2, q ) then
      Error( "canonical forms of <form> and <wanted> differ" );
    fi;
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    # Create the group w.r.t. the prescribed form.
    gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
    gg:= GroupWithGenerators( gens );

    if HasSize( g ) then
      SetSize( gg, Size( g ) );
    fi;

    if HasName( g ) then
      SetName( gg, Name( g ) );
    fi;

    SetInvariantSesquilinearForm( gg, rec( matrix:= form!.matrix ) );
    if HasIsFullSubgroupGLorSLRespectingSesquilinearForm( g ) then
      SetIsFullSubgroupGLorSLRespectingSesquilinearForm( gg,
          IsFullSubgroupGLorSLRespectingSesquilinearForm( g ) );
    fi;

    return gg;
end );


#############################################################################
##
#O  SymplecticGroupCons( <filter>, <form> )
#O  SymplecticGroupCons( <filter>, <d>, <q>, <form> )
#O  SymplecticGroupCons( <filter>, <d>, <R>, <form> )
##
##  'SymplecticGroup' is a plain function that is defined in the GAP
##  library.
##  It calls 'SymplecticGroupCons',
##  thus we have to declare the variants involving a quadratic form,
##  and install the corresponding methods.
##
Perform(
    [ IsMatrixOrMatrixObj, IsBilinearForm, IsGroup and HasInvariantBilinearForm ],
    function( obj )
      DeclareConstructor( "SymplecticGroupCons", [ IsGroup, obj ] );
      DeclareConstructor( "SymplecticGroupCons",
        [ IsGroup, IsPosInt, IsPosInt, obj ] );
      DeclareConstructor( "SymplecticGroupCons",
        [ IsGroup, IsPosInt, IsRing, obj ] );
    end );


#############################################################################
##
#M  SymplecticGroupCons( <filt>, <form> )
##
InstallMethod( SymplecticGroupCons,
    "matrix group for matrix of form",
    [ IsMatrixGroup and IsFinite, IsMatrixOrMatrixObj ],
    { filt, mat } -> SymplecticGroupCons( filt,
                       BilinearFormByMatrix( mat, BaseDomain( mat ) ) ) );

InstallMethod( SymplecticGroupCons,
    "matrix group for group with form",
    [ IsMatrixGroup and IsFinite, IsGroup and HasInvariantBilinearForm ],
    { filt, G } -> SymplecticGroupCons( filt,
                     BilinearFormByMatrix(
                       InvariantBilinearForm( G ).matrix,
                       FieldOfMatrixGroup( G ) ) ) );

InstallMethod( SymplecticGroupCons,
    "matrix group for form",
    [ IsMatrixGroup and IsFinite, IsBilinearForm ],
    { filt, form } -> SymplecticGroupCons( filt,
                        NumberRows( form!.matrix ),
                        Size( form!.basefield ), form ) );


#############################################################################
##
#M  SymplecticGroupCons( <filt>, <d>, <q>, <form> )
##
InstallMethod( SymplecticGroupCons,
    "matrix group for dimension, finite field size, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsMatrixOrMatrixObj ],
    { filt, d, q, mat } -> SymplecticGroupCons( filt, d, q,
                             BilinearFormByMatrix( mat, GF(q) ) ) );

InstallMethod( SymplecticGroupCons,
    "matrix group for dimension, finite field size, group with form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsGroup and HasInvariantBilinearForm ],
    { filt, d, q, G } -> SymplecticGroupCons( filt, d, q,
                           BilinearFormByMatrix(
                             InvariantBilinearForm( G ).matrix, GF(q) ) ) );

InstallMethod( SymplecticGroupCons,
    "matrix group for dimension, finite field size, form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsPosInt,
      IsBilinearForm ],
    function( filt, d, q, form )
    local g, stored, wanted, mat1, mat2, mat, matinv, gens, gg;

    # Create the default generators and form.
    g:= SymplecticGroupCons( filt, d, q );
    stored:= InvariantBilinearForm( g ).matrix;

    # If the prescribed form fits then just return.
    if stored = form!.matrix then
      return g;
    fi;

    # Compute a base change matrix.
    # (Check that the canonical forms are equal.)
    wanted:= BilinearFormByMatrix( stored, GF(q) );
    mat1:= BaseChangeToCanonical( form );
    mat2:= BaseChangeToCanonical( wanted );
    if mat1 * form!.matrix * TransposedMat( mat1 ) <>
       mat2 * stored * TransposedMat( mat2 ) then
      Error( "canonical forms of <form> and <wanted> differ" );
    fi;
    mat:= mat2^-1 * mat1;
    matinv:= mat^-1;

    # Create the group w.r.t. the prescribed form.
    gens:= List( GeneratorsOfGroup( g ), x -> matinv * x * mat );
    gg:= GroupWithGenerators( gens );

    if HasSize( g ) then
      SetSize( gg, Size( g ) );
    fi;

    if HasName( g ) then
      SetName( gg, Name( g ) );
    fi;

    SetInvariantBilinearForm( gg, rec( matrix:= form!.matrix ) );
    if HasIsFullSubgroupGLorSLRespectingBilinearForm( g ) then
      SetIsFullSubgroupGLorSLRespectingBilinearForm( gg,
          IsFullSubgroupGLorSLRespectingBilinearForm( g ) );
    fi;

    return gg;
end );


#############################################################################
##
#M  SymplecticGroupCons( <filt>, <d>, <R>, <form> )
##
InstallMethod( SymplecticGroupCons,
    "matrix group for dimension, finite field, matrix of form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsField and IsFinite,
      IsMatrixOrMatrixObj ],
    { filt, d, F, form } -> SymplecticGroupCons( filt, d, Size( F ),
                              BilinearFormByMatrix( form, F ) ) );

InstallMethod( SymplecticGroupCons,
    "matrix group for dimension, finite field, group with form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsField and IsFinite,
      IsGroup and HasInvariantBilinearForm ],
    { filt, d, F, G } -> SymplecticGroupCons( filt, d, Size( F ),
                           BilinearFormByMatrix(
                             InvariantBilinearForm( G ).matrix, F ) ) );

InstallMethod( SymplecticGroupCons,
    "matrix group for dimension, finite field, form",
    [ IsMatrixGroup and IsFinite,
      IsPosInt,
      IsField and IsFinite,
      IsBilinearForm ],
    { filt, d, F, form } -> SymplecticGroupCons( filt, d, Size( F ), form ) );
