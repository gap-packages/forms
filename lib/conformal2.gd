#############################################################################
##
##  conformal2.gd         'Forms' package
##
##  Declare the variants of 'ConformalOrthogonalGroupCons' and
##  'ConformalUnitaryGroupCons' involving a form as an argument.
##  GAP provides these constructors from version 4.17 on.
##

#############################################################################
##
#O  ConformalOrthogonalGroupCons( <filter>, <form> )
#O  ConformalOrthogonalGroupCons( <filter>, <e>, <d>, <q>, <form> )
#O  ConformalOrthogonalGroupCons( <filter>, <e>, <d>, <R>, <form> )
##
Perform(
    [ IsMatrixOrMatrixObj, IsQuadraticForm, IsGroup ],
    function( obj )
      DeclareConstructor( "ConformalOrthogonalGroupCons", [ IsGroup, obj ] );
      DeclareConstructor( "ConformalOrthogonalGroupCons",
        [ IsGroup, IsInt, IsPosInt, IsPosInt, obj ] );
      DeclareConstructor( "ConformalOrthogonalGroupCons",
        [ IsGroup, IsInt, IsPosInt, IsRing, obj ] );
    end );


#############################################################################
##
#O  ConformalUnitaryGroupCons( <filter>, <form> )
#O  ConformalUnitaryGroupCons( <filter>, <d>, <q>, <form> )
##
Perform(
    [ IsMatrixOrMatrixObj, IsHermitianForm, IsGroup ],
    function( obj )
      DeclareConstructor( "ConformalUnitaryGroupCons", [ IsGroup, obj ] );
      DeclareConstructor( "ConformalUnitaryGroupCons",
        [ IsGroup, IsPosInt, IsPosInt, obj ] );
    end );
