#############################################################################
##
##  conformal.gd          'Forms' package
##

#############################################################################
##
#O  ConformalSymplecticGroupCons( <filter>, <form> )
#O  ConformalSymplecticGroupCons( <filter>, <matrix> )
#O  ConformalSymplecticGroupCons( <filter>, <G> )
#O  ConformalSymplecticGroupCons( <filter>, <d>, <R>, <form> )
#O  ConformalSymplecticGroupCons( <filter>, <d>, <R>, <matrix> )
#O  ConformalSymplecticGroupCons( <filter>, <d>, <R>, <G> )
#O  ConformalSymplecticGroupCons( <filter>, <d>, <q>, <form> )
#O  ConformalSymplecticGroupCons( <filter>, <d>, <q>, <matrix> )
#O  ConformalSymplecticGroupCons( <filter>, <d>, <q>, <G> )
##
##  Declare the variants involving a bilinear form as an argument.
##
Perform(
    [ IsMatrixOrMatrixObj, IsBilinearForm, IsGroup ],
    function( obj )
      DeclareConstructor( "ConformalSymplecticGroupCons", [ IsGroup, obj ] );
      DeclareConstructor( "ConformalSymplecticGroupCons",
        [ IsGroup, IsPosInt, IsRing, obj ] );
      DeclareConstructor( "ConformalSymplecticGroupCons",
        [ IsGroup, IsPosInt, IsPosInt, obj ] );
    end );
