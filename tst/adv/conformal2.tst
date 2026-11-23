#@local is_equal, q, F, F2, d, es, e, g, stored, pi, permmat, form, gg, go
#@local gu

gap> START_TEST( "Forms: conformal2.tst" );

# Provide an auxiliary function (until GAP's '=' gets fast).
gap> is_equal:= function( G1, G2 )
>      return IsSubset( G1, GeneratorsOfGroup( G2 ) ) and
>             IsSubset( G2, GeneratorsOfGroup( G1 ) );
>    end;;

# Test the creation of conformal orthogonal groups.
gap> for q in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25 ] do
>      F:= GF(q);
>      for d in [ 2 .. 6 ] do
>        if IsEvenInt( d ) then
>          es:= [ -1, 1 ];
>        else
>          es:= [ 0 ];
>        fi;
>        if d = 2 then
>          pi:= PermutationMat( (1,2), d, F );
>        else
>          pi:= PermutationMat( (1,2,3), d, F );
>        fi;
>        for e in es do
>          g:= ConformalOrthogonalGroup( e, d, q );
>          stored:= InvariantQuadraticFormUpToScalars( g ).matrix;
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= QuadraticFormByMatrix( stored, F );
>          gg:= ConformalOrthogonalGroup( e, d, q, permmat );
>          if not ( is_equal( g, ConformalOrthogonalGroup( g ) ) and
>                   ( is_equal( g, ConformalOrthogonalGroup( stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, ConformalOrthogonalGroup( form ) ) and
>                   is_equal( g, ConformalOrthogonalGroup( e, d, q, g ) ) and
>                   is_equal( g, ConformalOrthogonalGroup( e, d, q, stored ) ) and
>                   is_equal( g, ConformalOrthogonalGroup( e, d, q, form ) ) and
>                   is_equal( g, ConformalOrthogonalGroup( e, d, F, g ) ) and
>                   is_equal( g, ConformalOrthogonalGroup( e, d, F, stored ) ) and
>                   is_equal( g, ConformalOrthogonalGroup( e, d, F, form ) ) and
>                   is_equal( g, ConformalOrthogonalGroup( e, d, q,
>                                    PrimitiveRoot( F ) * stored ) ) and
>                   Size( gg ) = Size( g ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with CO(", e, ",", d, ",", q, ")" );
>          fi;
>          if e = 0 then
>            if not ( is_equal( g, ConformalOrthogonalGroup( d, q, g ) ) and
>                     is_equal( g, ConformalOrthogonalGroup( d, q, stored ) ) and
>                     is_equal( g, ConformalOrthogonalGroup( d, q, form ) ) and
>                     is_equal( g, ConformalOrthogonalGroup( d, F, g ) ) and
>                     is_equal( g, ConformalOrthogonalGroup( d, F, stored ) ) and
>                     is_equal( g, ConformalOrthogonalGroup( d, F, form ) ) ) then
>              Error( "problem with CO(", d, ",", q, ")" );
>            fi;
>          fi;
>          go:= GeneralOrthogonalGroup( e, d, q, permmat );
>          if not IsSubset( ConformalOrthogonalGroup( go ), go ) then
>            Error( "problem with CO(", e, ",", d, ",", q, ")" );
>          fi;
>        od;
>      od;
>    od;

# Test the creation of conformal unitary groups.
gap> for q in [ 2, 3, 4, 5, 7, 8, 9 ] do
>      F:= GF(q);
>      F2:= GF(q^2);
>      for d in [ 2 .. 6 ] do
>        g:= ConformalUnitaryGroup( d, q );
>        stored:= InvariantSesquilinearFormUpToScalars( g ).matrix;
>        pi:= PermutationMat( (1,2), d, F );
>        permmat:= pi * stored * TransposedMat( pi );
>        form:= HermitianFormByMatrix( stored, F2 );
>        gg:= ConformalUnitaryGroup( d, q, permmat );
>        if not ( is_equal( g, ConformalUnitaryGroup( g ) ) and
>                 ( is_equal( g, ConformalUnitaryGroup( stored ) ) or
>                   BaseDomain( stored ) <> F or IsSquareInt(q) ) and
>                 is_equal( g, ConformalUnitaryGroup( form ) ) and
>                 is_equal( g, ConformalUnitaryGroup( d, q, g ) ) and
>                 is_equal( g, ConformalUnitaryGroup( d, q, stored ) ) and
>                 is_equal( g, ConformalUnitaryGroup( d, q, form ) ) and
>                 Size( gg ) = Size( g ) and
>                 IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                 IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>          Error( "problem with CU(", d, ",", q, ")" );
>        fi;
>        gu:= GeneralUnitaryGroup( d, q, permmat );
>        if not IsSubset( ConformalUnitaryGroup( gu ), gu ) then
>          Error( "problem with CU(", d, ",", q, ")" );
>        fi;
>      od;
>    od;

##
gap> STOP_TEST( "Forms: conformal2.tst" );
