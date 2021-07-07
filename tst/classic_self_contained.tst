#@local is_equal, q, F, d, es, e, g, stored, pi, permmat, form, gg, F2

gap> START_TEST( "Forms: classic_self_contained.tst" );

# Test the methods for constructing classical groups w.r.t. prescribed forms,
# but without assuming that the global functions in the GAP library
# already delegate to these methods.

# Provide an auxiliary function (until GAP's '=' gets fast).
gap> is_equal:= function( G1, G2 )
>      return IsSubset( G1, GeneratorsOfGroup( G2 ) ) and
>             IsSubset( G2, GeneratorsOfGroup( G1 ) );
>    end;;

# Test the creation of orthogonal groups.
gap> for q in [ 2, 3, 4, 5, 7, 8 ] do
>      F:= GF(q);
>      for d in [ 3 .. 5 ] do
>        if IsEvenInt( d ) then
>          es:= [ -1, 1 ];
>        else
>          es:= [ 0 ];
>        fi;
>        for e in es do
>          # GO(e,d,q)
>          g:= GeneralOrthogonalGroup( e, d, q );
>          stored:= InvariantQuadraticForm( g ).matrix;
>          pi:= PermutationMat( (1,2,3), d, F );
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= QuadraticFormByMatrix( stored, F );
>          gg:= GeneralOrthogonalGroupCons( IsMatrixGroup, e, d, q, permmat );
>          if not ( is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, g ) ) and
>                   ( is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, form ) ) and
>                   is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, e, d, q, g ) ) and
>                   is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, e, d, q, stored ) ) and
>                   is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, e, d, q, form ) ) and
>                   is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, e, d, F, g ) ) and
>                   is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, e, d, F, stored ) ) and
>                   is_equal( g, GeneralOrthogonalGroupCons( IsMatrixGroup, e, d, F, form ) ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with GeneralOrthogonalGroupCons( IsMatrixGroup, ", e, ",", d, ",", q, ")" );
>          fi;
>          # SO(e,d,q)
>          g:= SpecialOrthogonalGroup( e, d, q );
>          stored:= InvariantQuadraticForm( g ).matrix;
>          pi:= PermutationMat( (1,2,3), d, F );
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= QuadraticFormByMatrix( stored, F );
>          gg:= SpecialOrthogonalGroupCons( IsMatrixGroup, e, d, q, permmat );
>          if not ( is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, g ) ) and
>                   ( is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, form ) ) and
>                   is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, e, d, q, g ) ) and
>                   is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, e, d, q, stored ) ) and
>                   is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, e, d, q, form ) ) and
>                   is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, e, d, F, g ) ) and
>                   is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, e, d, F, stored ) ) and
>                   is_equal( g, SpecialOrthogonalGroupCons( IsMatrixGroup, e, d, F, form ) ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with SpecialOrthogonalGroupCons( IsMatrixGroup, ", e, ",", d, ",", q, ")" );
>          fi;
>          # Omega(e,d,q)
>          g:= Omega( e, d, q );
>          stored:= InvariantQuadraticForm( g ).matrix;
>          pi:= PermutationMat( (1,2,3), d, F );
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= QuadraticFormByMatrix( stored, F );
>          gg:= Omega( e, d, q, permmat );
>          if not ( is_equal( g, Omega( g ) ) and
>                   ( is_equal( g, Omega( stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, Omega( form ) ) and
>                   is_equal( g, Omega( e, d, q, g ) ) and
>                   is_equal( g, Omega( e, d, q, stored ) ) and
>                   is_equal( g, Omega( e, d, q, form ) ) and
>                   is_equal( g, Omega( e, d, F, g ) ) and
>                   is_equal( g, Omega( e, d, F, stored ) ) and
>                   is_equal( g, Omega( e, d, F, form ) ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with Omega(", e, ",", d, ",", q, ")" );
>          fi;
>          if e = 0 then
>            if not ( is_equal( g, Omega( d, q, g ) ) and
>                     is_equal( g, Omega( d, q, stored ) ) and
>                     is_equal( g, Omega( d, q, form ) ) and
>                     is_equal( g, Omega( d, F, g ) ) and
>                     is_equal( g, Omega( d, F, stored ) ) and
>                     is_equal( g, Omega( d, F, form ) ) ) then
>              Error( "problem with Omega", d, ",", q, ")" );
>            fi;
>          fi;
>        od;
>      od;
>    od;

# Test the creation of unitary groups.
gap> for q in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25 ] do
>      F:= GF(q);
>      F2:= GF(q^2);
>      for d in [ 2 .. 8 ] do
>        # GU(d,q)
>        g:= GeneralUnitaryGroup( d, q );
>        stored:= InvariantSesquilinearForm( g ).matrix;
>        pi:= PermutationMat( (1,2), d, F );
>        permmat:= pi * stored * TransposedMat( pi );
>        form:= HermitianFormByMatrix( stored, F2 );
>        gg:= GeneralUnitaryGroupCons( IsMatrixGroup, d, q, permmat );
>        if not ( is_equal( g, GeneralUnitaryGroupCons( IsMatrixGroup, g ) ) and
>                 ( is_equal( g, GeneralUnitaryGroupCons( IsMatrixGroup, stored ) ) or
>                   BaseDomain( stored ) <> F ) and
>                 is_equal( g, GeneralUnitaryGroupCons( IsMatrixGroup, form ) ) and
>                 is_equal( g, GeneralUnitaryGroupCons( IsMatrixGroup, d, q, g ) ) and
>                 is_equal( g, GeneralUnitaryGroupCons( IsMatrixGroup, d, q, stored ) ) and
>                 is_equal( g, GeneralUnitaryGroupCons( IsMatrixGroup, d, q, form ) ) and
>                 IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                 IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>          Error( "problem with GeneralUnitaryGroupCons( IsMatrixGroup, ", d, ",", q, ")" );
>        fi;
>        # SU(d,q)
>        g:= SpecialUnitaryGroup( d, q );
>        stored:= InvariantSesquilinearForm( g ).matrix;
>        pi:= PermutationMat( (1,2), d, F );
>        permmat:= pi * stored * TransposedMat( pi );
>        form:= HermitianFormByMatrix( stored, F2 );
>        gg:= SpecialUnitaryGroupCons( IsMatrixGroup, d, q, permmat );
>        if not ( is_equal( g, SpecialUnitaryGroupCons( IsMatrixGroup, g ) ) and
>                 ( is_equal( g, SpecialUnitaryGroupCons( IsMatrixGroup, stored ) ) or
>                   BaseDomain( stored ) <> F ) and
>                 is_equal( g, SpecialUnitaryGroupCons( IsMatrixGroup, form ) ) and
>                 is_equal( g, SpecialUnitaryGroupCons( IsMatrixGroup, d, q, g ) ) and
>                 is_equal( g, SpecialUnitaryGroupCons( IsMatrixGroup, d, q, stored ) ) and
>                 is_equal( g, SpecialUnitaryGroupCons( IsMatrixGroup, d, q, form ) ) and
>                 IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                 IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>          Error( "problem with SpecialUnitaryGroupCons( IsMatrixGroup, ", d, ",", q, ")" );
>        fi;
>      od;
>    od;

# Test the creation of symplectic groups.
gap> for q in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25 ] do
>      F:= GF(q);
>      for d in [ 2, 4 .. 8 ] do
>        g:= SymplecticGroup( d, q );
>        stored:= InvariantBilinearForm( g ).matrix;
>        pi:= PermutationMat( (1,2), d, F );
>        permmat:= pi * stored * TransposedMat( pi );
>        form:= BilinearFormByMatrix( stored, F );
>        gg:= SymplecticGroupCons( IsMatrixGroup, d, q, permmat );
>        if not ( is_equal( g, SymplecticGroupCons( IsMatrixGroup, g ) ) and
>                 ( is_equal( g, SymplecticGroupCons(IsMatrixGroup,  stored ) ) or
>                   BaseDomain( stored ) <> F ) and
>                 is_equal( g, SymplecticGroupCons( IsMatrixGroup, form ) ) and
>                 is_equal( g, SymplecticGroupCons( IsMatrixGroup, d, q, g ) ) and
>                 is_equal( g, SymplecticGroupCons( IsMatrixGroup, d, q, stored ) ) and
>                 is_equal( g, SymplecticGroupCons( IsMatrixGroup, d, q, form ) ) and
>                 is_equal( g, SymplecticGroupCons( IsMatrixGroup, d, F, g ) ) and
>                 is_equal( g, SymplecticGroupCons( IsMatrixGroup, d, F, stored ) ) and
>                 is_equal( g, SymplecticGroupCons( IsMatrixGroup, d, F, form ) ) and
>                 IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                 IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>          Error( "problem with SymplecticGroupCons( IsMatrixGroup, ", d, ",", q, ")" );
>        fi;
>      od;
>    od;

##
gap> STOP_TEST( "classic_self_contained.tst" );
