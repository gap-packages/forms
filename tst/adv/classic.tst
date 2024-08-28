#@local is_equal, q, F, d, es, e, g, stored, pi, permmat, form, gg, F2, mat

gap> START_TEST( "Forms: classic.tst" );

# Test the methods for constructing classical groups w.r.t. prescribed forms,
# by calling the global functions in the GAP library that delegate to these
# methods.  Thus we also test these functions.

# Provide an auxiliary function (until GAP's '=' gets fast).
gap> is_equal:= function( G1, G2 )
>      return IsSubset( G1, GeneratorsOfGroup( G2 ) ) and
>             IsSubset( G2, GeneratorsOfGroup( G1 ) );
>    end;;

# Test the creation of orthogonal groups.
gap> for q in [ 2, 3, 4, 5, 8 ] do
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
>          gg:= GeneralOrthogonalGroup( e, d, q, permmat );
>          if not ( is_equal( g, GeneralOrthogonalGroup( g ) ) and
>                   ( is_equal( g, GeneralOrthogonalGroup( stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, GeneralOrthogonalGroup( form ) ) and
>                   is_equal( g, GeneralOrthogonalGroup( e, d, q, g ) ) and
>                   is_equal( g, GeneralOrthogonalGroup( e, d, q, stored ) ) and
>                   is_equal( g, GeneralOrthogonalGroup( e, d, q, form ) ) and
>                   is_equal( g, GeneralOrthogonalGroup( e, d, F, g ) ) and
>                   is_equal( g, GeneralOrthogonalGroup( e, d, F, stored ) ) and
>                   is_equal( g, GeneralOrthogonalGroup( e, d, F, form ) ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with GO(", e, ",", d, ",", q, ")" );
>          fi;
>          if e = 0 then
>            if not ( is_equal( g, GeneralOrthogonalGroup( d, q, g ) ) and
>                     is_equal( g, GeneralOrthogonalGroup( d, q, stored ) ) and
>                     is_equal( g, GeneralOrthogonalGroup( d, q, form ) ) and
>                     is_equal( g, GeneralOrthogonalGroup( d, F, g ) ) and
>                     is_equal( g, GeneralOrthogonalGroup( d, F, stored ) ) and
>                     is_equal( g, GeneralOrthogonalGroup( d, F, form ) ) ) then
>              Error( "problem with GO(", d, ",", q, ")" );
>            fi;
>          fi;
>          # SO(e,d,q)
>          g:= SpecialOrthogonalGroup( e, d, q );
>          stored:= InvariantQuadraticForm( g ).matrix;
>          pi:= PermutationMat( (1,2,3), d, F );
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= QuadraticFormByMatrix( stored, F );
>          gg:= SpecialOrthogonalGroup( e, d, q, permmat );
>          if not ( is_equal( g, SpecialOrthogonalGroup( g ) ) and
>                   ( is_equal( g, SpecialOrthogonalGroup( stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, SpecialOrthogonalGroup( form ) ) and
>                   is_equal( g, SpecialOrthogonalGroup( e, d, q, g ) ) and
>                   is_equal( g, SpecialOrthogonalGroup( e, d, q, stored ) ) and
>                   is_equal( g, SpecialOrthogonalGroup( e, d, q, form ) ) and
>                   is_equal( g, SpecialOrthogonalGroup( e, d, F, g ) ) and
>                   is_equal( g, SpecialOrthogonalGroup( e, d, F, stored ) ) and
>                   is_equal( g, SpecialOrthogonalGroup( e, d, F, form ) ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with SO(", e, ",", d, ",", q, ")" );
>          fi;
>          if e = 0 then
>            if not ( is_equal( g, SpecialOrthogonalGroup( d, q, g ) ) and
>                     is_equal( g, SpecialOrthogonalGroup( d, q, stored ) ) and
>                     is_equal( g, SpecialOrthogonalGroup( d, q, form ) ) and
>                     is_equal( g, SpecialOrthogonalGroup( d, F, g ) ) and
>                     is_equal( g, SpecialOrthogonalGroup( d, F, stored ) ) and
>                     is_equal( g, SpecialOrthogonalGroup( d, F, form ) ) ) then
>              Error( "problem with SO(", d, ",", q, ")" );
>            fi;
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
gap> for q in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25 ] do
>      F:= GF(q);
>      F2:= GF(q^2);
>      for d in [ 2 .. 8 ] do
>        # GU(d,q)
>        g:= GeneralUnitaryGroup( d, q );
>        stored:= InvariantSesquilinearForm( g ).matrix;
>        pi:= PermutationMat( (1,2), d, F );
>        permmat:= pi * stored * TransposedMat( pi );
>        form:= HermitianFormByMatrix( stored, F2 );
>        gg:= GeneralUnitaryGroup( d, q, permmat );
>        if not ( is_equal( g, GeneralUnitaryGroup( g ) ) and
>                 ( is_equal( g, GeneralUnitaryGroup( stored ) ) or
>                   BaseDomain( stored ) <> F or IsSquareInt(q) ) and
>                 is_equal( g, GeneralUnitaryGroup( form ) ) and
>                 is_equal( g, GeneralUnitaryGroup( d, q, g ) ) and
>                 is_equal( g, GeneralUnitaryGroup( d, q, stored ) ) and
>                 is_equal( g, GeneralUnitaryGroup( d, q, form ) ) and
>                 IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                 IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>          Error( "problem with GU(", d, ",", q, ")" );
>        fi;
>        # SU(d,q)
>        g:= SpecialUnitaryGroup( d, q );
>        stored:= InvariantSesquilinearForm( g ).matrix;
>        pi:= PermutationMat( (1,2), d, F );
>        permmat:= pi * stored * TransposedMat( pi );
>        form:= HermitianFormByMatrix( stored, F2 );
>        gg:= SpecialUnitaryGroup( d, q, permmat );
>        if not ( is_equal( g, SpecialUnitaryGroup( g ) ) and
>                 ( is_equal( g, SpecialUnitaryGroup( stored ) ) or
>                   BaseDomain( stored ) <> F or IsSquareInt(q) ) and
>                 is_equal( g, SpecialUnitaryGroup( form ) ) and
>                 is_equal( g, SpecialUnitaryGroup( d, q, g ) ) and
>                 is_equal( g, SpecialUnitaryGroup( d, q, stored ) ) and
>                 is_equal( g, SpecialUnitaryGroup( d, q, form ) ) and
>                 IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                 IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>          Error( "problem with SU(", d, ",", q, ")" );
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
>        gg:= SymplecticGroup( d, q, permmat );
>        if not ( is_equal( g, SymplecticGroup( g ) ) and
>                 ( is_equal( g, SymplecticGroup( stored ) ) or
>                   BaseDomain( stored ) <> F ) and
>                 is_equal( g, SymplecticGroup( form ) ) and
>                 is_equal( g, SymplecticGroup( d, q, g ) ) and
>                 is_equal( g, SymplecticGroup( d, q, stored ) ) and
>                 is_equal( g, SymplecticGroup( d, q, form ) ) and
>                 is_equal( g, SymplecticGroup( d, F, g ) ) and
>                 is_equal( g, SymplecticGroup( d, F, stored ) ) and
>                 is_equal( g, SymplecticGroup( d, F, form ) ) and
>                 IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                 IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>          Error( "problem with Sp(", d, ",", q, ")" );
>        fi;
>      od;
>    od;

# Test a form given by a matrix that cannot be transformed to the
# default matrix with a base change.
gap> mat:= IdentityMat( 3, GF(7) );;
gap> g:= GO( 3, 7, mat );;
gap> InvariantQuadraticForm( g ).matrix = mat;
true
gap> mat:= IdentityMat( 3, GF(257) );;
gap> g:= GO( 3, 257, mat );;
gap> InvariantQuadraticForm( g ).matrix = mat;
true

# Increase the code coverage.
gap> mat:= IdentityMat( 3, GF(5) );;
gap> _IsEqualModScalars( mat, Z(5) * mat );
true
gap> _IsEqualModScalars( mat, Zero( mat ) );
false
gap> _IsEqualModScalars( Zero( mat ), Zero( mat ) );
true
gap> _IsEqualModScalars( mat, IdentityMat( 2, GF(5) ) );
false
gap> _IsEqualModScalars( mat, NullMat( 3, 2, GF(5) ) );
false

##
gap> STOP_TEST( "classic.tst" );
