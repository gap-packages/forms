#@local is_equal, q, F, d, es, e, g, filters, filt, stored, pi, permmat
#@local form, gg, F2, mat

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
gap> if IsBound( ConformalSymplecticGroup ) then
>      # Support for matrix objects was added together with this function,
>      # https://github.com/gap-system/gap/pull/6213.
>      # Once we know a GAP version that decides the availability,
>      # the version number can be used for the distinction.
>      filters:= [ IsPlistRep, IsPlistMatrixRep ];;
>    else
>      filters:= [ IsPlistRep ];;
>    fi;
gap> for q in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25 ] do
>      F:= GF(q);
>      for d in [ 2, 4 .. 8 ] do
>        for filt in filters do
>          PushOptions( rec( ConstructingFilter:= filt ) );
> 
>          g:= SymplecticGroup( d, q );
>          if filt <> IsPlistRep and not filt( One( g ) ) then
>            Error( "wrong repres. of matrices", [ q, d, filt ] );
>          fi;
>          stored:= InvariantBilinearForm( g ).matrix;
>          if filt <> IsPlistRep and not filt( stored ) then
>            Error( "wrong repres. of matrices" );
>          fi;
>          pi:= Matrix( PermutationMat( (1,2), d, F ), stored );
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= BilinearFormByMatrix( stored, F );
>          gg:= SymplecticGroup( d, q, permmat );
>          if filt <> IsPlistRep and not filt( One( gg ) ) then
>            Error( "wrong repres. of matrices" );
>          fi;
>          if not ( is_equal( g, SymplecticGroup( g ) ) and
>                   ( is_equal( g, SymplecticGroup( stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, SymplecticGroup( form ) ) and
>                   is_equal( g, SymplecticGroup( d, q, g ) ) and
>                   is_equal( g, SymplecticGroup( d, q, stored ) ) and
>                   is_equal( g, SymplecticGroup( d, q, form ) ) and
>                   is_equal( g, SymplecticGroup( d, F, g ) ) and
>                   is_equal( g, SymplecticGroup( d, F, stored ) ) and
>                   is_equal( g, SymplecticGroup( d, F, form ) ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with Sp(", d, ",", q, ") for ", filt );
>          fi;
> 
>          PushOptions( rec( ConstructingFilter:= fail ) );
>        od;
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

# Test inconsistent fields of definition.
gap> g:= GeneralOrthogonalGroup( 5, GF(3) );;
gap> mat:= InvariantQuadraticForm( g ).matrix;;
gap> form:= QuadraticFormByMatrix( mat, GF(9) );;
gap> GeneralOrthogonalGroup( 5, GF(3), form );;
Error, the defining field of <form> does not fit to <gf>
gap> g:= GeneralUnitaryGroup( 4, 2 );;
gap> mat:= InvariantSesquilinearForm( g ).matrix;;
gap> form:= HermitianFormByMatrix( mat, GF(16) );;
gap> GeneralUnitaryGroup( 4, 2, form );
Error, the defining field of <form> does not fit to <q>
gap> g:= SymplecticGroup( 4, GF(2) );;
gap> mat:= InvariantBilinearForm( g ).matrix;;
gap> form:= BilinearFormByMatrix( mat, GF(4) );;
gap> SymplecticGroup( 4, GF(2), form );
Error, the defining field of <form> does not fit to <q>

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
gap> STOP_TEST( "Forms: classic.tst" );
