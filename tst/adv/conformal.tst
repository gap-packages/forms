#@local is_equal, q, F, d, filt, g, stored, pi, permmat, form, gg, pg

gap> START_TEST( "Forms: conformal.tst" );

# Provide an auxiliary function (until GAP's '=' gets fast).
gap> is_equal:= function( G1, G2 )
>      return IsSubset( G1, GeneratorsOfGroup( G2 ) ) and
>             IsSubset( G2, GeneratorsOfGroup( G1 ) );
>    end;;

# Test the creation of conformal symplectic groups.
gap> for q in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25 ] do
>      F:= GF(q);
>      for d in [ 2, 4 .. 8 ] do
>        for filt in [ IsPlistRep, IsPlistMatrixRep ] do
>          PushOptions( rec( ConstructingFilter:= filt ) );
> 
>          g:= ConformalSymplecticGroup( d, q );
>          stored:= InvariantBilinearFormUpToScalars( g ).matrix;
>          pi:= Matrix( PermutationMat( (1,2), d, F ), stored );
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= BilinearFormByMatrix( stored, F );
>          gg:= ConformalSymplecticGroup( d, q, permmat );
>          if not ( is_equal( g, ConformalSymplecticGroup( g ) ) and
>                   ( is_equal( g, ConformalSymplecticGroup( stored ) ) or
>                     BaseDomain( stored ) <> F ) and
>                   is_equal( g, ConformalSymplecticGroup( d, F ) ) and
>                   is_equal( g, ConformalSymplecticGroup( form ) ) and
>                   is_equal( g, ConformalSymplecticGroup( d, q, g ) ) and
>                   is_equal( g, ConformalSymplecticGroup( d, q, stored ) ) and
>                   is_equal( g, ConformalSymplecticGroup( d, q, form ) ) and
>                   is_equal( g, ConformalSymplecticGroup( d, F, g ) ) and
>                   is_equal( g, ConformalSymplecticGroup( d, F, stored ) ) and
>                   is_equal( g, ConformalSymplecticGroup( d, F, form ) ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with CSp(", d, ",", q, ")" );
>          fi;
> 
>          if Size( g ) < 10^7 and filt = IsPlistRep then
> #TODO: Make this work for `IsPlistMatrixRep`
>            pg:= ConformalSymplecticGroup( IsPermGroup, d, q );
>            if Size( g ) <> Size( pg ) then
>              Error( "problem with CSp(IsPermGroup, ", d, ",", q, ")" );
>            fi;
>          fi;
> 
> #TODO: Once `SymplecticGroup` supports matrix objects,
> #      the following will work.
> #        sp:= SymplecticGroup( d, q, permmat );
> #        g:= ConformalSymplecticGroup( sp );
> #        if not IsSubset( g, sp ) then
> #          Error( "problem with CSp(", d, ",", q, ")" );
> #        fi;
> 
>          PushOptions( rec( ConstructingFilter:= fail ) );
>        od;
>      od;
>    od;

##
gap> STOP_TEST( "conformal.tst" );

