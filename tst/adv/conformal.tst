#@local q, F, d, filts, filt, g, stored, pi, permmat, form, gg, pg, sp

gap> START_TEST( "Forms: conformal.tst" );

# Test the creation of conformal symplectic groups.
gap> for q in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25 ] do
>      F:= GF(q);
>      for d in [ 2, 4 .. 8 ] do
>        filts:= [ IsPlistRep, IsPlistMatrixRep ];
>        if q = 2 then
>          Add( filts, IsGF2MatrixRep );
>        elif q <= 256 then
>          Add( filts, Is8BitMatrixRep );
>        fi;
>        for filt in filts do
>          PushOptions( rec( ConstructingFilter:= filt ) );
> 
>          g:= ConformalSymplecticGroup( d, F );
>          if filt <> IsPlistRep and not filt( One( g ) ) then
>            Error( "wrong repres. of matrices", [ q, d, filt ] );
>          fi;
>          stored:= InvariantBilinearFormUpToScalars( g ).matrix;
>          if filt <> IsPlistRep and not filt( stored ) then
>            Error( "wrong repres. of matrices" );
>          fi;
>          pi:= Matrix( PermutationMat( (1,2), d, F ), stored );
>          permmat:= pi * stored * TransposedMat( pi );
>          form:= BilinearFormByMatrix( stored, F );
>          gg:= ConformalSymplecticGroup( d, F, permmat );
>          if filt <> IsPlistRep and not filt( One( gg ) ) then
>            Error( "wrong repres. of matrices" );
>          fi;
>          if not ( g = ConformalSymplecticGroup( g ) and
>                   ( g = ConformalSymplecticGroup( stored ) or
>                     BaseDomain( stored ) <> F ) and
>                   g = ConformalSymplecticGroup( d, q ) and
>                   g = ConformalSymplecticGroup( form ) and
>                   g = ConformalSymplecticGroup( d, q, g ) and
>                   g = ConformalSymplecticGroup( d, q, stored ) and
>                   g = ConformalSymplecticGroup( d, q, form ) and
>                   g = ConformalSymplecticGroup( d, F, g ) and
>                   g = ConformalSymplecticGroup( d, F, stored ) and
>                   g = ConformalSymplecticGroup( d, F, form ) and
>                   IsSubset( gg, GeneratorsOfGroup( gg ) ) and
>                   IsSubset( g, List( GeneratorsOfGroup( gg ), x -> x^pi ) ) ) then
>            Error( "problem with CSp(", d, ",", q, ")" );
>          fi;
> 
>          if Size( g ) < 10^7 then
>            pg:= ConformalSymplecticGroup( IsPermGroup, d, q );
>            if Size( g ) <> Size( pg ) then
>              Error( "problem with CSp(IsPermGroup, ", d, ",", q, ")" );
>            fi;
>          fi;
> 
>          sp:= SymplecticGroup( d, q, permmat );
>          if filt <> IsPlistRep and not filt( One( sp ) ) then
>            Error( "wrong repres. of matrices" );
>          fi;
>          g:= ConformalSymplecticGroup( sp );
>          if not IsSubset( g, sp ) then
>            Error( "problem with CSp(", d, ",", q, ")" );
>          fi;
> 
>          PushOptions( rec( ConstructingFilter:= fail ) );
>        od;
>      od;
>    od;

##
gap> STOP_TEST( "Forms: conformal.tst" );

