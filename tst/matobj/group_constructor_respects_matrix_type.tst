#@local q, d, e, F, C, g, stored, form, calls, x
gap> START_TEST("Forms: /matobj/group_constructor_respects_matrix_type.tst");
gap> ReadPackage("forms", "tst/formspace/custom_test_functions.g");;
gap> # Also AI SLOP, but proved very usefull. This test tests the different constructors of classical groups to make sure that if they are provided with a IsPlistMatrixRep matrix group/form/matrix the produced group is again consisting of IsPlistMatrixRep matrices. It would be nice as a future change to also make the base group constructors maybe take a filter as an input to allow for specifying what type of group should be created.
gap> q:= 2;; d:= 4;; e:= -1;; F:= GF(q);;
gap> C:= GeneralOrthogonalGroup;;
gap> g:= TestForceMatrixObjGens( C( e, d, q ), IsPlistMatrixRep );;
gap> stored:= InvariantQuadraticForm( g ).matrix;;
gap> form:= QuadraticFormByMatrix( stored, F );;
gap> calls:= [ [ "C(g)",            {} -> C( g ) ],
>              [ "C(stored)",       {} -> C( stored ) ],
>              [ "C(form)",         {} -> C( form ) ],
>              [ "C(e,d,q,g)",      {} -> C( e, d, q, g ) ],
>              [ "C(e,d,q,stored)", {} -> C( e, d, q, stored ) ],
>              [ "C(e,d,q,form)",   {} -> C( e, d, q, form ) ],
>              [ "C(e,d,F,g)",      {} -> C( e, d, F, g ) ],
>              [ "C(e,d,F,stored)", {} -> C( e, d, F, stored ) ],
>              [ "C(e,d,F,form)",   {} -> C( e, d, F, form ) ] ];;
gap> for x in calls do
>      Print( x[1], ": ",
>             ConstructingFilter( GeneratorsOfGroup( x[2]() )[1] ), "\n" );
>    od;
C(g): <Representation "IsPlistMatrixRep">
C(stored): <Representation "IsPlistMatrixRep">
C(form): <Representation "IsPlistMatrixRep">
C(e,d,q,g): <Representation "IsPlistMatrixRep">
C(e,d,q,stored): <Representation "IsPlistMatrixRep">
C(e,d,q,form): <Representation "IsPlistMatrixRep">
C(e,d,F,g): <Representation "IsPlistMatrixRep">
C(e,d,F,stored): <Representation "IsPlistMatrixRep">
C(e,d,F,form): <Representation "IsPlistMatrixRep">
gap> ##Same for IsGenericMatrixRep
gap> g:= TestForceMatrixObjGens( C( e, d, q ), IsGenericMatrixRep );;
gap> stored:= InvariantQuadraticForm( g ).matrix;;
gap> form:= QuadraticFormByMatrix( stored, F );;
gap> calls:= [ [ "C(g)",            {} -> C( g ) ],
>              [ "C(stored)",       {} -> C( stored ) ],
>              [ "C(form)",         {} -> C( form ) ],
>              [ "C(e,d,q,g)",      {} -> C( e, d, q, g ) ],
>              [ "C(e,d,q,stored)", {} -> C( e, d, q, stored ) ],
>              [ "C(e,d,q,form)",   {} -> C( e, d, q, form ) ],
>              [ "C(e,d,F,g)",      {} -> C( e, d, F, g ) ],
>              [ "C(e,d,F,stored)", {} -> C( e, d, F, stored ) ],
>              [ "C(e,d,F,form)",   {} -> C( e, d, F, form ) ] ];;
gap> for x in calls do
>      Print( x[1], ": ",
>             ConstructingFilter( GeneratorsOfGroup( x[2]() )[1] ), "\n" );
>    od;
C(g): <Representation "IsGenericMatrixRep">
C(stored): <Representation "IsGenericMatrixRep">
C(form): <Representation "IsGenericMatrixRep">
C(e,d,q,g): <Representation "IsGenericMatrixRep">
C(e,d,q,stored): <Representation "IsGenericMatrixRep">
C(e,d,q,form): <Representation "IsGenericMatrixRep">
C(e,d,F,g): <Representation "IsGenericMatrixRep">
C(e,d,F,stored): <Representation "IsGenericMatrixRep">
C(e,d,F,form): <Representation "IsGenericMatrixRep">
gap> ############ Now for SpecialOrthogonalGroup and ############
gap> q:= 2;; d:= 4;; e:= -1;; F:= GF(q);;
gap> C:= SpecialOrthogonalGroup;;
gap> g:= TestForceMatrixObjGens( C( e, d, q ), IsPlistMatrixRep );;
gap> stored:= InvariantQuadraticForm( g ).matrix;;
gap> form:= QuadraticFormByMatrix( stored, F );;
gap> calls:= [ [ "C(g)",            {} -> C( g ) ],
>              [ "C(stored)",       {} -> C( stored ) ],
>              [ "C(form)",         {} -> C( form ) ],
>              [ "C(e,d,q,g)",      {} -> C( e, d, q, g ) ],
>              [ "C(e,d,q,stored)", {} -> C( e, d, q, stored ) ],
>              [ "C(e,d,q,form)",   {} -> C( e, d, q, form ) ],
>              [ "C(e,d,F,g)",      {} -> C( e, d, F, g ) ],
>              [ "C(e,d,F,stored)", {} -> C( e, d, F, stored ) ],
>              [ "C(e,d,F,form)",   {} -> C( e, d, F, form ) ] ];;
gap> for x in calls do
>      Print( x[1], ": ",
>             ConstructingFilter( GeneratorsOfGroup( x[2]() )[1] ), "\n" );
>    od;
C(g): <Representation "IsPlistMatrixRep">
C(stored): <Representation "IsPlistMatrixRep">
C(form): <Representation "IsPlistMatrixRep">
C(e,d,q,g): <Representation "IsPlistMatrixRep">
C(e,d,q,stored): <Representation "IsPlistMatrixRep">
C(e,d,q,form): <Representation "IsPlistMatrixRep">
C(e,d,F,g): <Representation "IsPlistMatrixRep">
C(e,d,F,stored): <Representation "IsPlistMatrixRep">
C(e,d,F,form): <Representation "IsPlistMatrixRep">
gap> ## Same for IsGenericMatrixRep
gap> g:= TestForceMatrixObjGens( C( e, d, q ), IsGenericMatrixRep );;
gap> stored:= InvariantQuadraticForm( g ).matrix;;
gap> form:= QuadraticFormByMatrix( stored, F );;
gap> calls:= [ [ "C(g)",            {} -> C( g ) ],
>              [ "C(stored)",       {} -> C( stored ) ],
>              [ "C(form)",         {} -> C( form ) ],
>              [ "C(e,d,q,g)",      {} -> C( e, d, q, g ) ],
>              [ "C(e,d,q,stored)", {} -> C( e, d, q, stored ) ],
>              [ "C(e,d,q,form)",   {} -> C( e, d, q, form ) ],
>              [ "C(e,d,F,g)",      {} -> C( e, d, F, g ) ],
>              [ "C(e,d,F,stored)", {} -> C( e, d, F, stored ) ],
>              [ "C(e,d,F,form)",   {} -> C( e, d, F, form ) ] ];;
gap> for x in calls do
>      Print( x[1], ": ",
>             ConstructingFilter( GeneratorsOfGroup( x[2]() )[1] ), "\n" );
>    od;
C(g): <Representation "IsGenericMatrixRep">
C(stored): <Representation "IsGenericMatrixRep">
C(form): <Representation "IsGenericMatrixRep">
C(e,d,q,g): <Representation "IsGenericMatrixRep">
C(e,d,q,stored): <Representation "IsGenericMatrixRep">
C(e,d,q,form): <Representation "IsGenericMatrixRep">
C(e,d,F,g): <Representation "IsGenericMatrixRep">
C(e,d,F,stored): <Representation "IsGenericMatrixRep">
C(e,d,F,form): <Representation "IsGenericMatrixRep">
gap> ############### Now for Omega
gap> q:= 2;; d:= 4;; e:= -1;; F:= GF(q);;
gap> C:= Omega;;
gap> g:= TestForceMatrixObjGens( C( e, d, q ), IsPlistMatrixRep );;
gap> stored:= InvariantQuadraticForm( g ).matrix;;
gap> form:= QuadraticFormByMatrix( stored, F );;
gap> calls:= [ [ "C(g)",            {} -> C( g ) ],
>              [ "C(stored)",       {} -> C( stored ) ],
>              [ "C(form)",         {} -> C( form ) ],
>              [ "C(e,d,q,g)",      {} -> C( e, d, q, g ) ],
>              [ "C(e,d,q,stored)", {} -> C( e, d, q, stored ) ],
>              [ "C(e,d,q,form)",   {} -> C( e, d, q, form ) ],
>              [ "C(e,d,F,g)",      {} -> C( e, d, F, g ) ],
>              [ "C(e,d,F,stored)", {} -> C( e, d, F, stored ) ],
>              [ "C(e,d,F,form)",   {} -> C( e, d, F, form ) ] ];;
gap> for x in calls do
>      Print( x[1], ": ",
>             ConstructingFilter( GeneratorsOfGroup( x[2]() )[1] ), "\n" );
>    od;
C(g): <Representation "IsPlistMatrixRep">
C(stored): <Representation "IsPlistMatrixRep">
C(form): <Representation "IsPlistMatrixRep">
C(e,d,q,g): <Representation "IsPlistMatrixRep">
C(e,d,q,stored): <Representation "IsPlistMatrixRep">
C(e,d,q,form): <Representation "IsPlistMatrixRep">
C(e,d,F,g): <Representation "IsPlistMatrixRep">
C(e,d,F,stored): <Representation "IsPlistMatrixRep">
C(e,d,F,form): <Representation "IsPlistMatrixRep">
gap> #Same for IsGenericMatrixRep
gap> g:= TestForceMatrixObjGens( C( e, d, q ), IsGenericMatrixRep );;
gap> stored:= InvariantQuadraticForm( g ).matrix;;
gap> form:= QuadraticFormByMatrix( stored, F );;
gap> calls:= [ [ "C(g)",            {} -> C( g ) ],
>              [ "C(stored)",       {} -> C( stored ) ],
>              [ "C(form)",         {} -> C( form ) ],
>              [ "C(e,d,q,g)",      {} -> C( e, d, q, g ) ],
>              [ "C(e,d,q,stored)", {} -> C( e, d, q, stored ) ],
>              [ "C(e,d,q,form)",   {} -> C( e, d, q, form ) ],
>              [ "C(e,d,F,g)",      {} -> C( e, d, F, g ) ],
>              [ "C(e,d,F,stored)", {} -> C( e, d, F, stored ) ],
>              [ "C(e,d,F,form)",   {} -> C( e, d, F, form ) ] ];;
gap> for x in calls do
>      Print( x[1], ": ",
>             ConstructingFilter( GeneratorsOfGroup( x[2]() )[1] ), "\n" );
>    od;
C(g): <Representation "IsGenericMatrixRep">
C(stored): <Representation "IsGenericMatrixRep">
C(form): <Representation "IsGenericMatrixRep">
C(e,d,q,g): <Representation "IsGenericMatrixRep">
C(e,d,q,stored): <Representation "IsGenericMatrixRep">
C(e,d,q,form): <Representation "IsGenericMatrixRep">
C(e,d,F,g): <Representation "IsGenericMatrixRep">
C(e,d,F,stored): <Representation "IsGenericMatrixRep">
C(e,d,F,form): <Representation "IsGenericMatrixRep">
gap> STOP_TEST("Forms: /matobj/group_constructor_respects_matrix_type.tst");
