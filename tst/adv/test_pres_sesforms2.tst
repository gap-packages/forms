gap> START_TEST("Forms: test_pres_sesforms2.tst");
gap> a := [ [ -1, 0, 0, -1, 0, 1 ], [ 0, -1, -1, 0, 0, 1 ], 
>        [ -1, 0, 0, 1, 0, 0 ],  [ 0, -1, 1, 0, 0, -1 ], 
>        [ 0, 0, 0, 0, 0, -1 ], [ 0, -1, -1, 1, 1, 1 ] ] * One(GF(3));;
gap> b := [ [ 1, -1, 1, -1, 1, -1 ], [ 1, 1, -1, 1, 1, 0 ], 
>        [ -1, 0, 1, 0, 0, 0 ], [ 0, -1, 0, 0, 0, 1 ], 
>        [ 1, 1, 1, 1, 1, 1 ], [ -1, 1, 1, 1, -1, 0 ] ] * One(GF(3));;
gap> g := Group( a, b );
<matrix group with 2 generators>
gap> forms := PreservedSesquilinearForms( g );;
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> m := GModuleByMats( [a,b], GF(3) );;
gap> usemeataxe := MTX.InvariantBilinearForm(m);
fail
gap> STOP_TEST("test_pres_sesforms2.tst", 10000 );
