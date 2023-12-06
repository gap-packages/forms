gap> START_TEST("Forms: pres_sesforms2.tst");
gap> #Preserved Sesquilinear forms: example 2
gap> a := [ [ -1, 0, 0, -1, 0, 1 ], [ 0, -1, -1, 0, 0, 1 ], 
>        [ -1, 0, 0, 1, 0, 0 ],  [ 0, -1, 1, 0, 0, -1 ], 
>        [ 0, 0, 0, 0, 0, -1 ], [ 0, -1, -1, 1, 1, 1 ] ] * One(GF(3));;
gap> b := [ [ 1, -1, 1, -1, 1, -1 ], [ 1, 1, -1, 1, 1, 0 ], 
>        [ -1, 0, 1, 0, 0, 0 ], [ 0, -1, 0, 0, 0, 1 ], 
>        [ 1, 1, 1, 1, 1, 1 ], [ -1, 1, 1, 1, -1, 0 ] ] * One(GF(3));;
gap> g := Group( a, b );
<matrix group with 2 generators>
gap> forms := PreservedSesquilinearForms( g );
[ < bilinear form > ]
gap> Display( forms[1] );
Bilinear form
Gram Matrix:
 . 1 . . . .
 1 . . . . .
 . . . 1 . .
 . . 1 . . .
 . . . . . 1
 . . . . 1 .
gap> m := GModuleByMats( [a,b], GF(3) );;
gap> usemeataxe := MTX.InvariantBilinearForm(m);
fail
gap> STOP_TEST("pres_sesforms2.tst", 10000 );
