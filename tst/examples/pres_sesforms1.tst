gap> START_TEST("Forms: pres_sesforms1.tst");
gap> #Preserved Sesquilinear forms: example 1
gap> g := SU(4,3);
SU(4,3)
gap> forms := PreservedSesquilinearForms(g);
[ < hermitian form > ]
gap> Display( forms[1] );
Hermitian form
Gram Matrix:
 . . . 1
 . . 1 .
 . 1 . .
 1 . . .
gap> STOP_TEST("pres_sesforms1.tst", 10000 );
