gap> START_TEST("Forms: test_pres_sesforms1.tst");
gap> g := SU(4,3);
SU(4,3)
gap> forms := PreservedSesquilinearForms(g);;
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> Display( forms[1] );
Hermitian form
Gram Matrix:
 . . . 1
 . . 1 .
 . 1 . .
 1 . . .
gap> STOP_TEST("test_pres_sesforms1.tst", 10000 );
