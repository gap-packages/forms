gap> START_TEST("Forms: test_pres_sesforms1.tst");
gap> g := SU(4,3);
SU(4,3)
gap> forms := PreservedSesquilinearForms(g);;
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> STOP_TEST("test_pres_sesforms1.tst", 10000 );
