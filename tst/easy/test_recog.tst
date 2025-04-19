gap> START_TEST("Forms: test_recog.tst");
gap> g := Sp(6,3);
Sp(6,3)
gap> forms := PreservedSesquilinearForms(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SU(4,3);
SU(4,3)
gap> forms := PreservedSesquilinearForms(g);
[ < hermitian form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(1,4,3);
SO(+1,4,3)
gap> forms := PreservedSesquilinearForms(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(1,4,4);
GO(+1,4,4)
gap> forms := PreservedSesquilinearForms(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(-1,4,3);
SO(-1,4,3)
gap> forms := PreservedSesquilinearForms(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(5,3);
SO(0,5,3)
gap> forms := PreservedSesquilinearForms(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> STOP_TEST("test_recog.tst", 10000 );
