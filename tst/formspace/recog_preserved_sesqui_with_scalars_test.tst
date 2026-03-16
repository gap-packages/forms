gap> START_TEST("Forms: recog_preserved_sesqui_with_scalars_test.tst"); # Adapted from test_recog.tst for PreservedSesquilinearFormsWithScalars
gap> g := Sp(6,3);
Sp(6,3)
gap> forms := PreservedSesquilinearFormsWithScalars(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SU(4,3);
SU(4,3)
gap> forms := PreservedSesquilinearFormsWithScalars(g);
[ < hermitian form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(1,4,3);
SO(+1,4,3)
gap> forms := PreservedSesquilinearFormsWithScalars(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(1,4,4);
GO(+1,4,4)
gap> forms := PreservedSesquilinearFormsWithScalars(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(-1,4,3);
SO(-1,4,3)
gap> forms := PreservedSesquilinearFormsWithScalars(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(5,3);
SO(0,5,3)
gap> forms := PreservedSesquilinearFormsWithScalars(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> STOP_TEST("recog_preserved_sesqui_with_scalars_test.tst", 10000 );
