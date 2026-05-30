gap> START_TEST("Forms: recog_preserved_with_scalars_test.tst"); # Adapted from test_recog.tst for PreservedFormsWithScalars
gap> g := Sp(6,3);
Sp(6,3)
gap> forms := PreservedFormsWithScalars(g);
[ < bilinear form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SU(4,3);
SU(4,3)
gap> forms := PreservedFormsWithScalars(g);
[ < hermitian form > ]
gap> TestPreservedSesquilinearForms(g,forms);
true
gap> g := SO(1,4,3);
SO(+1,4,3)
gap> forms := PreservedFormsWithScalars(g);
[ < quadratic form > ]
gap> TestPreservedSesquilinearForms(g,[BilinearFormByQuadraticForm(forms[1])]);
true
gap> g := SO(1,4,4);
GO(+1,4,4)
gap> forms := PreservedFormsWithScalars(g); # Because of Char = 2 BilinearFormByQuadraticForm(forms[1]) will not work
[ < quadratic form > ]
gap> g := SO(-1,4,3);
SO(-1,4,3)
gap> forms := PreservedFormsWithScalars(g);
[ < quadratic form > ]
gap> TestPreservedSesquilinearForms(g,[BilinearFormByQuadraticForm(forms[1])]);
true
gap> g := SO(5,3);
SO(0,5,3)
gap> forms := PreservedFormsWithScalars(g);
[ < quadratic form > ]
gap> TestPreservedSesquilinearForms(g,[BilinearFormByQuadraticForm(forms[1])]);
true
gap> STOP_TEST("recog_preserved_with_scalars_test.tst", 10000 );
