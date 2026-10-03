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

# PossibleClassicalForms must decide on elements with nonzero traces too:
# over a large field these are almost all elements.
gap> g := DiagonalMat(List([2,3,1,1,1,1/6], x -> x*Z(257)^0));;
gap> forms := rec(maybeDual := true, maybeFrobenius := false, field := GF(257));;
gap> PossibleClassicalForms(SL(6,257), g, forms);;
gap> forms.maybeDual;
false
gap> forms := rec(maybeDual := true, maybeFrobenius := true, field := GF(257^2));;
gap> PossibleClassicalForms(SL(6,257^2), g, forms);;
gap> [forms.maybeDual, forms.maybeFrobenius];
[ false, false ]

# but must not rule out forms which are preserved
gap> ForAll(GeneratorsOfGroup(Sp(6,257)), function(g)
>      local forms;
>      forms := rec(maybeDual := true, maybeFrobenius := false, field := GF(257));
>      PossibleClassicalForms(Sp(6,257), g, forms);
>      return forms.maybeDual;
>    end);
true
gap> ForAll(GeneratorsOfGroup(GU(4,257)), function(g)
>      local forms;
>      forms := rec(maybeDual := true, maybeFrobenius := true, field := GF(257^2));
>      PossibleClassicalForms(GU(4,257), g, forms);
>      return forms.maybeFrobenius;
>    end);
true
gap> STOP_TEST("test_recog.tst", 10000 );
