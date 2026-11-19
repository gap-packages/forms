gap> START_TEST("Formspace: preserved_quadratic_forms.tst");
gap> # symplectic groups preserve no quadratic form, neither in odd nor in even characteristic
gap> PreservedQuadraticForms(Sp(4, 3));
[  ]
gap> PreservedQuadraticForms(Sp(4, 4));
[  ]
gap> # orthogonal groups preserve exactly one quadratic form (up to scalars)
gap> L := PreservedQuadraticForms(GO(1, 4, 3));
[ < quadratic form > ]
gap> L := PreservedQuadraticForms(GO(-1, 4, 4));
[ < quadratic form > ]
gap> STOP_TEST("Formspace: preserved_quadratic_forms.tst", 10000 );
