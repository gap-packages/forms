gap> START_TEST("Formspace: preserved_quadratic_forms_scaled_generators.tst");
gap> # char 2 groups preserving a quadratic form modulo scalars that are not 1
gap> s := DiagonalMat([Z(4), 1, Z(4), 1] * Z(4)^0);;
gap> H := Group(Concatenation(GeneratorsOfGroup(SO(1, 4, 4)), [s]));;
gap> L := PreservedQuadraticForms(H);
[ < quadratic form > ]
gap> ScalarsOfPreservedForm(H, L[1]);
[ Z(2)^0, Z(2)^0, Z(2^2) ]
gap> H := Group(List(GeneratorsOfGroup(SO(1, 4, 4)), g -> g * s));;
gap> L := PreservedQuadraticForms(H);
[ < quadratic form > ]
gap> ScalarsOfPreservedForm(H, L[1]);
[ Z(2^2), Z(2^2) ]
gap> s := [[Z(4), 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, Z(4)^2]] * Z(4)^0;;
gap> H := Group(Concatenation(GeneratorsOfGroup(SO(-1, 4, 4)), [s]));;
gap> L := PreservedQuadraticForms(H);
[ < quadratic form > ]
gap> ScalarsOfPreservedForm(H, L[1]);
[ Z(2)^0, Z(2)^0, Z(2^2) ]
gap> s := DiagonalMat([Z(8)^3, 1, Z(8)^3, 1, Z(8)^3, 1] * Z(8)^0);;
gap> H := Group(Concatenation(GeneratorsOfGroup(SO(1, 6, 8)), [s]));;
gap> L := PreservedQuadraticForms(H);
[ < quadratic form > ]
gap> ScalarsOfPreservedForm(H, L[1]);
[ Z(2)^0, Z(2)^0, Z(2)^0, Z(2^3)^3 ]
gap> STOP_TEST("Formspace: preserved_quadratic_forms_scaled_generators.tst", 10000);
