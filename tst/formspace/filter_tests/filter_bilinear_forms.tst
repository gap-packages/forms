gap> START_TEST("Formspace: Filter Bilinear Forms"); # todo test char 2
gap> F := GF(5^2);;
gap> Forms := [[[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]*One(GF(5)), [[0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0]]*One(GF(5))];; # This should give a symmetric form and a syplectic form
gap> L := FORMS_FilterBilinearForms(Forms, F, 4);;
gap> Size(L);
2
gap> Size(L[1]);
1
gap> Size(L[2]);
1
gap> FORMS_IsSymmetricMatrix(L[1][1]);
true
gap> FORMS_IsSymplecticMatrix(L[2][1], GF(5^2));
true
gap> STOP_TEST("Formspace: Filter Bilinear Forms");
