gap> START_TEST("Formspace: preserved_forms_non_trivial_scalars.tst"); 
gap> G := Sp(4, 5);;
gap> H := Group(G.1 * Z(5), G.2 * Z(5)^3);
<matrix group with 2 generators>
gap> L := PreservedForms(H);
[ < bilinear form > ]
gap> ScalarsOfPreservedForm(H, L[1]);
[ Z(5)^2, Z(5)^2 ]
gap> TestPreservedSesquilinearForms(H,L);
true
gap> G := SU(4, 5);;
gap> H := Group(G.1 * Z(5), G.2 * Z(5)^3);
<matrix group with 2 generators>
gap> L := PreservedForms(H);
[ < hermitian form > ]
gap> ScalarsOfPreservedForm(H, L[1]);
[ Z(5)^2, Z(5)^2 ]
gap> TestPreservedSesquilinearForms(H,L);
true
gap> STOP_TEST("Formspace: preserved_forms_non_trivial_scalars.tst", 10000 );
