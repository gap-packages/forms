#example of ScalarsOfPreservedForm
G := Sp(4,5);
group := Group(G.1 * Z(5), G.2 * Z(5)^3);
forms := PreservedForms(group);
ScalarsOfPreservedForm(group, forms[1]);
form := BilinearFormByMatrix(IdentityMat(4, GF(5)), GF(5));
ScalarsOfPreservedForm(group, form);
quit;
