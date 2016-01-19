#attributes and properties of a trivial form
mat := [[0,0,0],[0,0,0],[0,0,0]]*Z(11)^0;
form := QuadraticFormByMatrix(mat,GF(121));
IsReflexiveForm(form);
IsAlternatingForm(form);
IsSymmetricForm(form);
IsOrthogonalForm(form);
IsPseudoForm(form);
IsSymplecticForm(form);
IsDegenerateForm(form);
IsSingularForm(form);
BaseField(form);
GramMatrix(form);
RadicalOfForm(form);
quit;
