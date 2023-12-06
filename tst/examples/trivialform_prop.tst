gap> START_TEST("Forms: trivialform_prop.tst");
gap> #attributes and properties of a trivial form
gap> mat := [[0,0,0],[0,0,0],[0,0,0]]*Z(11)^0;
[ [ 0*Z(11), 0*Z(11), 0*Z(11) ], [ 0*Z(11), 0*Z(11), 0*Z(11) ], 
  [ 0*Z(11), 0*Z(11), 0*Z(11) ] ]
gap> form := QuadraticFormByMatrix(mat,GF(121));
< trivial form >
gap> IsReflexiveForm(form);
true
gap> IsAlternatingForm(form);
true
gap> IsSymmetricForm(form);
true
gap> IsOrthogonalForm(form);
false
gap> IsPseudoForm(form);
false
gap> IsSymplecticForm(form);
true
gap> IsDegenerateForm(form);
true
gap> IsSingularForm(form);
true
gap> BaseField(form);
GF(11^2)
gap> GramMatrix(form);
[ [ 0*Z(11), 0*Z(11), 0*Z(11) ], [ 0*Z(11), 0*Z(11), 0*Z(11) ], 
  [ 0*Z(11), 0*Z(11), 0*Z(11) ] ]
gap> Dimension(RadicalOfForm(form));
3
gap> STOP_TEST("trivialform_prop.tst", 10000 );
