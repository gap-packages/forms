#technical test of ScalarsOfPreservedForm
go := GO(5, 5);
x := 
[ [ Z(5)^0, Z(5)^3, 0*Z(5), Z(5)^3, Z(5)^3 ], 
  [ Z(5)^2, Z(5)^3, 0*Z(5), Z(5)^2, Z(5) ], 
  [ Z(5)^2, Z(5)^2, Z(5)^0, Z(5), Z(5)^3 ],
  [ Z(5)^0, Z(5)^3, Z(5), Z(5)^0, Z(5)^3 ], 
  [ Z(5)^3, 0*Z(5), Z(5)^0, 0*Z(5), Z(5) ] 
 ];;
grp := go^x;
forms := PreservedSesquilinearForms( grp );;
ScalarsOfPreservedForm(grp,forms[1]);
ScalarsOfPreservedForm(go,forms[1]);
quad := QuadraticFormByBilinearForm(forms[1]);
ScalarsOfPreservedForm(grp,quad);
ScalarsOfPreservedForm(go,quad);
quit;

