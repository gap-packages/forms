#What is the form preserved by this group?
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
TestPreservedSesquilinearForms(grp,forms);
quit;
