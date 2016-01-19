#What is the form preserved by this group?
go := GO(5, 5);
x := 
[ [ Z(5)^0, Z(5)^3, 0*Z(5), Z(5)^3, Z(5)^3 ], 
  [ Z(5)^2, Z(5)^3, 0*Z(5), Z(5)^2, Z(5) ], 
  [ Z(5)^2, Z(5)^2, Z(5)^0, Z(5), Z(5)^3 ],
  [ Z(5)^0, Z(5)^3, Z(5), Z(5)^0, Z(5)^3 ], 
  [ Z(5)^3, 0*Z(5), Z(5)^0, 0*Z(5), Z(5) ] 
 ];;
go2 := go^x;
forms := PreservedSesquilinearForms( go2 );
Display( forms[1] );
quit;
