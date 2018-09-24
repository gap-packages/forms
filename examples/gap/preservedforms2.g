#example (due to John) of more than one preserved form (after three calls).
gens := [ [ [ Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5) ], 
	[ 0*Z(5), 0*Z(5), Z(5)^3, Z(5^2)^21 ],
	[ 0*Z(5), Z(5), Z(5), Z(5^2)^3 ], 
     	[ 0*Z(5), Z(5^2)^21, Z(5^2)^15, Z(5)^2 ] ], 
  [ [ Z(5)^3, Z(5^2)^7, Z(5^2)^16, Z(5^2)^15 ], 
      	[ 0*Z(5), Z(5)^0, Z(5^2)^4, Z(5)^3 ], 
      	[ Z(5^2)^22, Z(5^2)^10, Z(5^2)^21, Z(5)^2 ], 
      	[ Z(5^2)^7, Z(5^2)^23, Z(5^2)^9, Z(5^2)^11 ] ], 
  [ [ 0*Z(5), Z(5^2), 0*Z(5), 0*Z(5) ], 
	[ Z(5^2)^5, 0*Z(5), 0*Z(5), 0*Z(5) ], 
      	[ 0*Z(5), 0*Z(5), Z(5)^0, Z(5^2)^4 ], 
      	[ 0*Z(5), 0*Z(5), Z(5^2)^20, Z(5)^2 ] ] ];
group := Group(gens);
PreservedForms(group);
PreservedForms(group);
PreservedForms(group);
PreservedForms(group);
quit;
