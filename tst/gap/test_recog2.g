#case 1: DeltaOMinus(6,GF(9));

matrices := [ [ [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), Z(3^2)^3, Z(3^2), Z(3)^0, Z(3)^0 ], 
      [ 0*Z(3), 0*Z(3), Z(3^2), Z(3^2)^3, Z(3), Z(3) ], 
      [ 0*Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3^2)^3, Z(3^2)^5 ], 
      [ 0*Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3^2)^5, Z(3^2)^3 ] ], 
  [ [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3^2)^6, Z(3^2)^2 ], 
      [ Z(3^2)^7, Z(3^2), Z(3^2)^5, Z(3^2), Z(3^2), Z(3^2) ], 
      [ Z(3)^0, Z(3^2)^6, Z(3^2)^6, Z(3^2)^2, Z(3^2), Z(3^2) ], 
      [ Z(3^2)^7, Z(3^2)^6, Z(3^2)^5, Z(3^2), Z(3^2)^2, Z(3^2)^2 ], 
      [ Z(3^2)^6, 0*Z(3), Z(3^2)^5, Z(3^2)^3, 0*Z(3), 0*Z(3) ], 
      [ Z(3^2)^6, 0*Z(3), Z(3^2)^7, Z(3^2), 0*Z(3), 0*Z(3) ] ], 
  [ [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ Z(3^2), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), Z(3^2), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ], 
      [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3^2), 0*Z(3) ] ] ];
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

#case 2: DeltaOMinus(6,GF(49));

matrices := [ [ [ Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7^2)^47, 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7^2) ] ], 
  [ [ Z(7^2)^3, Z(7^2)^13, Z(7^2)^44, Z(7^2)^44, Z(7^2)^30, 0*Z(7) ], 
      [ Z(7^2)^3, Z(7^2)^13, Z(7^2)^20, Z(7^2)^20, Z(7^2)^30, 0*Z(7) ], 
      [ Z(7^2)^3, Z(7^2)^37, Z(7^2)^20, Z(7^2)^44, Z(7^2)^13, 0*Z(7) ], 
      [ Z(7^2)^27, Z(7^2)^13, Z(7^2)^20, Z(7^2)^44, Z(7^2)^37, 0*Z(7) ], 
      [ Z(7^2)^12, Z(7^2)^36, 0*Z(7), 0*Z(7), Z(7)^3, Z(7)^0 ], 
      [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7) ] ], 
  [ [ Z(7^2), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), Z(7^2), 0*Z(7), 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7^2), 0*Z(7) ], 
      [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ] ] ];
      
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

grp := DeltaOminus(4,GF(9));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

grp := DeltaOminus(8,GF(9));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

grp := DeltaOplus(6,GF(9));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

otherform := BilinearFormByQuadraticForm(form);
ScalarsOfPreservedForm(grp,otherform);


grp := Spdesargues(6,GF(11^4));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

grp := SOdesargues(-1,6,GF(9));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

grp := SOdesargues(0,7,GF(11));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

grp := SOdesargues(1,8,GF(16));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);


grp := SOdesargues(0,7,GF(11));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
otherform := BilinearFormByQuadraticForm(form);
ScalarsOfPreservedForm(grp,otherform);

grp := SUdesargues(4,GF(9));
gens := GeneratorsOfGroup(grp); 
matrices := List(gens,x->Unpack(x!.mat));
grp := Group(matrices);
forms := PreservedForms(grp);
form := forms[1];
ScalarsOfPreservedForm(grp,form);

otherform := BilinearFormByQuadraticForm(form);
ScalarsOfPreservedForm(grp,otherform);

