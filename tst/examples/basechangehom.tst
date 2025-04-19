gap> START_TEST("Forms: basechangehom.tst");
gap> gl:=GL(3,3);
GL(3,3)
gap> go:=GO(3,3);
GO(0,3,3)
gap> form := PreservedSesquilinearForms(go)[1]; 
< bilinear form >
gap> gram := GramMatrix( form );  
[ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ Z(3)^0, 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), Z(3) ] ]
gap> b := BaseChangeToCanonical(form);;
gap> hom := BaseChangeHomomorphism(b, GF(3));
^[ [ 0*Z(3), 0*Z(3), Z(3) ], [ Z(3), Z(3), Z(3)^0 ], [ Z(3), 0*Z(3), Z(3) ] ]
gap> newgo := Image(hom, go); 
Group(
[ 
  [ [ Z(3)^0, 0*Z(3), Z(3) ], [ Z(3)^0, Z(3), Z(3)^0 ], 
      [ 0*Z(3), 0*Z(3), Z(3) ] ], 
  [ [ Z(3), Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3), 0*Z(3) ], 
      [ Z(3)^0, Z(3)^0, Z(3) ] ] ])
gap> gens := GeneratorsOfGroup(newgo);;
gap> canonical := b * gram * TransposedMat(b);
[ [ Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), 0*Z(3), Z(3)^0 ], 
  [ 0*Z(3), Z(3)^0, 0*Z(3) ] ]
gap> ForAll(gens, y -> y * canonical * TransposedMat(y) = canonical);
true
gap> STOP_TEST("basechangehom.tst", 10000 );
