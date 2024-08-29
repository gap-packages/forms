gap> START_TEST("Forms: test_tech1.tst");
gap> g := GL(4,5);
GL(4,5)
gap> frob := FrobeniusAutomorphism(GF(5));
IdentityMapping( GF(5) )
gap> gens := GeneratorsOfGroup(g);
[ [ [ Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], [ 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5) ], 
      [ 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5) ], [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0 ] ]
    , 
  [ [ Z(5)^2, 0*Z(5), 0*Z(5), Z(5)^0 ], [ Z(5)^2, 0*Z(5), 0*Z(5), 0*Z(5) ], 
      [ 0*Z(5), Z(5)^2, 0*Z(5), 0*Z(5) ], [ 0*Z(5), 0*Z(5), Z(5)^2, 0*Z(5) ] 
     ] ]
gap> m := gens[1];
[ [ Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], [ 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5) ], [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0 ] ]
gap> ClassicalForms_PossibleScalarsSesquilinear(GF(5),m,frob);
false
gap> ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob);
false
gap> m := gens[2];
[ [ Z(5)^2, 0*Z(5), 0*Z(5), Z(5)^0 ], [ Z(5)^2, 0*Z(5), 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), Z(5)^2, 0*Z(5), 0*Z(5) ], [ 0*Z(5), 0*Z(5), Z(5)^2, 0*Z(5) ] ]
gap> ClassicalForms_PossibleScalarsSesquilinear(GF(5),m,frob);
false
gap> ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob);
false
gap> g := GL(5,9);
GL(5,9)
gap> frob := FrobeniusAutomorphism(GF(9));
FrobeniusAutomorphism( GF(3^2) )
gap> gens := GeneratorsOfGroup(g);
[ [ [ Z(3^2), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ] ], 
  [ [ Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ], 
      [ Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3) ] ] ]
gap> m := gens[1];
[ [ Z(3^2), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ] ]
gap> ClassicalForms_PossibleScalarsSesquilinear(GF(9),m,frob);
false
gap> ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob);
false
gap> ClassicalForms_PossibleScalarsSesquilinear(GF(9),m,frob^0);
false
gap> ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob^0);
false
gap> STOP_TEST("test_tech1.tst", 10000 );
