#technical test, to make sure code is better covered
g := GL(4,5);
frob := FrobeniusAutomorphism(GF(5));
gens := GeneratorsOfGroup(g);
m := gens[1];
ClassicalForms_PossibleScalarsSesquilinear(GF(5),m,frob);
ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob);
m := gens[2];
ClassicalForms_PossibleScalarsSesquilinear(GF(5),m,frob);
ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob);
g := GL(5,9);
frob := FrobeniusAutomorphism(GF(9));
gens := GeneratorsOfGroup(g);
m := gens[1];
ClassicalForms_PossibleScalarsSesquilinear(GF(9),m,frob);
ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob);
ClassicalForms_PossibleScalarsSesquilinear(GF(9),m,frob^0);
ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(Group(m),frob^0);
quit;
