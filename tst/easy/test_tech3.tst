gap> START_TEST("Forms: test_tech3.tst");
gap> g := Sp(6,3);
Sp(6,3)
gap> module := GModuleByMats(GeneratorsOfGroup(g),GF(3));;
gap> dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(g);;
gap> ClassicalForms_InvariantFormDual(module,dmodule);;
gap> g := GO(3,3);
GO(0,3,3)
gap> module := GModuleByMats(GeneratorsOfGroup(g),GF(3));;
gap> dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(g);;
gap> ClassicalForms_InvariantFormDual(module,dmodule);;
gap> STOP_TEST("test_tech3.tst", 10000 );
