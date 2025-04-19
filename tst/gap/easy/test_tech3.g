#technical test that will become obsolete in the future, output is even not important
g := Sp(6,3);
module := GModuleByMats(GeneratorsOfGroup(g),GF(3));;
dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(g);;
ClassicalForms_InvariantFormDual(module,dmodule);;
g := GO(3,3);
module := GModuleByMats(GeneratorsOfGroup(g),GF(3));;
dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(g);;
ClassicalForms_InvariantFormDual(module,dmodule);;
quit;

