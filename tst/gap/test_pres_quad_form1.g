for q in [2,3,4,5,7,8,9] do
group := GO(-1,4,q);
pres_forms := PreservedQuadraticForms(group);
form := pres_forms[1];
gram := GramMatrix(form);
inv := InvariantQuadraticForm(group).matrix;
Print("q = ",q," equality: ",gram=inv,"\n");
#GramMatrix(pres_forms[1]) = InvariantQuadraticForm(group).matrix;
od;
