for q in [2,3,4,5,7,8,9] do
group := GO(-1,4,q);
pres_forms := PreservedQuadraticForms(group);
form := pres_forms[1];
gram := GramMatrix(form);
inv := InvariantQuadraticForm(group).matrix;
Print("q = ",q," equality: ",gram=inv,"\n");
#GramMatrix(pres_forms[1]) = InvariantQuadraticForm(group).matrix;
od;

pres_forms_quad := PreservedQuadraticForms(group);
form := pres_forms_quad[1];
Display(GramMatrix(form));
inv := InvariantQuadraticForm(group).matrix;
Display(inv);

pres_forms_ses := PreservedSesquilinearForms(group);
form := pres_forms_ses[1];
Display(GramMatrix(form));
inv := InvariantQuadraticForm(group).matrix;
Display(inv);
