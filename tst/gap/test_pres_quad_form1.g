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

e := -1;
d := 4;
q := 3;

for e in [-1,1] do
for d in [2,4,6,8,10] do
for q in [2,4,8,16,32,64,128,256,512,1024] do
group := GO(e,d,q);
invariant_quad_form := QuadraticFormByMatrix(InvariantQuadraticForm(group).matrix);
pres_forms := PreservedQuadraticForms(group);
form := pres_forms[1];
Print("e = ",e,", d = ",d,", q = ",q," equality: ",GramMatrix(form)=GramMatrix(invariant_quad_form),"\n");
od;
od;
od;

