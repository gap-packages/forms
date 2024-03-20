test := function(groups)
local grp, field, frob, output, newgens, bool;
bool := true;
for grp in groups do 
field := DefaultFieldOfMatrixGroup(grp);
frob := FrobeniusAutomorphism(field)^0;
output := ClassicalForms_GeneratorsWithBetterScalarsSesquilinear(grp,frob);
newgens := output[1];
bool := bool and Size(Group(newgens))= Size(grp);
od;
return bool;
end;

qs := [3,4,5,7,8,9];
ds := [4,6,8,10];

parameters := Cartesian(ds,qs);

for param in parameters do 
groups := ClassicalMaximalsGeneric("S",param[1],param[2]);
groups := List(groups,x->Group(GeneratorsOfGroup(x)));
Print(param," ",test(groups),"\n");
#test(groups);
od;

qs := [3,4,5];
ds := [4,6,8,10];

parameters := Cartesian(ds,qs);

for param in parameters do
for q in qs do
for d in ds do
groups := ClassicalMaximalsGeneric("S",d,q);
forms := [];
for grp in groups do
Add(forms,PreservedForms(grp));
od;
od;
od;

matrices := [];
for fs in forms do
if not IsEmpty(fs) then
Add(matrices,List(fs,x->GramMatrix(x)));
fi;
od;

#make a list of groups

qs := [3,4,5];
ds := [4,6,8,10];

#parameters := Cartesian(ds,qs);

many_groups := [];
for q in qs do
for d in ds do
Add(many_groups,ClassicalMaximalsGeneric("S",d,q));
od;
od;

groups := ClassicalMaximalsGeneric("S",4,5);  

test_preserved := function(grp,form)
local gens, frob, gram, field, p, i, e, q;
gens := GeneratorsOfGroup(grp);
gram := GramMatrix(form);
field := BaseField(form);
if IsHermitianForm(form) then
q := Size(field);
i := Sqrt(q);
p := Characteristic(field);
e := LogInt(i,p);
frob := FrobeniusAutomorphism(field)^e;
else
frob := FrobeniusAutomorphism(field)^0;
fi;
if IsSesquilinearForm(form) then
scalar_matrices := List(gens,x->gram^(-1)*x*gram*TransposedMat(x)^frob);
if not ForAll(scalar_matrices,IsDiagonalMat) then
    Print("not all matrices are diagonal\n");
    return false;
elif ForAll(List(scalar_matrices,x->Length(Set(DiagonalOfMat(x)))=1),x->x=true) then
    Print("all matrices are scalar\n");
    return true;
else
    Print("not all matrices are scalar\n");
    return false;
fi;
else
    Print("quadratic forms not yet ready\n");
    return false;
fi;
end;


group_forms_pairs := function(groups)
local pairs, grp, forms;
pairs := [];
for grp in groups do
forms := PreservedForms(grp);
Add(pairs,[grp,forms]);
od;
return pairs;
end;

test_pairs := function(pairs)
local grp,forms,pair,bad,test,form;
bad := [];
for pair in pairs do
grp := pair[1];
forms := pair[2];
if not Length(forms) = 0 then
for form in forms do
test := test_preserved(grp,form);
if not test then
Add(bad,[grp,form]);
fi;
od;
fi;
od;
return bad;
end;

groups := ClassicalMaximalsGeneric("S",4,3);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",4,4);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",4,5);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",4,7);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",4,8);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",4,9);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",4,11);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",4,16);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",6,3);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",6,4);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",6,5);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",6,7);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",6,8);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",6,11);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",6,16);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",8,4);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("S",10,4);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

#unitary gruups

groups := ClassicalMaximalsGeneric("U",3,3);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",3,4);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",3,5);  
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",3,7);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",3,8);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",3,9);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",3,11);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",4,3);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",4,4);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",4,5);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",4,7);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",4,8);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("U",4,9);
pairs := group_forms_pairs(groups);
bad := test_pairs(pairs);

group_ses_forms_pairs := function(groups)
local pairs, grp, forms;
pairs := [];
for grp in groups do
forms := PreservedSesquilinearForms(grp);
Add(pairs,[grp,forms]);
od;
return pairs;
end;


#orthogonal gruups

groups := ClassicalMaximalsGeneric("O",7,3);
pairs := group_ses_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("O-",8,3);
pairs := group_ses_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("O+",8,3);
pairs := group_ses_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("O+",8,9);
pairs := group_ses_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("O+",8,25);
pairs := group_ses_forms_pairs(groups);
bad := test_pairs(pairs);

groups := ClassicalMaximalsGeneric("O",7,5);
pairs := group_ses_forms_pairs(groups);
bad := test_pairs(pairs);




mat := [[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]]*Z(2)^0;

check_scalar_matrix := function(T,F)
local i,lambda;
i := PositionNonZero(T[1]);
if i <> PositionNonZero(F[1]) then
return false;
fi;
lambda := T[1][i] / F[1][i];
if T <> lambda*F then
return false;
else
return lambda;
fi;
end;






####stuff on induced action of module


grp := Group([ [ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), Z(2^3), 0*Z(2) ],
>       [ 0*Z(2), 0*Z(2), Z(2^3)^6 ] ],
>   [ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ Z(2)^0, Z(2)^0, Z(2)^0 ],
>       [ 0*Z(2), Z(2)^0, 0*Z(2) ] ] ]
> );

field := GF(8);

module := GModuleByMats(GeneratorsOfGroup(grp), field);

BaseSteinitzVectors(localouter, [localinner]).factorspace;
#localouter:= localouter := BasisVectors(Basis(ps!.vectorspace)); -> in our case always the identity matrix.



basis := MTX.BasesMinimalSubmodules(module);
basis := basis[1];
B := [];
B[1] := basis[1];
stein := BaseSteinitzVectors(One(grp), basis).factorspace;
inducedaction := MTX.InducedActionFactorModule(module,sub,stein);
subgens := inducedaction.generators;


B := Concatenation(B,BaseSteinitzVectors(One(grp), basis).factorspace);


sub := MTX.SubGModule(module,basis);
inducedaction := MTX.InducedActionFactorModule(module,sub);
subgens := inducedaction.generators;
grp_induced := Group(subgens);
forms := PreservedForms(grp_induced);
form := forms[1];
subformmat := GramMatrix(form);
zeromat := MutableCopyMat(One(grp)*0);
zeromat{[2..3]}{[2..3]} := subformmat;
newgrammat := Inverse(B)*zeromat*TransposedMat(Inverse(B));
newform := BilinearFormByMatrix(newgrammat,GF(8));
test_preserved(grp,form);

newgrammat := B*zeromat*TransposedMat(B);
newform := BilinearFormByMatrix(newgrammat,GF(8));
test_preserved(grp,form);

newgrammat := Inverse(B)*zeromat*TransposedMat(Inverse(B));
newform := BilinearFormByMatrix(newgrammat,GF(8));
test_preserved(grp,form);

newgrammat := TransposedMat(B)*zeromat*B;
newform := BilinearFormByMatrix(newgrammat,GF(8));
test_preserved(grp,form);

newgrammat := Inverse(TransposedMat(B))*zeromat*Inverse(B);
newform := BilinearFormByMatrix(newgrammat,GF(8));
test_preserved(grp,form);


