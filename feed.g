field := GF(64);
grp := GO(1,6,GF(64));
gens := GeneratorsOfGroup(grp);

matrixses := [ [ 0*Z(2), Z(2^6)^23, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ],
  [ Z(2^6)^23, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ],
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^6)^23, 0*Z(2), 0*Z(2) ],
  [ 0*Z(2), 0*Z(2), Z(2^6)^23, 0*Z(2), 0*Z(2), 0*Z(2) ],
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2^6)^23 ],
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2^6)^23, 0*Z(2) ] ];
  
List([1..Length(gens)],i->gens[i]*matrixses*TransposedMat(gens[i]) - matrixses = ZeroOp(matrixses));

mus := [Z(2^4)^5,Z(2^4)^23,Z(2^4)^87];
#scalars2 := [1,1,1]*Z(2)^0;
newgens := [];
for i in [1..Length(gens)] do
mu := mus[i];
Add(newgens,Unpack(gens[i])*mu);
od;
scalars := List(mus,x->x^2);

List([1..Length(newgens)],i->newgens[i]*matrixses*TransposedMat(newgens[i]) - scalars[i]*matrixses = ZeroOp(matrixses));
newgens[1]*matrixses*TransposedMat(newgens[1]);
gens[1]*matrixses*TransposedMat(gens[1]);


gens := newgens;

matrix := [ [ 0*Z(2), Z(2^6)^60, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ],
  [ Z(2^6)^60, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ],
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2^6)^60, 0*Z(2), 0*Z(2) ],
  [ 0*Z(2), 0*Z(2), Z(2^6)^60, 0*Z(2), 0*Z(2), 0*Z(2) ],
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2^6)^60 ],
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2^6)^60, 0*Z(2) ] ];
  
form := matrix;

quad := ClassicalForms_QuadraticForm2(field,matrix,newgens,scalars);






##### to reproduce bug

grp := GO(1,6,GF(64));
gens := GeneratorsOfGroup(grp);
g1 := gens[1];

mu := Z(4);
TransposedMat(g1*mu);
Display(last);
TransposedMat(g1)*mu;
Display(last);


mus := [Z(2^4)^5,Z(2^4)^23,Z(2^4)^87];
TransposedMat(g1*mus[1]);
Display(last);
TransposedMat(g1)*mus[1];
Display(last);


