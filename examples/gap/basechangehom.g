#morphisms: BaseChangeHomomorphism
gl:=GL(3,3);
go:=GO(3,3);
form := PreservedSesquilinearForms(go)[1]; 
gram := GramMatrix( form );  
b := BaseChangeToCanonical(form);;
hom := BaseChangeHomomorphism(b, GF(3));
newgo := Image(hom, go); 
gens := GeneratorsOfGroup(newgo);;
canonical := b * gram * TransposedMat(b);
ForAll(gens, y -> y * canonical * TransposedMat(y) = canonical);
quit;
