#Discriminant of form: DiscriminantOfForm
gram := InvariantQuadraticForm(GO(-1,4,5))!.matrix;;
qform := QuadraticFormByMatrix(gram, GF(5));
DiscriminantOfForm( qform );
quit;
