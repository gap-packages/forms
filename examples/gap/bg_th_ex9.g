#Background theory: example 9
mat := [ [ Z(2^2), Z(2^2), Z(2^2), Z(2^2), Z(2^2) ], 
   [ 0*Z(2), Z(2^2), Z(2^2)^2, 0*Z(2), Z(2)^0 ], 
   [ 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2)^0 ], 
   [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0 ], 
   [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ] ];;
qform := QuadraticFormByMatrix(mat,GF(4));
IsSingularForm(qform);
IsDegenerateForm(qform);
biform := AssociatedBilinearForm(qform);
Display(biform);
IsDegenerateForm(biform);
quit;
