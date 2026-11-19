gap> START_TEST("Forms: matobj/bg_th_ex9.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([ [ Z(2^2), Z(2^2), Z(2^2), Z(2^2), Z(2^2) ], 
>    [ 0*Z(2), Z(2^2), Z(2^2)^2, 0*Z(2), Z(2)^0 ], 
>    [ 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2)^0 ], 
>    [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0 ], 
>    [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ] ], GF(4));;
gap> qform := QuadraticFormByMatrix(mat,GF(4));
< quadratic form >
gap> IsSingularForm(qform);
false
gap> IsDegenerateForm(qform);
#I  Testing degeneracy of the *associated bilinear form*
true
gap> biform := AssociatedBilinearForm(qform);
< bilinear form >
gap> Display(biform);
Bilinear form
Gram Matrix:
<immutable 5x5 plist matrix over GF(2^2):
[[ 0*Z(2), Z(2^2), Z(2^2), Z(2^2), Z(2^2) ]
 [ Z(2^2), 0*Z(2), Z(2^2)^2, 0*Z(2), Z(2)^0 ]
 [ Z(2^2), Z(2^2)^2, 0*Z(2), Z(2)^0, Z(2)^0 ]
 [ Z(2^2), 0*Z(2), Z(2)^0, 0*Z(2), Z(2)^0 ]
 [ Z(2^2), Z(2)^0, Z(2)^0, Z(2)^0, 0*Z(2) ]
]>
gap> IsDegenerateForm(biform);
true
gap> STOP_TEST("matobj/bg_th_ex9.tst", 10000 );
