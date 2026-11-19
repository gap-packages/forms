gap> START_TEST("Forms: matobj/bg_th_ex4.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,0,0,1],
>         [0,0,0,0,1,0],[0,0,0,1,0,0],[0,0,1,0,0,0]]*Z(16)^0, GF(16));
<6x6 plist matrix over GF(2^4)>
gap> form := BilinearFormByMatrix(mat,GF(16));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
<immutable 6x6 plist matrix over GF(2^4):
[[ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ]
 [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ]
]>
gap> IsSymmetricForm(form);
true
gap> IsAlternatingForm(form);
true
gap> IsDegenerateForm(form);
false
gap> WittIndex(form);
3
gap> STOP_TEST("matobj/bg_th_ex4.tst", 10000 );
