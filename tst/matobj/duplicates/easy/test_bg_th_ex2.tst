gap> START_TEST("Forms: matobj/test_bg_th_ex2.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]]*Z(9)^0, GF(9));
<4x4 plist matrix over GF(3^2)>
gap> form := BilinearFormByMatrix(mat,GF(9));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
<immutable 4x4 plist matrix over GF(3^2):
[[ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3) ]
]>
gap> IsReflexiveForm(form);
true
gap> IsSymmetricForm(form);
true
gap> IsAlternatingForm(form);
false
gap> r := RadicalOfForm(form);;
gap> Dimension(r);
0
gap> STOP_TEST("matobj/test_bg_th_ex2.tst", 10000 );
