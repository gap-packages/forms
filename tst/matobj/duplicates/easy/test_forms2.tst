gap> START_TEST("Forms: matobj/test_forms2.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> f := GF(7);
GF(7)
gap> gram := ToMatObj([[-3,0,0,0,0,0],[0,0,3,0,0,0],[0,3,0,0,0,0],[0,0,0,0,0,-1/2],[0,0,0,0,1,0],[0,0,0,-1/2,0,0]]*Z(7)^0, f);
<6x6 plist matrix over GF(7)>
gap> form := BilinearFormByMatrix(gram,f);
< bilinear form >
gap> IsEllipticForm(form);
true
gap> TypeOfForm(form);
-1
gap> Display(form);
Elliptic bilinear form
Gram Matrix:
<immutable 6x6 plist matrix over GF(7):
[[ Z(7)^4, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ]
 [ 0*Z(7), 0*Z(7), Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ]
 [ 0*Z(7), Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ]
 [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7) ]
 [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7) ]
 [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7), 0*Z(7), 0*Z(7) ]
]>
Witt Index: 2
gap> STOP_TEST("matobj/test_forms2.tst", 10000 );
