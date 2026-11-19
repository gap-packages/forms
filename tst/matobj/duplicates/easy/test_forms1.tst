gap> START_TEST("Forms: matobj/test_forms1.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> f := GF(3);;
gap> gram := ToMatObj([[0,0,0,0,0,2],[0,0,0,0,2,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,2,0,0,0,0],[2,0,0,0,0,0]]*Z(3)^0, f);
<6x6 plist matrix over GF(3)>
gap> form := BilinearFormByMatrix(gram,f);
< bilinear form >
gap> Display(BaseChangeToCanonical(form));
<immutable 6x6 plist matrix over GF(3):
[[ Z(3)^0, 0*Z(3), Z(3)^0, Z(3)^0, 0*Z(3), Z(3)^0 ]
 [ Z(3)^0, 0*Z(3), Z(3), Z(3), 0*Z(3), Z(3)^0 ]
 [ Z(3)^0, Z(3), Z(3), Z(3)^0, 0*Z(3), Z(3) ]
 [ Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3), Z(3)^0 ]
 [ Z(3), 0*Z(3), Z(3)^0, Z(3), 0*Z(3), Z(3)^0 ]
 [ 0*Z(3), Z(3)^0, Z(3), Z(3)^0, Z(3), 0*Z(3) ]
]>
gap> Display(form);
Hyperbolic bilinear form
Gram Matrix:
<immutable 6x6 plist matrix over GF(3):
[[ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
]>
Witt Index: 3
gap> TypeOfForm(form);
1
gap> STOP_TEST("matobj/test_forms1.tst", 10000 );
