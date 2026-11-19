gap> START_TEST("Forms: matobj/test_forms5.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> f := GF(9);
GF(3^2)
gap> gram := ToMatObj([[0,0,0,Z(9)^2],[0,0,1,0],[0,1,0,0],[-Z(9)^2,0,0,0]]*Z(9)^0, f);
<4x4 plist matrix over GF(3^2)>
gap> form := HermitianFormByMatrix(gram,f);
< hermitian form >
gap> TypeOfForm(form);
1/2
gap> BaseChangeToCanonical(form);
<immutable 4x4 plist matrix over GF(3^2)>
gap> iso := IsometricCanonicalForm(form);
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
<immutable 4x4 plist matrix over GF(3^2):
[[ 0*Z(3), 0*Z(3), 0*Z(3), Z(3^2)^2 ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ Z(3^2)^6, 0*Z(3), 0*Z(3), 0*Z(3) ]
]>
Witt Index: 2
gap> Display(iso);
Hermitian form
Gram Matrix:
<immutable 4x4 plist matrix over GF(3^2):
[[ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
]>
Witt Index: 2
gap> STOP_TEST("matobj/test_forms5.tst", 10000 );
