gap> START_TEST("Forms: matobj/w53.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> f := GF(3);
GF(3)
gap> gram := ToMatObj([
> [0,0,0,1,0,0], 
> [0,0,0,0,1,0],
> [0,0,0,0,0,1],
> [-1,0,0,0,0,0],
> [0,-1,0,0,0,0],
> [0,0,-1,0,0,0]] * One(f), f);;
gap> form := BilinearFormByMatrix( gram, f );
< bilinear form >
gap> IsSymplecticForm( form );
true
gap> Display( form );
Symplectic form
Gram Matrix:
<immutable 6x6 plist matrix over GF(3):
[[ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
 [ Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
]>
gap> b := BaseChangeToCanonical( form );
<immutable 6x6 plist matrix over GF(3)>
gap> Display( b );
<immutable 6x6 plist matrix over GF(3):
[[ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ]
 [ 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
]>
gap> Display( b * gram * TransposedMat(b) );
<6x6 plist matrix over GF(3):
[[ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3) ]
]>
gap> STOP_TEST("matobj/w53.tst", 10000 );
