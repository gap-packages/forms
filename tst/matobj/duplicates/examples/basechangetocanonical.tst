gap> START_TEST("Forms: matobj/basechangetocanonical.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> gf := GF(3);
GF(3)
gap> gram := ToMatObj([
> [0,0,0,1,0,0], 
> [0,0,0,0,1,0],
> [0,0,0,0,0,1],
> [-1,0,0,0,0,0],
> [0,-1,0,0,0,0],
> [0,0,-1,0,0,0]] * One(gf), gf);;
gap> form := BilinearFormByMatrix( gram, gf );
< bilinear form >
gap> b := BaseChangeToCanonical( form );;
gap> Display( b * gram * TransposedMat(b) );
<6x6 plist matrix over GF(3):
[[ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3) ]
]>
gap> STOP_TEST("matobj/basechangetocanonical.tst", 10000 );
