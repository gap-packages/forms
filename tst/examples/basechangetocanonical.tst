gap> START_TEST("Forms: basechangetocanonical.tst");
gap> #morphisms: BaseChangeToCanonical
gap> gf := GF(3);
GF(3)
gap> gram := [
> [0,0,0,1,0,0], 
> [0,0,0,0,1,0],
> [0,0,0,0,0,1],
> [-1,0,0,0,0,0],
> [0,-1,0,0,0,0],
> [0,0,-1,0,0,0]] * One(gf);;
gap> form := BilinearFormByMatrix( gram, gf );
< bilinear form >
gap> b := BaseChangeToCanonical( form );;
gap> Display( b * gram * TransposedMat(b) );
 . 1 . . . .
 2 . . . . .
 . . . 1 . .
 . . 2 . . .
 . . . . . 1
 . . . . 2 .
gap> STOP_TEST("basechangetocanonical.tst", 10000 );
