#morphisms: BaseChangeToCanonical
gf := GF(3);
gram := [
[0,0,0,1,0,0], 
[0,0,0,0,1,0],
[0,0,0,0,0,1],
[-1,0,0,0,0,0],
[0,-1,0,0,0,0],
[0,0,-1,0,0,0]] * One(gf);;
form := BilinearFormByMatrix( gram, gf );
b := BaseChangeToCanonical( form );;
Display( b * gram * TransposedMat(b) );
quit;
