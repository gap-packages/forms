#A form for W(5,3)
f := GF(3);
gram := [
[0,0,0,1,0,0], 
[0,0,0,0,1,0],
[0,0,0,0,0,1],
[-1,0,0,0,0,0],
[0,-1,0,0,0,0],
[0,0,-1,0,0,0]] * One(f);;
form := BilinearFormByMatrix( gram, f );
IsSymplecticForm( form );
Display( form );
b := BaseChangeToCanonical( form );
Display( b );
Display( b * gram * TransposedMat(b) );
quit;
