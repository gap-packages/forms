#test_forms2.g
f := GF(7);
gram := [[-3,0,0,0,0,0],[0,0,3,0,0,0],[0,3,0,0,0,0],[0,0,0,0,0,-1/2],[0,0,0,0,1,0],[0,0,0,-1/2,0,0]]*Z(7)^0;
form := BilinearFormByMatrix(gram,f);
IsEllipticForm(form);
TypeOfForm(form);
Display(form);
quit;
