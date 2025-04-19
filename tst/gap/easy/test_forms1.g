#test_forms1.g
f := GF(3);;
gram := [[0,0,0,0,0,2],[0,0,0,0,2,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,2,0,0,0,0],[2,0,0,0,0,0]]*Z(3)^0;
form := BilinearFormByMatrix(gram,f);
Display(BaseChangeToCanonical(form));
Display(form);
TypeOfForm(form);
quit;
