#Behaviour of a trivial form
mat := [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]*Z(3)^0;
form := BilinearFormByMatrix(mat,GF(3));
v := Random(GF(3)^4);
[v,v]^form;
v^form;
quit;
