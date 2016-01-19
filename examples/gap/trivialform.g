#Constructing a trivial form
mat := [[0,0,0],[0,0,0],[0,0,0]]*Z(7)^0;
form1 := BilinearFormByMatrix(mat,GF(7));
form2 := QuadraticFormByMatrix(mat,GF(7));
form1 = form2;
IsQuadraticForm(form1);
IsSesquilinearForm(form1);
mat := [[0,0],[0,0]]*Z(4)^0;
form3 := BilinearFormByMatrix(mat,GF(4));
form3 = form1;
quit;
