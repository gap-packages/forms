#Background theory: example 4
mat := [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,0,0,1],
        [0,0,0,0,1,0],[0,0,0,1,0,0],[0,0,1,0,0,0]]*Z(16)^0;
form := BilinearFormByMatrix(mat,GF(16));
Display(form);
IsSymmetricForm(form);
IsAlternatingForm(form);
IsDegenerateForm(form);
WittIndex(form);
quit;
