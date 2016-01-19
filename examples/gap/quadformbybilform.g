#Constructing form: QuadraticFormByBilinearForm
mat := [ [ Z(3^2)^7, Z(3)^0, Z(3^2)^2, 0*Z(3), Z(3^2)^5 ], 
   [ Z(3)^0, Z(3^2)^7, Z(3^2)^6, Z(3^2)^5, Z(3^2)^2 ], 
   [ Z(3^2)^2, Z(3^2)^6, Z(3^2)^7, Z(3^2)^2, Z(3^2)^2 ], 
   [ 0*Z(3), Z(3^2)^5, Z(3^2)^2, Z(3^2)^6, Z(3^2)^7 ], 
   [ Z(3^2)^5, Z(3^2)^2, Z(3^2)^2, Z(3^2)^7, Z(3) ] ];
form := BilinearFormByMatrix(mat,GF(9));
Q := QuadraticFormByBilinearForm(form);
Display(form);
Display(Q);
Set(List(GF(9)^5),x->[x,x]^form=x^Q);
PolynomialOfForm(form);
PolynomialOfForm(Q);
quit;
