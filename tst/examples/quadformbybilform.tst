gap> START_TEST("Forms: quadformbybilform.tst");
gap> mat := [ [ Z(3^2)^7, Z(3)^0, Z(3^2)^2, 0*Z(3), Z(3^2)^5 ], 
>    [ Z(3)^0, Z(3^2)^7, Z(3^2)^6, Z(3^2)^5, Z(3^2)^2 ], 
>    [ Z(3^2)^2, Z(3^2)^6, Z(3^2)^7, Z(3^2)^2, Z(3^2)^2 ], 
>    [ 0*Z(3), Z(3^2)^5, Z(3^2)^2, Z(3^2)^6, Z(3^2)^7 ], 
>    [ Z(3^2)^5, Z(3^2)^2, Z(3^2)^2, Z(3^2)^7, Z(3) ] ];
[ [ Z(3^2)^7, Z(3)^0, Z(3^2)^2, 0*Z(3), Z(3^2)^5 ], 
  [ Z(3)^0, Z(3^2)^7, Z(3^2)^6, Z(3^2)^5, Z(3^2)^2 ], 
  [ Z(3^2)^2, Z(3^2)^6, Z(3^2)^7, Z(3^2)^2, Z(3^2)^2 ], 
  [ 0*Z(3), Z(3^2)^5, Z(3^2)^2, Z(3^2)^6, Z(3^2)^7 ], 
  [ Z(3^2)^5, Z(3^2)^2, Z(3^2)^2, Z(3^2)^7, Z(3) ] ]
gap> form := BilinearFormByMatrix(mat,GF(9));
< bilinear form >
gap> Q := QuadraticFormByBilinearForm(form);
< quadratic form >
gap> Display(form);
Bilinear form
Gram Matrix:
z = Z(9)
 z^7   1 z^2   . z^5
   1 z^7 z^6 z^5 z^2
 z^2 z^6 z^7 z^2 z^2
   . z^5 z^2 z^6 z^7
 z^5 z^2 z^2 z^7   2
gap> Display(Q);
Quadratic form
Gram Matrix:
z = Z(9)
 z^7   2 z^6   . z^1
   . z^7 z^2 z^1 z^6
   .   . z^7 z^6 z^6
   .   .   . z^6 z^3
   .   .   .   .   2
gap> Set(List(GF(9)^5),x->[x,x]^form=x^Q);
[ true ]
gap> PolynomialOfForm(form);
Z(3^2)^7*x_1^2-x_1*x_2+Z(3^2)^6*x_1*x_3+Z(3^2)*x_1*x_5+Z(3^2)^7*x_2^2+Z(3^2)^2\
*x_2*x_3+Z(3^2)*x_2*x_4+Z(3^2)^6*x_2*x_5+Z(3^2)^7*x_3^2+Z(3^2)^6*x_3*x_4+Z(3^2\
)^6*x_3*x_5+Z(3^2)^6*x_4^2+Z(3^2)^3*x_4*x_5-x_5^2
gap> PolynomialOfForm(Q);
Z(3^2)^7*x_1^2-x_1*x_2+Z(3^2)^6*x_1*x_3+Z(3^2)*x_1*x_5+Z(3^2)^7*x_2^2+Z(3^2)^2\
*x_2*x_3+Z(3^2)*x_2*x_4+Z(3^2)^6*x_2*x_5+Z(3^2)^7*x_3^2+Z(3^2)^6*x_3*x_4+Z(3^2\
)^6*x_3*x_5+Z(3^2)^6*x_4^2+Z(3^2)^3*x_4*x_5-x_5^2
gap> STOP_TEST("quadformbybilform.tst", 10000 );
