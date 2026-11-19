gap> START_TEST("Forms: matobj/quadformbybilform.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([ [ Z(3^2)^7, Z(3)^0, Z(3^2)^2, 0*Z(3), Z(3^2)^5 ], 
>    [ Z(3)^0, Z(3^2)^7, Z(3^2)^6, Z(3^2)^5, Z(3^2)^2 ], 
>    [ Z(3^2)^2, Z(3^2)^6, Z(3^2)^7, Z(3^2)^2, Z(3^2)^2 ], 
>    [ 0*Z(3), Z(3^2)^5, Z(3^2)^2, Z(3^2)^6, Z(3^2)^7 ], 
>    [ Z(3^2)^5, Z(3^2)^2, Z(3^2)^2, Z(3^2)^7, Z(3) ] ], GF(9));
<5x5 plist matrix over GF(3^2)>
gap> form := BilinearFormByMatrix(mat,GF(9));
< bilinear form >
gap> Q := QuadraticFormByBilinearForm(form);
< quadratic form >
gap> Display(form);
Bilinear form
Gram Matrix:
<immutable 5x5 plist matrix over GF(3^2):
[[ Z(3^2)^7, Z(3)^0, Z(3^2)^2, 0*Z(3), Z(3^2)^5 ]
 [ Z(3)^0, Z(3^2)^7, Z(3^2)^6, Z(3^2)^5, Z(3^2)^2 ]
 [ Z(3^2)^2, Z(3^2)^6, Z(3^2)^7, Z(3^2)^2, Z(3^2)^2 ]
 [ 0*Z(3), Z(3^2)^5, Z(3^2)^2, Z(3^2)^6, Z(3^2)^7 ]
 [ Z(3^2)^5, Z(3^2)^2, Z(3^2)^2, Z(3^2)^7, Z(3) ]
]>
gap> Display(Q);
Quadratic form
Gram Matrix:
<immutable 5x5 plist matrix over GF(3^2):
[[ Z(3^2)^7, Z(3), Z(3^2)^6, 0*Z(3), Z(3^2) ]
 [ 0*Z(3), Z(3^2)^7, Z(3^2)^2, Z(3^2), Z(3^2)^6 ]
 [ 0*Z(3), 0*Z(3), Z(3^2)^7, Z(3^2)^6, Z(3^2)^6 ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3^2)^6, Z(3^2)^3 ]
 [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3) ]
]>
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
gap> STOP_TEST("matobj/quadformbybilform.tst", 10000 );
