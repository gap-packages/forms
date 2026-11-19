gap> START_TEST("Forms: matobj/bg_th_ex5.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> ToVecObj := {v, F} -> Vector(F, v);;
gap> mat := ToMatObj([[0*Z(5),0*Z(5),0*Z(25),Z(25)^3],[0*Z(5),0*Z(5),Z(25)^3,0*Z(25)],
>         [0*Z(5),-Z(25)^3,0*Z(5),0*Z(5)],[-Z(25)^3,0*Z(5),0*Z(25),0*Z(25)]], GF(25));
<4x4 plist matrix over GF(5^2)>
gap> form := HermitianFormByMatrix(mat,GF(25));
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
<immutable 4x4 plist matrix over GF(5^2):
[[ 0*Z(5), 0*Z(5), 0*Z(5), Z(5^2)^3 ]
 [ 0*Z(5), 0*Z(5), Z(5^2)^3, 0*Z(5) ]
 [ 0*Z(5), Z(5^2)^15, 0*Z(5), 0*Z(5) ]
 [ Z(5^2)^15, 0*Z(5), 0*Z(5), 0*Z(5) ]
]>
gap> WittIndex(form);
2
gap> form2 := BilinearFormByMatrix(mat,GF(25));
< bilinear form >
gap> Display(form2);
Bilinear form
Gram Matrix:
<immutable 4x4 plist matrix over GF(5^2):
[[ 0*Z(5), 0*Z(5), 0*Z(5), Z(5^2)^3 ]
 [ 0*Z(5), 0*Z(5), Z(5^2)^3, 0*Z(5) ]
 [ 0*Z(5), Z(5^2)^15, 0*Z(5), 0*Z(5) ]
 [ Z(5^2)^15, 0*Z(5), 0*Z(5), 0*Z(5) ]
]>
gap> IsAlternatingForm(form2);
true
gap> Display(IsometricCanonicalForm(form));
Hermitian form
Gram Matrix:
<immutable 4x4 plist matrix over GF(5^2):
[[ Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5) ]
 [ 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5) ]
 [ 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5) ]
 [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0 ]
]>
Witt Index: 2
gap> Display(IsometricCanonicalForm(form2));
Bilinear form
Gram Matrix:
<immutable 4x4 plist matrix over GF(5^2):
[[ 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5) ]
 [ Z(5)^2, 0*Z(5), 0*Z(5), 0*Z(5) ]
 [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0 ]
 [ 0*Z(5), 0*Z(5), Z(5)^2, 0*Z(5) ]
]>
Witt Index: 2
gap> V := GF(25)^4;
( GF(5^2)^4 )
gap> u := ToVecObj([Z(5)^0,Z(5^2)^11,Z(5)^3,Z(5^2)^13 ], GF(25));
[ Z(5)^0, Z(5^2)^11, Z(5)^3, Z(5^2)^13 ]
gap> [u,u]^form;
0*Z(5)
gap> v := ToVecObj([Z(5)^0,Z(5^2)^5,Z(5^2),Z(5^2)^13 ], GF(25));
[ Z(5)^0, Z(5^2)^5, Z(5^2), Z(5^2)^13 ]
gap> [v,v]^form;                                     
0*Z(5)
gap> [u,v]^form;
Z(5^2)^7
gap> ([v,u]^form)^5;
Z(5^2)^7
gap> w := ToVecObj([Z(5^2)^21,Z(5^2)^19,Z(5^2)^4,Z(5)^3 ], GF(25));
[ Z(5^2)^21, Z(5^2)^19, Z(5^2)^4, Z(5)^3 ]
gap> [w,w]^form;
Z(5)
gap> v := ToVecObj([Z(5)^0,Z(5^2)^10,Z(5^2)^15,Z(5^2)^3 ], GF(25));
[ Z(5)^0, Z(5^2)^10, Z(5^2)^15, Z(5^2)^3 ]
gap> u := ToVecObj([Z(5)^3,Z(5^2)^9,Z(5^2)^4,Z(5^2)^16 ], GF(25));
[ Z(5)^3, Z(5^2)^9, Z(5^2)^4, Z(5^2)^16 ]
gap> w := ToVecObj([Z(5)^2,Z(5^2)^9,Z(5^2)^23,Z(5^2)^11 ], GF(25));
[ Z(5)^2, Z(5^2)^9, Z(5^2)^23, Z(5^2)^11 ]
gap> [u,v]^form;
0*Z(5)
gap> [u,w]^form;
0*Z(5)
gap> [v,w]^form;
0*Z(5)
gap> s := Subspace(V,[v,u,w]);
<vector space over GF(5^2), with 3 generators>
gap> Dimension(s);
2
gap> WittIndex(form);
2
gap> STOP_TEST("matobj/bg_th_ex5.tst", 10000 );
