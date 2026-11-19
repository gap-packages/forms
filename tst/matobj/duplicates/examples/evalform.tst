gap> START_TEST("Forms: matobj/evalform.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[Z(8),0,0,0],[0,0,Z(8)^4,0],[0,0,0,1],[0,0,0,0]]*Z(8)^0, GF(8));;
gap> form := QuadraticFormByMatrix(mat,GF(8));
< quadratic form >
gap> u := [ Z(2^3)^4, Z(2^3)^4, Z(2)^0, Z(2^3)^3 ];
[ Z(2^3)^4, Z(2^3)^4, Z(2)^0, Z(2^3)^3 ]
gap> EvaluateForm( form, u );
Z(2^3)^6
gap> u^form;
Z(2^3)^6
gap> gram := ToMatObj([[0,0,0,0,0,2],[0,0,0,0,2,0],[0,0,0,1,0,0],
>               [0,0,1,0,0,0],[0,2,0,0,0,0],[2,0,0,0,0,0]]*Z(3)^0, GF(3));;
gap> form := BilinearFormByMatrix(gram,GF(3));
< bilinear form >
gap> u := [ [ Z(3)^0, 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), Z(3)^0 ], 
>   [ 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3), 0*Z(3) ] ];;
gap> v := [ [ Z(3)^0, 0*Z(3), Z(3)^0, Z(3), 0*Z(3), Z(3) ], 
>   [ 0*Z(3), Z(3)^0, 0*Z(3), Z(3), Z(3), Z(3) ] ];;
gap> # gap> EvaluateForm( form, u, v); TODO: this does not work currently (by choice) to not deal with the whole what to do if i get a matrix and a matrix object?? should it work
gap> # [ [ Z(3)^0, Z(3)^0 ], [ 0*Z(3), 0*Z(3) ] ]
gap> #[u,v]^form;
gap> # [ [ Z(3)^0, Z(3)^0 ], [ 0*Z(3), 0*Z(3) ] ]
gap> #STOP_TEST("matobj/evalform.tst", 10000 );
