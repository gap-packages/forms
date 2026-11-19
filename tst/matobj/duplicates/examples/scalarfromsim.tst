gap> START_TEST("Forms: matobj/scalarfromsim.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> gram := ToMatObj([ [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ] ], GF(3));;
gap> form := BilinearFormByMatrix( gram, GF(3) );
< bilinear form >
gap> m := ToMatObj([ [ Z(3)^0, Z(3)^0, Z(3), 0*Z(3), Z(3)^0, Z(3) ], 
>   [ Z(3), Z(3), Z(3)^0, 0*Z(3), Z(3)^0, Z(3) ], 
>   [ 0*Z(3), Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), Z(3), Z(3)^0, Z(3), Z(3), Z(3) ], 
>   [ Z(3)^0, Z(3)^0, Z(3), Z(3), Z(3)^0, Z(3)^0 ], 
>   [ Z(3)^0, 0*Z(3), Z(3), Z(3)^0, Z(3), Z(3) ] ], GF(3));;
gap> ScalarOfSimilarity( m, form );
Z(3)
gap> STOP_TEST("matobj/scalarfromsim.tst", 10000 );
