gap> START_TEST("Forms: scalarfromsim.tst");
gap> #ScalarOfSimilairty example
gap> gram := [ [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ] ];;
gap> form := BilinearFormByMatrix( gram, GF(3) );
< bilinear form >
gap> m := [ [ Z(3)^0, Z(3)^0, Z(3), 0*Z(3), Z(3)^0, Z(3) ], 
>   [ Z(3), Z(3), Z(3)^0, 0*Z(3), Z(3)^0, Z(3) ], 
>   [ 0*Z(3), Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), Z(3), Z(3)^0, Z(3), Z(3), Z(3) ], 
>   [ Z(3)^0, Z(3)^0, Z(3), Z(3), Z(3)^0, Z(3)^0 ], 
>   [ Z(3)^0, 0*Z(3), Z(3), Z(3)^0, Z(3), Z(3) ] ];;
gap> ScalarOfSimilarity( m, form );
Z(3)
gap> STOP_TEST("scalarfromsim.tst", 10000 );
