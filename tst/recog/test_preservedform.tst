gap> START_TEST("Forms: test_preservedform.tst");
gap> go := GO(5, 5);
GO(0,5,5)
gap> x := 
> [ [ Z(5)^0, Z(5)^3, 0*Z(5), Z(5)^3, Z(5)^3 ], 
>   [ Z(5)^2, Z(5)^3, 0*Z(5), Z(5)^2, Z(5) ], 
>   [ Z(5)^2, Z(5)^2, Z(5)^0, Z(5), Z(5)^3 ],
>   [ Z(5)^0, Z(5)^3, Z(5), Z(5)^0, Z(5)^3 ], 
>   [ Z(5)^3, 0*Z(5), Z(5)^0, 0*Z(5), Z(5) ] 
>  ];;
gap> grp := go^x;
<matrix group of size 18720000 with 2 generators>
gap> forms := PreservedSesquilinearForms( grp );;
gap> TestPreservedSesquilinearForms(grp,forms);
true
gap> STOP_TEST("test_preservedform.tst", 10000 );
