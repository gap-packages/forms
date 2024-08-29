gap> START_TEST("Forms: test_tech2.tst");
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
gap> ScalarsOfPreservedForm(grp,forms[1]);
[ Z(5)^0, Z(5)^0 ]
gap> ScalarsOfPreservedForm(go,forms[1]);
false
gap> quad := QuadraticFormByBilinearForm(forms[1]);
< quadratic form >
gap> ScalarsOfPreservedForm(grp,quad);
[ Z(5)^0, Z(5)^0 ]
gap> ScalarsOfPreservedForm(go,quad);
false
gap> STOP_TEST("test_tech2.tst", 10000 );
