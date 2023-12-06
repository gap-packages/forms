gap> START_TEST("Forms: preservedform.tst");
gap> #What is the form preserved by this group?
gap> go := GO(5, 5);
GO(0,5,5)
gap> x := 
> [ [ Z(5)^0, Z(5)^3, 0*Z(5), Z(5)^3, Z(5)^3 ], 
>   [ Z(5)^2, Z(5)^3, 0*Z(5), Z(5)^2, Z(5) ], 
>   [ Z(5)^2, Z(5)^2, Z(5)^0, Z(5), Z(5)^3 ],
>   [ Z(5)^0, Z(5)^3, Z(5), Z(5)^0, Z(5)^3 ], 
>   [ Z(5)^3, 0*Z(5), Z(5)^0, 0*Z(5), Z(5) ] 
>  ];;
gap> go2 := go^x;
<matrix group of size 18720000 with 2 generators>
gap> forms := PreservedSesquilinearForms( go2 );
[ < bilinear form > ]
gap> Display( forms[1] );
Bilinear form
Gram Matrix:
 2 1 2 4 4
 1 1 1 4 4
 2 1 4 3 2
 4 4 3 1 2
 4 4 2 2 4
gap> STOP_TEST("preservedform.tst", 10000 );
