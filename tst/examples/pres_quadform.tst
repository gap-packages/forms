gap> START_TEST("Forms: pres_quadform.tst");
gap> group := GO(-1,4,8);
GO(-1,4,8)
gap> pres_forms := PreservedQuadraticForms(group);
[ < quadratic form > ]
gap> group := GO(1,6,9);
GO(+1,6,9)
gap> pres_forms := PreservedQuadraticForms(group);
[ < quadratic form > ]
gap> STOP_TEST("pres_quadform.tst", 10000 );
