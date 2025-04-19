gap> START_TEST("Forms: bg_th_ex3.tst");
gap> mat := [[0,0,-2],[0,0,1],[2,-1,0]]*Z(7)^0;
[ [ 0*Z(7), 0*Z(7), Z(7)^5 ], [ 0*Z(7), 0*Z(7), Z(7)^0 ], 
  [ Z(7)^2, Z(7)^3, 0*Z(7) ] ]
gap> form := BilinearFormByMatrix(mat,GF(7));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
 . . 5
 . . 1
 2 6 .
gap> IsSymmetricForm(form);
false
gap> IsAlternatingForm(form);
true
gap> r := RadicalOfForm(form);
<vector space of dimension 1 over GF(7)>
gap> Dimension(r);
1
gap> STOP_TEST("bg_th_ex3.tst", 10000 );
