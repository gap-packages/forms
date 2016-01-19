#computing subspace orthogonal to given vector
mat := [[Z(8),0,0,0],[0,0,Z(8)^4,0],[0,0,0,1],[0,0,0,0]]*Z(8)^0;;
form := QuadraticFormByMatrix(mat,GF(8));
u := [ Z(2^3)^4, Z(2^3)^4, Z(2)^0, Z(2^3)^3 ];
EvaluateForm( form, u );
u^form;
gram := [[0,0,0,0,0,2],[0,0,0,0,2,0],[0,0,0,1,0,0],
              [0,0,1,0,0,0],[0,2,0,0,0,0],[2,0,0,0,0,0]]*Z(3)^0;;
form := BilinearFormByMatrix(gram,GF(3));
u := [ [ Z(3)^0, 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), Z(3)^0 ], 
  [ 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3), 0*Z(3) ] ];;
v := [ [ Z(3)^0, 0*Z(3), Z(3)^0, Z(3), 0*Z(3), Z(3) ], 
  [ 0*Z(3), Z(3)^0, 0*Z(3), Z(3), Z(3), Z(3) ] ];;
EvaluateForm( form, u, v);
[u,v]^form;
quit;
