#example of MatrixIntertwiners
A := [[1,1,0],[0,1,0],[0,0,2]]*Z(3)^0;
B := [[2,0,0],[0,1,1],[0,0,1]]*Z(3)^0;
C := MatrixIntertwiners([A], [A], GF(3));;
Size(C);
C := MatrixIntertwiners([A], [B], GF(3));
ForAll(C, c -> A*c = c*B);
quit;
