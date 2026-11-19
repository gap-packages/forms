gap> START_TEST("Forms: matobj/form_eval_with_list_vectors.tst");
gap> F := GF(25);;
gap> ToMatObj := {m} -> Matrix(IsPlistMatrixRep, F, m);;
gap> ToVecObj := {v} -> Vector(IsPlistVectorRep, F, v);;
gap> mat := ToMatObj([[0*Z(5),0*Z(5),0*Z(25),Z(25)^3],[0*Z(5),0*Z(5),Z(25)^3,0*Z(25)],
>         [0*Z(5),-Z(25)^3,0*Z(5),0*Z(5)],[-Z(25)^3,0*Z(5),0*Z(25),0*Z(25)]]);;
gap> form := HermitianFormByMatrix(mat, F);
< hermitian form >
gap> u := ([Z(5)^0,Z(5^2)^11,Z(5)^3,Z(5^2)^13]);;
gap> v := ([Z(5)^0,Z(5^2)^5,Z(5^2),Z(5^2)^13]);;
gap> w := ([Z(5^2)^21,Z(5^2)^19,Z(5^2)^4,Z(5)^3]);;
gap> IsVectorObj(u) and IsVectorObj(v) and IsVectorObj(w);
false
gap> [u,u]^form; # TODO: This should be a question if this is even allowed when the gram matrix is a matrix object??. For now i have decided it should.
0*Z(5)
gap> [v,v]^form;
0*Z(5)
gap> [u,v]^form;
Z(5^2)^7
gap> ([v,u]^form)^5;
Z(5^2)^7
gap> [w,w]^form;
Z(5)
gap> v := ([Z(5)^0,Z(5^2)^10,Z(5^2)^15,Z(5^2)^3]);;
gap> u := ([Z(5)^3,Z(5^2)^9,Z(5^2)^4,Z(5^2)^16]);;
gap> w := ([Z(5)^2,Z(5^2)^9,Z(5^2)^23,Z(5^2)^11]);;
gap> [u,v]^form;
0*Z(5)
gap> [u,w]^form;
0*Z(5)
gap> [v,w]^form;
0*Z(5)
gap> STOP_TEST("matobj/form_eval_with_list_vectors.tst", 10000);
