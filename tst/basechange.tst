gap> START_TEST("basechange.tst");

#
gap> TestBaseChangeOrthogonalBilinear := function(dim,q)
>        local F, mat, bc, m;
>        F := GF(q);
>        # produce a symmetric matrix
>        mat:=RandomInvertibleMat(dim,F);
>        mat:=mat+TransposedMat(mat);
>        # compute base change
>        bc:=BaseChangeOrthogonalBilinear(mat, F);
>        # verify the base change produces a matrix with the right properties:
>        m:=bc[1]*mat*TransposedMat(bc[1]);
>        Assert(0, m = TransposedMat(m)); # symmetric
>        # TODO: verify it is canonical
>    end;;

#
gap> for dim in [2,10,20,80] do
>      for q in [3, 5, 7, 9, 25, 27, 17^2] do
>        for i in [0..3] do
>          TestBaseChangeOrthogonalBilinear(dim+i, q);
>        od;
>      od;
>    od;

#
gap> TestBaseChangeOrthogonalQuadratic := function(dim,q)
>        local F, mat, i, j, bc, m;
>        F := GF(q);
>        # produce a symmetric matrix
>        mat:=RandomMat(dim,dim,F);
>        for i in [2..dim] do
>          for j in [1..i-1] do
>            mat[i,j] := Zero(F);
>          od;
>        od;
>        #mat:=mat+TransposedMat(mat);
>        # compute base change
>        bc:=BaseChangeOrthogonalQuadratic(mat, F);
>        # verify the base change produces a matrix with the right properties:
>        m:=bc[1]*mat*TransposedMat(bc[1]);
>        #Assert(0, m=TransposedMat(m)); # symmetric
>        # TODO: verify it is canonical
>    end;;
gap> for dim in [3..10] do
>      for q in [2, 4, 8, 16, 2^9] do
>        for i in [1..5] do
>          TestBaseChangeOrthogonalQuadratic(dim, q);
>        od;
>      od;
>    od;

#
gap> TestBaseChangeHermitian := function(dim,q)
>        local F, mat, bc, m;
>        F := GF(q^2);
>        # produce a symmetric matrix
>        mat:=RandomInvertibleMat(dim,F);
>        mat:=mat+Forms_HERM_CONJ(mat, dim, q);
>        # compute base change
>        bc:=BaseChangeHermitian(mat, F);
>        # verify the base change produces a matrix with the right properties:
>        m:=bc[1]*mat*Forms_HERM_CONJ(bc[1], dim, q);
>        Assert(0, m = Forms_HERM_CONJ(m, dim, q));
>        # TODO: verify it is canonical
>    end;;

#
gap> for dim in [2..10] do
>      for q in [2, 3, 4, 5, 7, 9, 16, 25, 27, 17^2] do
>        for i in [1..5] do
>          TestBaseChangeHermitian(dim, q);
>        od;
>      od;
>    od;

#
gap> STOP_TEST("basechange.tst", 0);
