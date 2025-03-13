# don't invoke START_TEST, which resets the random number generators,
# so that we can re-read this .tst file multiple times to test more cases.
#gap> START_TEST("basechange.tst");

#
gap> MakeCanonicalOrthogonalBilinear := function(dim,F,rank,w)
>        local A, e, i, offset;
>        A:=NullMat(dim,dim,F);
>        e:=One(F);
>        if w = 0 then
>          A[1,1] := e;
>          if Size(F) mod 4 = 3 then
>            A[2,2] := e;
>          else
>            A[2,2] := PrimitiveRoot(F);
>          fi;
>        elif w = 1 then
>          A[1,1] := e;
>        else
>          Assert(0, w = 2);
>        fi;
>        offset := 2-w;
>        for i in [1..QuoInt(rank - offset,2)] do
>          A[offset+2*i,offset+2*i-1] := e/2;
>          A[offset+2*i-1,offset+2*i] := e/2;
>        od;
>        return A;
>    end;;
gap> TestBaseChangeOrthogonalBilinear := function(dim,q)
>        local F, mat, bc, m, c;
>        F := GF(q);
>        # produce a symmetric matrix
>        mat:=RandomMat(dim,dim,F);
>        mat:=mat+TransposedMat(mat);
>        if IsZero(mat) then return; fi;
>        # compute base change
>        bc:=BaseChangeOrthogonalBilinear(mat, F);
>        Assert(0, dim = RankMat(bc[1]));
>        Assert(0, bc[2]+1 = RankMat(mat));  # FIXME: why +1?
>        # verify the base change really yields the canonical form
>        m:=bc[1]*mat*TransposedMat(bc[1]);
>        c := MakeCanonicalOrthogonalBilinear(dim,F,bc[2]+1,bc[3]);
>        # HACK: 'normalize' both matrices
>        Assert(0, m/First(m[1],x->not IsZero(x)) = c/First(c[1],x->not IsZero(x)));
>    end;;

#
gap> for dim in [1,2,10,20] do
>      for q in [3, 5, 7, 9, 25, 27, 17^2] do
>        for i in [0..4] do
>          TestBaseChangeOrthogonalBilinear(dim+i, q);
>        od;
>      od;
>    od;

#
gap> MakeCanonicalOrthogonalQuadratic := function(dim,F,rank,w)
>        local A, e, i, offset;
>        A:=NullMat(dim,dim,F);
>        e:=One(F);
>        if w = 0 then
>          A[1,1] := e;
>          A[1,2] := e;
>          A[2,2] := Forms_C1(F, DegreeOverPrimeField(F));
>        elif w = 1 then
>          A[1,1] := e;
>        else
>          Assert(0, w = 2);
>        fi;
>        offset := 2-w;
>        for i in [1..QuoInt(rank - offset,2)] do
>          A[offset+2*i-1,offset+2*i] := e;
>        od;
>        return A;
>    end;;
gap> TestBaseChangeOrthogonalQuadratic := function(dim,q)
>        local F, mat, i, j, bc, m, c;
>        F := GF(q);
>        # produce upper triangular matrix
>        mat:=NullMat(dim,dim,F);
>        for i in [1..dim] do
>          for j in [i..dim] do
>            mat[i,j] := Random(F);
>          od;
>        od;
>        if IsZero(mat) then return; fi;
>        # compute base change
>        bc:=BaseChangeOrthogonalQuadratic(mat, F);
>        Assert(0, dim = RankMat(bc[1]));
>        #Assert(0, bc[2] = RankMat(mat));  # FIXME: what does bc[2] really mean?
>        # verify the base change really yields the canonical form
>        m:=bc[1]*mat*TransposedMat(bc[1]);
>        c:=MakeCanonicalOrthogonalQuadratic(dim,F,bc[2]+1,bc[3]);
>        Assert(0, Forms_RESET(m, dim, q) = c);
>    end;;

#
gap> for dim in [1,2,10,20] do
>      for q in [2, 4, 8, 16, 2^9] do
>        for i in [0..4] do
>          TestBaseChangeOrthogonalQuadratic(dim+1, q);
>        od;
>      od;
>    od;

#
gap> MakeCanonicalHermitian := function(dim,F,rank)
>        local A, e, i;
>        A:=NullMat(dim,dim,F);
>        e:=One(F);
>        for i in [1..rank] do
>          A[i,i] := e;
>        od;
>        return A;
>    end;;
gap> TestBaseChangeHermitian := function(dim,q)
>        local F, mat, bc, m, c;
>        F := GF(q^2);
>        # produce a symmetric matrix
>        mat:=RandomInvertibleMat(dim,F);
>        mat:=mat+Forms_HERM_CONJ(mat, q);
>        if IsZero(mat) then return; fi;
>        # compute base change
>        bc:=BaseChangeHermitian(mat, F);
>        Assert(0, dim = RankMat(bc[1]));
>        # verify the base change really yields the canonical form
>        m:=bc[1]*mat*Forms_HERM_CONJ(bc[1], q);
>        c:=MakeCanonicalHermitian(dim,F,RankMat(mat));
>        Assert(0, m = c);
>    end;;

#
gap> for dim in [1,2,10,20] do
>      for q in [2, 3, 4, 5, 7, 9, 16, 25, 27, 17^2] do
>        for i in [0..4] do
>          TestBaseChangeHermitian(dim+i, q);
>        od;
>      od;
>    od;

#
gap> MakeCanonicalSymplectic := function(dim,F,rank)
>        local A, e, i;
>        A:=NullMat(dim,dim,F);
>        e:=One(F);
>        for i in [1..QuoInt(rank,2)] do
>          A[2*i,2*i-1] := -e;
>          A[2*i-1,2*i] := e;
>        od;
>        return A;
>    end;;
gap> TestBaseChangeSymplectic := function(dim,q)
>        local F, mat, bc, m, c;
>        F := GF(q);
>        # produce an alternating matrix
>        mat:=RandomMat(dim,dim,F);
>        mat:=mat-TransposedMat(mat);
>        if IsZero(mat) then return; fi;
>        # compute base change
>        bc:=BaseChangeSymplectic(mat, F);
>        Assert(0, dim = RankMat(bc[1]));
>        # verify the base change really yields the canonical form
>        m:=bc[1]*mat*TransposedMat(bc[1]);
>        c:=MakeCanonicalSymplectic(dim,F,RankMat(mat));
>        Assert(0, m = c);
>    end;;

#
gap> for dim in [2,10,20] do
>      for q in [2, 3, 4, 5, 7, 8, 9, 25, 27, 17^2] do
>        for i in [0..4] do
>          TestBaseChangeSymplectic(dim, q);
>        od;
>      od;
>    od;

#
#gap> STOP_TEST("basechange.tst", 0);
