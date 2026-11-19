#@local f, mat, form
gap> START_TEST("Forms: matobj/gram_matrix_repr.tst");

# GramMatrix must return a matrix in the same representation as the
# matrix the form was created from.

# lists in lists (forms currently compresses these, so only the
# list-of-lists interface is required here)
gap> f := GF(5);
GF(5)
gap> mat := [[1,0],[0,1]] * One(f);
[ [ Z(5)^0, 0*Z(5) ], [ 0*Z(5), Z(5)^0 ] ]
gap> form := BilinearFormByMatrix(mat, f);
< bilinear form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat, f);
< quadratic form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := HermitianFormByMatrix(mat, GF(25));
< hermitian form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := BilinearFormByMatrix([[0,0],[0,0]] * One(f), f);
< trivial form >
gap> IsMatrix(GramMatrix(form));
true

# compressed matrices over GF(q), 2 < q <= 256
gap> mat := Matrix(Is8BitMatrixRep, f, [[1,0],[0,1]] * One(f));
[ [ Z(5)^0, 0*Z(5) ], [ 0*Z(5), Z(5)^0 ] ]
gap> form := BilinearFormByMatrix(mat, f);
< bilinear form >
gap> Is8BitMatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> Is8BitMatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat, f);
< quadratic form >
gap> Is8BitMatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> Is8BitMatrixRep(GramMatrix(form));
true
gap> form := HermitianFormByMatrix(Matrix(Is8BitMatrixRep, GF(25), [[1,0],[0,1]] * One(f)), GF(25));
< hermitian form >
gap> Is8BitMatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(ZeroMatrix(2, 2, mat), f);
< trivial form >
gap> Is8BitMatrixRep(GramMatrix(form));
true

# compressed matrices over GF(2)
gap> f := GF(2);
GF(2)
gap> mat := Matrix(IsGF2MatrixRep, f, [[1,0],[0,1]] * One(f));
<a 2x2 matrix over GF2>
gap> form := BilinearFormByMatrix(mat, f);
< bilinear form >
gap> IsGF2MatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> IsGF2MatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat, f);
< quadratic form >
gap> IsGF2MatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> IsGF2MatrixRep(GramMatrix(form));
true
gap> form := HermitianFormByMatrix(Matrix(Is8BitMatrixRep, GF(4), [[1,0],[0,1]] * One(f)), GF(4));
< hermitian form >
gap> Is8BitMatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(ZeroMatrix(2, 2, mat), f);
< trivial form >
gap> IsGF2MatrixRep(GramMatrix(form));
true

# matrix objects
gap> f := GF(5);
GF(5)
gap> mat := Matrix(IsPlistMatrixRep, f, [[1,0],[0,1]] * One(f));
<2x2 plist matrix over GF(5)>
gap> form := BilinearFormByMatrix(mat, f);
< bilinear form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat, f);
< quadratic form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := HermitianFormByMatrix(ChangedBaseDomain(mat, GF(25)), GF(25));
< hermitian form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(ZeroMatrix(2, 2, mat), f);
< trivial form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> STOP_TEST("matobj/gram_matrix_repr.tst", 10000);
