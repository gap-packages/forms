gap> START_TEST("Forms: matobj/istotallyisotropic.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> mat := ToMatObj([[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]]*Z(7)^0, GF(7));
<4x4 plist matrix over GF(7)>
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> sub:= ToMatObj([[Z(7)^0,0*Z(7),Z(7)^0,Z(7)],[0*Z(7),Z(7)^0,Z(7)^0,Z(7)^4]], GF(7)); # TODO: this right now also requires the basis of the given vector space to be the correct type of matrix object due to compatibility issues which in my opinion is kind of stupid. A better API might only accept actual vector space objects or something??
<2x4 plist matrix over GF(7)>
gap> IsTotallyIsotropicSubspace(form,sub);
true
gap> mat := ToMatObj(IdentityMat(6,GF(2)), GF(4));
<6x6 plist matrix over GF(2^2)>
gap> form := HermitianFormByMatrix(mat,GF(4));
< hermitian form >
gap> sub := [[Z(2)^0,0*Z(2),0*Z(2),Z(2)^0,Z(2)^0,Z(2)^0], 
>   [0*Z(2),Z(2)^0,0*Z(2),Z(2^2)^2,Z(2^2),Z(2)^0], 
>   [0*Z(2),0*Z(2),Z(2)^0,Z(2)^0,Z(2^2),Z(2^2)^2]];
[ [ Z(2)^0, 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2)^0 ], 
  [ 0*Z(2), Z(2)^0, 0*Z(2), Z(2^2)^2, Z(2^2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2^2), Z(2^2)^2 ] ]
gap> IsTotallyIsotropicSubspace(form,sub);
true
gap> STOP_TEST("matobj/istotallyisotropic.tst", 10000 );
