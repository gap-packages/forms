gap> START_TEST("Forms: test_conic.tst");
gap> gf := GF(8);
GF(2^3)
gap> vec := gf^3;
( GF(2^3)^3 )
gap> r := PolynomialRing( gf, 3);
GF(2^3)[x_1,x_2,x_3]
gap> poly := r.1^2 + r.2 * r.3;
x_1^2+x_2*x_3
gap> form := QuadraticFormByPolynomial( poly, r );
< quadratic form >
gap> Display( form );
Quadratic form
Gram Matrix:
 1 . .
 . . 1
 . . .
Polynomial: x_1^2+x_2*x_3

gap> IsDegenerateForm( form );
#I  Testing degeneracy of the *associated bilinear form*
true
gap> IsSingularForm( form );
false
gap> WittIndex( form );
1
gap> IsParabolicForm( form );
true
gap> r := RadicalOfForm( form );;
gap> Dimension(r);
0
gap> canonical := IsometricCanonicalForm( form );
< parabolic quadratic form >
gap> form = canonical;
true
gap> go := GO(3,8);
GO(0,3,8)
gap> mat := InvariantQuadraticForm( go )!.matrix;
[ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), Z(2)^0, 0*Z(2) ] ]
gap> gapform := QuadraticFormByMatrix( mat, GF(8) );
< quadratic form >
gap> b := BaseChangeToCanonical( gapform );
[ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), Z(2)^0, 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0 ] ]
gap> hom := BaseChangeHomomorphism( b, GF(8) );
^[ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), Z(2)^0, 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0 ] ]
gap> newgo := Image(hom, go);
Group(
[ 
  [ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), Z(2^3), 0*Z(2) ], 
      [ 0*Z(2), 0*Z(2), Z(2^3)^6 ] ], 
  [ [ Z(2)^0, 0*Z(2), 0*Z(2) ], [ Z(2)^0, Z(2)^0, Z(2)^0 ], 
      [ 0*Z(2), Z(2)^0, 0*Z(2) ] ] ])
gap> conic := Filtered(vec, x -> IsZero( x^form ));;
gap> Size(conic);
64
gap> orbs := Orbits(newgo, conic, OnRight);;
gap> List(orbs,Size);
[ 1, 63 ]
gap> STOP_TEST("test_conic.tst", 10000 );
