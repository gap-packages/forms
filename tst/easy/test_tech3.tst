gap> START_TEST("Forms: test_tech3.tst");
gap> g := Sp(6,3);
Sp(6,3)
gap> module := GModuleByMats(GeneratorsOfGroup(g),GF(3));
rec( IsOverFiniteField := true, dimension := 6, field := GF(3), 
  generators := [ < immutable compressed matrix 6x6 over GF(3) >, 
      < immutable compressed matrix 6x6 over GF(3) > ], isMTXModule := true )
gap> dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(g);
rec( IsAbsolutelyIrreducible := true, IsIrreducible := true, 
  IsOverFiniteField := true, dimension := 6, field := GF(3), 
  generators := [ < immutable compressed matrix 6x6 over GF(3) >, 
      < immutable compressed matrix 6x6 over GF(3) > ], isMTXModule := true, 
  smashMeataxe := 
    rec( 
      algebraElement := 
        [ [ [ 1, 2 ], [ 1, 2 ] ], [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ] ], 
      algebraElementMatrix := < immutable compressed matrix 6x6 over GF(3) >, 
      characteristicPolynomial := x_1^6+x_1^5-x_1^4+x_1^3-x_1^2+x_1+Z(3)^0, 
      charpolFactors := x_1-Z(3)^0, degreeFieldExt := 1, ndimFlag := 1, 
      nullspaceVector := [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3)^0 ] ) )
gap> ClassicalForms_InvariantFormDual(module,dmodule);
[ "symplectic", < immutable compressed matrix 6x6 over GF(3) >, 
  [ Z(3)^0, Z(3)^0 ] ]
gap> g := GO(3,3);
GO(0,3,3)
gap> module := GModuleByMats(GeneratorsOfGroup(g),GF(3));
rec( IsOverFiniteField := true, dimension := 3, field := GF(3), 
  generators := 
    [ 
      [ [ Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), Z(3), 0*Z(3) ], 
          [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], 
      [ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ Z(3)^0, Z(3), Z(3) ], 
          [ 0*Z(3), Z(3), Z(3)^0 ] ] ], isMTXModule := true )
gap> dmodule := ClassicalForms_GeneratorsWithoutScalarsDual(g);
rec( IsAbsolutelyIrreducible := true, IsIrreducible := true, 
  IsOverFiniteField := true, dimension := 3, field := GF(3), 
  generators := 
    [ 
      [ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ Z(3)^0, Z(3), Z(3)^0 ], 
          [ 0*Z(3), Z(3)^0, Z(3)^0 ] ], 
      [ [ Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), Z(3), 0*Z(3) ], 
          [ 0*Z(3), 0*Z(3), Z(3) ] ], 
      [ [ Z(3)^0, 0*Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0, 0*Z(3) ], 
          [ 0*Z(3), 0*Z(3), Z(3) ] ], 
      [ [ 0*Z(3), Z(3), 0*Z(3) ], [ Z(3), Z(3)^0, Z(3)^0 ], 
          [ 0*Z(3), Z(3), Z(3)^0 ] ] ], isMTXModule := true, 
  smashMeataxe := 
    rec( 
      algebraElement := 
        [ [ [ 1, 2 ], [ 3, 5 ] ], 
          [ Z(3), Z(3), Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ] ], 
      algebraElementMatrix := 
        [ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ Z(3)^0, Z(3), Z(3)^0 ], 
          [ 0*Z(3), 0*Z(3), Z(3) ] ], 
      characteristicPolynomial := x_1^3-x_1^2-Z(3)^0, 
      charpolFactors := x_1+Z(3)^0, degreeFieldExt := 1, ndimFlag := 1, 
      nullspaceVector := [ 0*Z(3), 0*Z(3), Z(3)^0 ] ) )
gap> ClassicalForms_InvariantFormDual(module,dmodule);
[ "orthogonalcircle", 
  [ [ 0*Z(3), Z(3), 0*Z(3) ], [ Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), Z(3)^0 ] ], [ Z(3)^0, Z(3)^0 ], 
  [ [ 0*Z(3), Z(3), 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3) ], 
      [ 0*Z(3), 0*Z(3), Z(3) ] ] ]
gap> STOP_TEST("test_tech3.tst", 10000 );
