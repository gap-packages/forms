gap> START_TEST("Formspace: Preserved Formspace SU(4,5)");
gap> G := SU(4,5);; # Some groups are defined in interesting_groups.g
gap> R := PseudoRandom(GL(4, GF(5^2)));;
gap> G := G^R;;
gap> f_space_expected_normal_d := 0;; # the dimensions of the expected formspaces 
gap> f_space_expected_unitary_d := 1;;
gap> L:=PreservedFormspace(G);;
gap> TestMatricesAreForms2(G, L); # found in custom_test_functions.g, Tests whether the matrices returned are forms preserved by G or not
[ "Ok", "Ok" ]
gap> Dimension(DefualtFieldOfMatrixGroup(G), Vectorspace(L[1]))=f_space_expected_normal_d;
true
gap> Dimension(DefualtFieldOfMatrixGroup(G), Vectorspace(L[2]))=f_space_expected_unitary_d;
true
gap> STOP_TEST("Formspace: Preserved Formspace SU(4,5)");
