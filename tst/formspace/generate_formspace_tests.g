# you do not need to run this. this creates .tst files that should already exist and be correct!

ReadPackage("forms", "tst/formspace/custom_test_functions.g");
ReadPackage("forms", "tst/interesting_groups.g");

## TODO: Tests with scalars that are not equal to one are desperately missing!!!

WriteTestFilePreservedFormspaceTest := function(path, name, G, n, F,f_space_expected_normal_d, f_space_expected_unitary_d)
    local full_name, start_test, end_test, stream, file_name, dir, full_path;
    full_name := StringFormatted("test_{}.tst", name);
    full_name := ReplacedString(full_name, ",", "_");
    full_path := StringFormatted("{}/{}", path, full_name);

    start_test := StringFormatted("gap> START_TEST(\"Formspace: Preserved Formspace {}\");\n", name);
    end_test := StringFormatted("gap> STOP_TEST(\"Formspace: Preserved Formspace {}\");\n", name);
    dir := DirectoriesPackageLibrary("forms", path);
    # create file with this
    PrintTo(full_path, "");

    file_name := Filename(dir, full_name);

    #write to file with this
    stream := OutputTextFile(file_name, true);
    WriteAll(stream, start_test);
    WriteAll(stream, StringFormatted("gap> G := {};; # Some groups are defined in interesting_groups.g\n", G[1]));
    WriteAll(stream, StringFormatted("gap> R := PseudoRandom(GL({}, {}));;\n", n, F));
    WriteAll(stream, "gap> G := G^R;;\n");
    WriteAll(stream, StringFormatted("gap> f_space_expected_normal_d := {};; # the dimensions of the expected formspaces \n", f_space_expected_normal_d));
    WriteAll(stream, StringFormatted("gap> f_space_expected_unitary_d := {};;\n", f_space_expected_unitary_d));
    WriteAll(stream, "gap> L:=PreservedFormspace(G);;\n");
    WriteAll(stream, "gap> TestMatricesAreForms2(G, L); # found in custom_test_functions.g, Tests whether the matrices returned are forms preserved by G or not\n");
    WriteAll(stream, "[ \"Ok\", \"Ok\" ]\n");
    WriteAll(stream, "gap> Dimension(DefaultFieldOfMatrixGroup(G), Vectorspace(L[1]))=f_space_expected_normal_d;\n");
    WriteAll(stream, "true\n");
    WriteAll(stream, "gap> Dimension(DefaultFieldOfMatrixGroup(G), Vectorspace(L[2]))=f_space_expected_unitary_d;\n");
    WriteAll(stream, "true\n");
    WriteAll(stream, end_test);
    CloseStream(stream);
end;

GenerateTestsForPreservedFormspace := function()
    local Groups, R, G, GG, n, F, Gens, path, lambdas, i, conjugated, f_space_expected_normal_d, f_space_expected_unitary_d;
    path := "tst/formspace/preserved_formspace";
    Groups := [["GO(5,3)", GO(5, 3)], ["SU(4,5)", SU(4, 5)], ["Sp(4,5)", Sp(4, 5)], ["Gtriv", Gtriv], ["G1", G1], ["G2", G2], ["G3", G3], ["G4", G4], ["G5", G5], ["G6", G6], ["G7", G7], ["G8", G8], ["GP22", GP22], ["Group(SP(4,5).1)", Group(SP(4,5).1)]];
    for GG in Groups do
        G := GG[2];
        Gens := GeneratorsOfGroup(G);
        n := DimensionsMat(Gens[1])[1];
        F := DefaultFieldOfMatrixGroup(G);
        R := PseudoRandom(GL(n, F));
        lambdas := [];
        for i in [1..n] do
            Add(lambdas, One(F));
        od;
        conjugated := GG[2]^R;
        f_space_expected_normal_d := Size(TestComputeFormspaceBruteForce(conjugated, lambdas, false));
        f_space_expected_unitary_d := Size(TestComputeFormspaceBruteForce(conjugated, lambdas, true));
        WriteTestFilePreservedFormspaceTest(path, GG[1], GG, n, F, f_space_expected_normal_d, f_space_expected_unitary_d);
    od;
end;


