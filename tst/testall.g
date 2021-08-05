LoadPackage("forms");

exclude := [];
if not IsBound(DescribesInvariantQuadraticForm) then
  # classic.tst should only run in 4.12
  Add( exclude, "classic.tst" );
fi;

TestDirectory(DirectoriesPackageLibrary("forms", "tst"),
    rec(
      exitGAP := true,
      exclude := exclude,
      #rewriteToFiles := true,  # enable this line to update tests
    ));
FORCE_QUIT_GAP(1);
