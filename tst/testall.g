LoadPackage("forms");

exclude := [];
if not CompareVersionNumbers( GAPInfo.Version, "4.12" ) then
  Add( exclude, "classic.tst" );
fi;

TestDirectory(DirectoriesPackageLibrary("forms", "tst"),
    rec(
      exitGAP := true,
      exclude := exclude,
      #rewriteToFiles := true,  # enable this line to update tests
    ));
FORCE_QUIT_GAP(1);
