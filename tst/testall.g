LoadPackage("forms");
TestDirectory(DirectoriesPackageLibrary("forms", "tst"),
    rec(
      exitGAP := true,
      #rewriteToFiles := true,  # enable this line to update tests
    ));
FORCE_QUIT_GAP(1);
