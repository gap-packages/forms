LoadPackage("forms");

exclude := [];
if not IsBound(DescribesInvariantQuadraticForm) then
  # classic.tst should only run in 4.12
  Add( exclude, "adv/classic.tst" );
fi;

if not IsBound(ConformalSymplecticGroup) then
  # conformal.tst cannot run
  Add( exclude, "adv/conformal.tst" );
fi;

if not IsBound(ConformalOrthogonalGroup) then
  # conformal2.tst cannot run
  Add( exclude, "adv/conformal2.tst" );
fi;

TestDirectory(DirectoriesPackageLibrary("forms", "tst"),
    rec(
      exitGAP := true,
      exclude := exclude,
      #rewriteToFile := true,  # enable this line to update tests
    ));
FORCE_QUIT_GAP(1);
