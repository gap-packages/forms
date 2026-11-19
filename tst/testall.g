LoadPackage("forms");

tstdir := DirectoriesPackageLibrary("forms", "tst")[1];
# all .tst files below tst/<sub>, as paths relative to tst/ (AI code)
TstFilesBelow := function(sub)
  local res, path, f, full;
  res := [];
  path := Filename(tstdir, sub);
  for f in DirectoryContents(path) do
    if f in [".", ".."] then continue; fi;
    full := Concatenation(sub, "/", f);
    if IsDirectoryPath(Filename(tstdir, full)) then
      Append(res, TstFilesBelow(full));
    elif EndsWith(f, ".tst") then
      Add(res, full);
    fi;
  od;
  return res;
end;

exclude := [];
if not IsBound(DescribesInvariantQuadraticForm) then
  # classic.tst should only run in 4.12
  Add( exclude, "adv/classic.tst" );
fi;

ReadPackage("forms", "tst/interesting_groups.g");
ReadPackage("forms", "tst/formspace/custom_test_functions.g");
if not IsBound(ConformalSymplecticGroup) then
  # conformal.tst cannot run
  Add( exclude, "adv/conformal.tst" );
fi;

# exclude all matrixobj tests for older gap versions. This is needed because the way that matrix objs are printed was changed.
if not CompareVersionNumbers(GAPInfo.Version, "4.17dev") then
  Append( exclude, TstFilesBelow("matobj") );
fi;

# Currently this test is horribly slow and can propabliy only really be run once it works well to check G1 = G2 for two groups made up of matrix objects.
Add(exclude, "matobj/duplicates/adv/classic.tst");

TestDirectory(DirectoriesPackageLibrary("forms", "tst"),
    rec(
      exitGAP := true,
      exclude := exclude,
      #rewriteToFile := true,  # enable this line to update tests
    ));
FORCE_QUIT_GAP(1);
