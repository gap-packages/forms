TestMyPackage := function( pkgname )
local pkgdir, testfiles, testresult, ff, fn;
LoadPackage( pkgname );
pkgdir := DirectoriesPackageLibrary( pkgname, "tst" );

# Arrange chapters as required
testfiles := [
"test_forms1.tst",
"test_forms2.tst",
"test_forms3.tst",
"test_forms4.tst",
"test_forms5.tst",
"test_forms6.tst",
"test_forms7.tst",
"test_forms8.tst",
"test_forms9.tst",
"test_forms10.tst",
"test_forms11.tst",
"test_recog.tst",
"test_forms12.tst",
"test_forms13.tst",
"test_forms14.tst",
"test_forms15.tst",
"test_forms16.tst"
];

testresult:=true;
for ff in testfiles do
  fn := Filename( pkgdir, ff );
  Print("#I  Testing ", fn, "\n");
  if not Test( fn, rec(compareFunction := "uptowhitespace") ) then
    testresult:=false;
  fi;
od;
if testresult then
  Print("#I  No errors detected while testing package ", pkgname, "\n");
else
  Print("#I  Errors detected while testing package ", pkgname, "\n");
fi;
end;

# Set the name of the package here
Print("EXecuting this file\n");

TestMyPackage( "forms" );

QUIT_GAP(0);

