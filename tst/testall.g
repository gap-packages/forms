LoadPackage("forms");
dir := DirectoriesPackageLibrary( "forms", "tst" )[1];

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
"test_recog.tst"
];

for f in testfiles do
	file := Filename(dir,f);
	ReadTest(file);
od;