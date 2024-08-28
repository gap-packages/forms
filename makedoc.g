##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex, dvips
##  
##  Call this with GAP.
##
## Note that this will write files in the doc directory of the forms package
## tree. Under decent operating systems, you need sufficient permissions to do.
## Executing this is NOT necessary for the installation of the package "forms".
## Run this from the directory of the package using "gap makedoc.g"

if fail = LoadPackage("AutoDoc", ">= 2019.04.10") then
    Error("AutoDoc 2019.04.10 or newer is required");
fi;
AutoDoc(rec(
    scaffold := rec( MainPage := false ),
    gapdoc := rec( main := "forms.xml" ),
));
QUIT;

