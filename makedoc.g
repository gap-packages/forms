##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex, dvips
##  
##  Call this with GAP.
##
## Note that this will write files in the doc directory of the forms package
## tree. Under decent operating systems, you need sufficient permissions to do.
## Executing this is NOT necessary for the installation of the package "forms".

if fail = LoadPackage("AutoDoc", ">= 2016.01.21") then
    Error("AutoDoc 2016.01.21 or newer is required");
fi;
AutoDoc("forms");
QUIT;

