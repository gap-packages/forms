##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex, dvips
##  
##  Call this with GAP.
##
## Note that this will write files in the doc directory of the forms package
## tree. Under decent operating systems, you need sufficient permissions to do.
## Executing this is NOT necessary for the installation of the package "forms".

# probably not necessary since GAPDoc is autoloading:
LoadPackage("GAPDoc");

#initialize the directory /doc/ directory in the package tree
docdir := DirectoriesPackageLibrary("forms","doc")[1];

MakeGAPDocDoc(docdir, "forms", [], "Forms", "MathJax");

GAPDocManualLab("forms");

quit;

