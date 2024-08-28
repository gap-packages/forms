#create output and include files from .g files in the examples/gap directory.
#Executing this file is NOT necessary for the installation of the package "forms", and the
#use of all documentation.
#
#Performing these steps, files will be written in the doc directory of the
#forms package tree. Under UNIX-like operating systems, you need sufficient
#permissions to do. Executing this is NOT necessary for the installation of the
#package "forms".
#
#Messy things happen when you do it, so don't try this at home kids!

#create workspace with packages
LoadPackage("forms");
SaveWorkspace("forms.ws");
quit;

#restart gap now.

#initialize filenames

files := ["conic", "w53", "preservedform", "bg_th_ex1", "bg_th_ex2","bg_th_ex3","bg_th_ex4",
            "bg_th_ex5", "bg_th_ex6", "bg_th_ex7", "bg_th_ex8", "bg_th_ex9", "bilformbymatrix",
            "quadformbymatrix", "hermitianformbymatrix", "bilformbypoly", "quadformbypoly",
            "hermitianformbypoly", "quadformbybilform", "bilformbyquadform", "assocbilform",
            "evalform", "radicalofform", "polyofform", "discofform", "pres_sesforms1",
            "pres_sesforms2", "pres_quadform", "preservedforms2", "basechangehom", "basechangetocanonical",
            "isometriccanonicalform", "quadformfields", "orthogonaltovector",
            "istotallysingular", "istotallyisotropic", "isisotropicvector",
            "issingularvector", "istotallysingular", "scalarfromsim", "trivialform",
            "trivialform_prop", "wittindex", "typeofform", "orthogonaltovector"];

#files := ["preservedforms2"];

homedir := DirectoryCurrent();
scriptfile := Filename(homedir,"generate_output_forms.sh");
PrintTo(scriptfile,"");

gapstart := "gap4r13.1"; #might be different on your computer
gap := Filename(Directory("/usr/local/bin/"),gapstart);
paths := JoinStringsWithSeparator(GAPInfo.RootPaths{[3,4]},";");
pathsstr := Concatenation("\"",paths,"\"");
exampledir := DirectoriesPackageLibrary("forms","examples/gap")[1];
outputdir := DirectoriesPackageLibrary("forms","examples/output")[1];

filename := "conic";
#cmd := ["gap4r13.1 -l "./;/opt/gap-4.13.1/" -L forms.ws -c "LogTo(\"test.out\");" < conic.g]
#gap4r13.1 -l "./;/opt/gap-4.13.1/" -L forms.ws -c "LogTo(\"test.out\");" < conic.g


for filename in files do
inputfile := Filename(exampledir,Concatenation(filename,".g"));
outputfile := Filename(outputdir,Concatenation(filename,".out"));
outputfilestr := Concatenation("\"LogTo(","\\\"",outputfile,"\\\"",");\"");
cmdlist := [gapstart,"-l",pathsstr,"-L","forms.ws","-o","4G","-c",outputfilestr,"<",inputfile,"\n"];
cmd := JoinStringsWithSeparator(cmdlist," ");
AppendTo(scriptfile,cmd);
od;

filename := "conic";

#create .include files
#for the include files, some characters will be translated to suitable xml
#codes, taking more then one character. Therefore we widen the screen a little bit.
#includir: directory containing the include files: ".../pkg/forms/examples/include"
SizeScreen([85,24]);
includedir := DirectoriesPackageLibrary("forms","examples/include")[1];
for filename in files do
  i := Filename(outputdir,Concatenation(filename,".out"));
  o := Filename(includedir,Concatenation(filename,".include"));
  PrintTo(o,"");
  input_stream := InputTextFile(i);
  #ReadLine(input_stream);
  ReadLine(input_stream);
  line := ReadLine(input_stream);
  while line <> "gap> quit;\n" do
    if line <> "\n" then
      line := ReplacedString(line,"\\\n","\n");
      AppendTo(o,ReplacedString(line,"<","&lt;"));
    fi;
    line := ReadLine(input_stream);
  od;
od;
SizeScreen([80,24]);


SizeScreen([85,24]);
includedir := DirectoriesPackageLibrary("forms","examples/include")[1];
for filename in files do
  i := Filename(outputdir,Concatenation(filename,".out"));
  o := Filename(includedir,Concatenation(filename,".include"));
  PrintTo(o,"");
  input_stream := InputTextFile(i);
  ReadLine(input_stream);
  list := [];
  line := ReadLine(input_stream);
  while line <> "gap> quit;\n" do
    if line <> "\n" then
      line := ReplacedString(line,"\\\n","\n");
      #AppendTo(o,ReplacedString(line,"<","&lt;"));
      Add(list,ReplacedString(line,"<","&lt;"));
    fi;
    line := ReadLine(input_stream);
  od;
  l := Length(list);
  new := ReplacedString(list[l],"\n","");
  list[l] := new;
  for line in list do
    AppendTo(o,line);
  od;
od;
SizeScreen([80,24]);


