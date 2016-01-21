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
            "pres_sesforms2", "basechangehom", "basechangetocanonical",
            "isometriccanonicalform", "quadformfields", "orthogonaltovector",
            "istotallysingular", "istotallyisotropic", "isisotropicvector",
            "issingularvector", "istotallysingular", "scalarfromsim", "trivialform",
            "trivialform_prop", "wittindex", "typeofform"];

#initialize directorynames
#exampledir = dir where .g files are located : ".../pkg/forms/examples/gap"
#preambledir = directory where 'preamble_sws.g is found' :  ".../pkg/forms/examples"
#outputdir = directory to write '.out' files: ".../pkg/forms/examples/output"
#name of script to start gap version. The user has to fill this in!

homedir := DirectoryCurrent();
exampledir := DirectoriesPackageLibrary("forms","examples/gap")[1];
preambledir := DirectoriesPackageLibrary("forms","examples/")[1];
outputdir := DirectoriesPackageLibrary("forms","examples/output")[1];
gap := Filename(Directory("/usr/bin/"),"gap4r7");
#paths := JoinStringsWithSeparator(GAPInfo.RootPaths{[2,3]},";");
paths := JoinStringsWithSeparator(GAPInfo.RootPaths{[3,4]},";");
args := JoinStringsWithSeparator(["-l",paths," -L forms.ws"," -o 4G"]," ");
args := ["-l",paths,"-L","forms.ws","-o","4G"];
extension := ".out\";";
cmddir := "dir \:\= DirectoriesPackageLibrary\(\"forms\"\,\"examples\/output\"\)\[1\]\;";

#create .out files using the saved workspace
#IMPORTANT: here we suppose that the script to start up our favorite version of
#GAP is called 'gap4r4', and is located in '/usr/bin'. Change the code if this is not true!
#you certainly now the name of the script, since you started gap. To find the
#dir, just issue in the gap session that is running:

#Exec("which gap4r4"); #for UNIX only

gapstart := "gap4r8"; #might be different on your computer
gap := Filename(Directory("/usr/bin/"),gapstart);

for filename in files do
  Print("Now converting file: ", filename, "\n");
  stream := InputOutputLocalProcess( homedir, gap, args);
  #cmd := Concatenation("file := \"",filename,".out\";");
  cmd := Concatenation("file := \"",filename,extension);
  WriteLine(stream,cmd);
  #cmd := "dir \:\= DirectoriesPackageLibrary\(\"forms\"\,\"examples\/output\"\)\[1\]\;";
  WriteLine(stream,cmddir);
  preamble := Filename(preambledir,"preamble_sws.g");
  preamble_stream := InputTextFile(preamble);
  cmds := ReadAll(preamble_stream);
  WriteLine(stream,cmds);
  repeat
    str := ReadLine(stream);
  until str = "true\n";
  inputfile := Filename(exampledir,Concatenation(filename,".g"));
  input_stream := InputTextFile(inputfile);
  cmd := ReadLine(input_stream);
  while cmd <> fail do
    WriteAll(stream,cmd);
    cmd := ReadLine(input_stream);
    ReadAll(stream);
  od;
  repeat until ReadAll(stream)=fail; #new since oct 2015.
od;

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


