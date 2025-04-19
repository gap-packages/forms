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
#This is the example generating file!


#create workspace with packages
LoadPackage("forms");
SaveWorkspace("forms.ws");
quit;

#restart gap now.

homedir := DirectoryCurrent();

#gapstart := "gap4r13.1"; #might be different on your computer
gapstart := "gap4r14"; #might be different on your computer
gap := Filename(Directory("/usr/local/bin/"),gapstart);
paths := JoinStringsWithSeparator(GAPInfo.RootPaths{[3,4]},";");
pathsstr := Concatenation("\"",paths,"\"");

#note that for the forms package, we use some of the examples as gap source. Therefore we adapted the generate_script:
#when sub = "examples", we point outside the ./tst/gap, instead we point to ./examples/gap
#We make sure that output is always in ./tst/output/sub directory.

generate_script := function(files)
local homedir,scriptfile,sourcedir,str,filename,inputfile,outputfile,outputfilestr,outputdir,cmdlist,cmd;
homedir := DirectoryCurrent();
str := Concatenation("generate_output_forms_examples",".sh");
scriptfile := Filename(homedir,str);
PrintTo(scriptfile,"");
str := "examples/gap/";
sourcedir := DirectoriesPackageLibrary("forms",str)[1];
str := "examples/output/";
outputdir := DirectoriesPackageLibrary("forms",str)[1];
for filename in files do
inputfile := Filename(sourcedir,Concatenation(filename,".g"));
outputfile := Filename(outputdir,Concatenation(filename,".out"));
outputfilestr := Concatenation("\"LogTo(","\\\"",outputfile,"\\\"",");\"");
cmdlist := [gapstart,"-l",pathsstr,"-L","forms.ws","-o","4G","-c",outputfilestr,"<",inputfile,"\n"];
cmd := JoinStringsWithSeparator(cmdlist," ");
AppendTo(scriptfile,cmd);
od;
end;

#initialize filenames. for each filename, there exists a .g file in ./pkg/forms/examples/gap/

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

generate_script(files);

outputdir := DirectoriesPackageLibrary("forms","examples/output")[1];

#Now there is a .sh file ready called generate_output_forms.sh. Make it executable and execute it, then for each filename
#in files an .out file will be generated in ./pkg/forms/examples/gap/output. Each of these .out files contains output almost suitable
#to include in the documentation, or to include in the test directory.

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
