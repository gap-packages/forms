#We need several files to reach a good code coverage. A file contains only gap commands and a quit;
#It seems not advisable to test Error messages as well.
#The objective of this file is to generate all necessary tst files from files containing gap commands.
#I experimented a lot with starting gap slaves and talking to them using streams, this is quite complicated.
#The new appraoch (since August 2024) is to generate a .sh script, starting one gap for each .g file, and
#creating .out files that can be processed into .tst files.
#Executing this file is NOT necessary for the installation of the package "forms" (including the documentation)

#to make sure tjere is no confusing with linebreaks, it seems better to have all input (in the .g files)
#in one (long) line.

#create workspace with packages
LoadPackage("forms");
SaveWorkspace("forms.ws");
quit;

#restart gap now.

homedir := DirectoryCurrent();
scriptfile := Filename(homedir,"generate_output_forms_testfiles.sh");
PrintTo(scriptfile,"");

#gapstart := "gap4r13.1"; #might be different on your computer
gapstart := "gap4r14"; #might be different on your computer
gap := Filename(Directory("/usr/local/bin/"),gapstart);
paths := JoinStringsWithSeparator(GAPInfo.RootPaths{[3,4]},";");
pathsstr := Concatenation("\"",paths,"\"");

#note that for the forms package, we use some of the examples as gap source. Therefore we adapted the generate_script:
#when sub = "examples", we point outside the ./tst/gap, instead we point to ./examples/gap
#We make sure that output is always in ./tst/output/sub directory.

generate_script := function(sub,files)
local homedir,scriptfile,sourcedir,str,filename,inputfile,outputfile,outputfilestr,outputdir,cmdlist,cmd;
homedir := DirectoryCurrent();
str := Concatenation("generate_output_forms_testfiles_",sub,".sh");
scriptfile := Filename(homedir,str);
PrintTo(scriptfile,"");
if sub = "examples" then
str := "examples/gap";
else
str := Concatenation("tst/gap/",sub);
fi;
sourcedir := DirectoriesPackageLibrary("forms",str)[1];
str := Concatenation("tst/output/",sub);
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

#We have several subdirectories with .g files to convert.

subs := ["examples", "easy", "adv" ];

#examples files
filesexamples := ["w53", "bg_th_ex3", "bg_th_ex4",
            "bg_th_ex5", "bg_th_ex6", "bg_th_ex8", "bg_th_ex9", "bilformbymatrix",
            "quadformbymatrix", "hermitianformbymatrix", "quadformbypoly",
            "hermitianformbypoly", "quadformbybilform", "bilformbyquadform", "assocbilform",
            "evalform", "polyofform", "discofform", "pres_quadform", "preservedforms2", "basechangehom", "basechangetocanonical",
            "isometriccanonicalform", "quadformfields", "orthogonaltovector",
            "istotallysingular", "istotallyisotropic", "isisotropicvector",
            "issingularvector", "istotallysingular", "scalarfromsim", "trivialform",
            "trivialform_prop", "wittindex", "orthogonaltovector"];

#easy test files
fileseasy := ["test_forms1", "test_forms2", "test_forms3", "test_forms4", "test_forms5",
            "test_forms6","test_forms7", "test_forms8", "test_forms9", "test_forms10",
            "test_forms11", "test_recog", "test_forms12", "test_forms13", "test_forms14",
            "test_forms15", "test_forms16", "test_tech1", "test_tech2", "test_tech3", "test_radicalofform",
            "test_bg_th_ex7", "test_bg_th_ex2", "test_conic" ];
           
#advanced tests   #"basechange", "classic",
filesadv := ["test_recog", "test_preservedform", "test_pres_sesforms1", "test_pres_sesforms2"];

generate_script("examples",filesexamples);
generate_script("easy",fileseasy);
generate_script("adv",filesadv);

#make your choice: for "easy" and "advanced":
sourcedir := DirectoriesPackageLibrary("forms","tst/gap")[1];

#for examples:
sourcedir := DirectoriesPackageLibrary("forms","examples/gap")[1];

#continue now
outputdir := DirectoriesPackageLibrary("forms","tst/output")[1];

for filename in files do
inputfile := Filename(sourcedir,Concatenation(filename,".g"));
outputfile := Filename(outputdir,Concatenation(filename,".out"));
outputfilestr := Concatenation("\"LogTo(","\\\"",outputfile,"\\\"",");\"");
cmdlist := [gapstart,"-l",pathsstr,"-L","forms.ws","-o","4G","-c",outputfilestr,"<",inputfile,"\n"];
cmd := JoinStringsWithSeparator(cmdlist," ");
AppendTo(scriptfile,cmd);
od;

preambledir := DirectoriesPackageLibrary("forms","examples/")[1];
outputdir := DirectoriesPackageLibrary("forms","tst/output")[1];
cmddir := "dir \:\= DirectoriesPackageLibrary\(\"forms\"\,\"tst\/output\"\)\[1\]\;";

create_tst_files := function(sub,files)
local includedir,str,filename,i,o,input_stream,line,outputdir;
str := Concatenation("tst/",sub);
includedir := DirectoriesPackageLibrary("forms",str)[1];
str := Concatenation("tst/output/",sub);
outputdir := DirectoriesPackageLibrary("forms",str)[1];
for filename in files do
  i := Filename(outputdir,Concatenation(filename,".out"));
  o := Filename(includedir,Concatenation(filename,".tst"));
  PrintTo(o,"");
  input_stream := InputTextFile(i);
  ReadLine(input_stream); #reads first line which is a line with a #
  #ReadLine(input_stream);
  AppendTo(o,Concatenation("gap> START_TEST(\"Forms: ",filename,".tst\");\n"));
  line := ReadLine(input_stream);
  while line <> "gap> quit;\n" do
    if Length(line) > 3 then
        if line{[1..4]} = "gap>" then
            RemoveCharacters(line,"\n");
            line := Concatenation(line,"\n");
        fi;
    fi;
    SizeScreen([500,24]);
    AppendTo(o,line);
    SizeScreen([80,24]);
    line := ReadLine(input_stream);
  od;
  AppendTo(o,Concatenation("gap> STOP_TEST(\"",filename,".tst\", 10000 );\n"));
od;
end;

subs := ["examples", "easy", "adv" ];

create_tst_files("examples",filesexamples);
create_tst_files("easy",fileseasy);
create_tst_files("adv",filesadv);
