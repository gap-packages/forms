#!/usr/bin/env bash

gap4r9 -l "./;/opt/gap-4.9.3/" -a 500M -m 500M -q -A --cover ./pkg/forms/tmp.json ~/pkg/forms/tst/testall.g

gap4r9 -l "./;/opt/gap-4.9.3/" -a 500M -m 500M <<GAPInput
if LoadPackage("profiling") <> true then
    Print("ERROR: could not load profiling package");
    FORCE_QUIT_GAP(1);
fi;
x := ReadLineByLineProfile("./pkg/forms/tmp.json");;
OutputAnnotatedCodeCoverageFiles(x,"./pkg/forms","./pkg/forms/tmp"); 
QUIT_GAP(0);
GAPInput

