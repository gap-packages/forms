# Update the examples in the manual to match the output of the GAP
# version used to process this file
example_tree := ExtractExamples( Directory("./doc/"), "forms.xml",[], 500);;
RunExamples(example_tree, rec(
    changeSources:=true,
    width := 80,
));
