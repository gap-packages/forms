cd ..
COPYFILE_DISABLE=1 tar --exclude=.* -cf forms-`cat ./forms/VERSION`.tar forms/doc forms/lib forms/init.g forms/makedoc.g forms/PackageInfo.g forms/read.g forms/tst forms/VERSION forms/README
zip -rq forms-`cat ./forms/VERSION`-win.zip forms/doc forms/lib forms/init.g forms/makedoc.g forms/PackageInfo.g forms/read.g forms/tst forms/VERSION forms/README
gzip -c9 forms-`cat ./forms/VERSION`.tar > forms-`cat ./forms/VERSION`.tar.gz
bzip2 -9 forms-`cat ./forms/VERSION`.tar
cd forms
