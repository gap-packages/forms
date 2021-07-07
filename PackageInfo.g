#############################################################################
##  
##  PackageInfo.g for the package `Forms'                 
##                                                              John Bamberg
##                                                              Jan De Beule
##
##  (created from Frank Luebeck's PackageInfo.g template file)
##  

SetPackageInfo( rec( 
  PackageName := "Forms", 
  Subtitle := "Sesquilinear and Quadratic",
  Version := "1.2.6",
  Date := "07/07/2021",
License := "GPL-2.0-or-later",

##  URL of the archive(s) of the current package release, but *without*
##  the format extension(s), like '.zoo', which are given next.
##  The archive file name *must be changed* with each version of the archive
##  (and probably somehow contain the package name and version).
##  The paths of the files in the archive must begin with the name of the
##  directory containing the package (in our "example" probably:
##  example/init.g, ...    or  example-1.3/init.g, ...  )
# 
ArchiveURL := "http://cage.ugent.be/geometry/software/forms/forms-1.2.6",
ArchiveFormats := ".tar.gz -win.zip .tar.bz2",
Persons := [
  rec( 
    LastName      := "Bamberg",
    FirstNames    := "John",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "bamberg@maths.uwa.edu.au",
    WWWHome       := "http://school.maths.uwa.edu.au/~bamberg/",
    PostalAddress := Concatenation( [
                       "John Bamberg\n",
                       "School of Mathematics and Statistics\n",
                       "The University of Western Australia\n",
                       "35 Stirling Highway\n",
                       "CrawleyY WA 6009, Perth\n",
                       "Australia" ] ),
    Place         := "Perth",
    Institution   := "The University of Western Australia",
  ),
  rec( 
    LastName      := "De Beule",
    FirstNames    := "Jan",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "jan@debeule.eu",
    WWWHome       := "http://www.debeule.eu",
    PostalAddress := Concatenation( [
                       "Jan De Beule\n",
                       "Department of Mathematics\n",
                       "Vrije Universiteit Brussel\n",
                       "Pleinlaan 2\n",
                       "B-1050 Brussel\n",
                       "Belgium" ] ),
    Place         := "Brussels",
    Institution   := "Vrije Universiteit Brussel",
  ),
],

Status := "accepted",
README_URL := "http://cage.ugent.be/geometry/software/forms/README",
PackageInfoURL := "http://cage.ugent.be/geometry/software/forms/PackageInfo.g",
AbstractHTML := "This package can be used for work with sesquilinear and quadratic forms on finite vector spaces; objects that are used to describe polar spaces and classical groups.",

PackageWWWHome := "http://cage.ugent.be/geometry/forms.php",
            
PackageDoc := rec(
  # use same as in GAP            
  BookName  := "Forms",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Forms - Sesquilinear and Quadratic",
),

Dependencies := rec(
  GAP := ">=4.8",
  NeededOtherPackages := [],
  SuggestedOtherPackages := [],
  ExternalConditions := []),

AvailabilityTest := ReturnTrue,

BannerString := Concatenation( 
  "---------------------------------------------------------------------\n",
  "Loading 'Forms' ", ~.Version," (", ~.Date,")", "\n",
  "by ", ~.Persons[1].FirstNames, " ", ~.Persons[1].LastName,
        " (", ~.Persons[1].WWWHome, ")\n",
  "   ", ~.Persons[2].FirstNames, " ", ~.Persons[2].LastName,
        " (", ~.Persons[2].WWWHome, ")\n",
   "For help, type: ?Forms \n",
  "---------------------------------------------------------------------\n" ),

TestFile := "tst/testall.g",

IssueTrackerURL := "https://github.com/gap-packages/forms/issues",
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/gap-packages/forms",
),

Keywords := ["Forms", "Sesquilinear", "Quadratic"],

CommunicatedBy := "Leonard Soicher (London)",
AcceptDate := "03/2009"

));


