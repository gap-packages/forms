#############################################################################
##
##  recognition.gd             'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##                                                              Frank Celler
##                      
##  Copyright 2012, Ghent University
##  Copyright 2012, The University of Western Austalia
##  Copyright 1996,  Lehrstuhl D fuer Mathematik, RWTH Aachen, Germany
##
##  This file contains the declarations (functions and operations)
##  for the recognition code.
##
##  *** Bamberg and De Beule are very grateful to Frank Celler for 
##  generously providing the bulk of this code.  ***
##
#############################################################################

#############################################################################
# Functions (not to be used by the user):
#############################################################################

DeclareGlobalFunction( "ClassicalForms_ScalarMultipleFrobenius" );
DeclareGlobalFunction( "ClassicalForms_GeneratorsWithoutScalarsFrobenius" );
DeclareGlobalFunction( "ClassicalForms_ScalarMultipleDual" );
DeclareGlobalFunction( "ClassicalForms_GeneratorsWithoutScalarsDual" );
DeclareGlobalFunction( "ClassicalForms_Signum2" );
DeclareGlobalFunction( "ClassicalForms_Signum" );
DeclareGlobalFunction( "ClassicalForms_QuadraticForm2" );
DeclareGlobalFunction( "ClassicalForms_QuadraticForm" );
DeclareGlobalFunction( "PossibleClassicalForms" );

#############################################################################
# Methods (to be used by the user):
#############################################################################

DeclareOperation( "PreservedSesquilinearForms", [ IsMatrixGroup ] );
DeclareOperation( "ScalarOfSimilarity", [ IsMatrix, IsSesquilinearForm ]);

InfoForms := NewInfoClass("InfoForms");;

