#############################################################################
##
##  recognition.gd             'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##                                                              Frank Celler
##
##  Copyright 2017, Vrije Universiteit Brussel
##  Copyright 2017, The University of Western Austalia
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

DeclareGlobalFunction( "DualFrobeniusGModule" );

DeclareGlobalFunction( "ClassicalForms_PossibleScalarsSesquilinear" );
DeclareGlobalFunction( "ClassicalForms_GeneratorsWithBetterScalarsSesquilinear" );
DeclareGlobalFunction( "ClassicalForms_InvariantForms" );

DeclareGlobalFunction( "ClassicalForms_Signum2" );
DeclareGlobalFunction( "ClassicalForms_Signum" );
DeclareGlobalFunction( "ClassicalForms_QuadraticForm2" );
DeclareGlobalFunction( "ClassicalForms_QuadraticForm" );
DeclareGlobalFunction( "PossibleClassicalForms" );

DeclareGlobalFunction( "ScalarsOfPreservedForm" );


#############################################################################
# Methods (to be used by the user):
#############################################################################

DeclareOperation( "ScalarOfSimilarity", [ IsMatrix, IsSesquilinearForm ]);
DeclareOperation( "PreservedFormsOp", [ IsMatrixGroup ] ); #jdb 19/09/2018: was PreservedForms.

DeclareOperation( "PreservedForms", [ IsMatrixGroup ] );
DeclareOperation( "PreservedQuadraticForms", [ IsMatrixGroup ] );
DeclareOperation( "PreservedSesquilinearForms", [ IsMatrixGroup ] );

InfoForms := NewInfoClass("InfoForms");;

