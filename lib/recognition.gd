#############################################################################
##
##  recognition.gd             'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##                                                              Frank Celler
##
##  Copyright 2024, Vrije Universiteit Brussel
##  Copyright 2024, The University of Western Austalia
##  Copyright 2024,  Lehrstuhl D fuer Mathematik, RWTH Aachen, Germany
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
DeclareGlobalFunction( "DualFrobeniusGModuleNew" );

DeclareGlobalFunction( "ClassicalForms_InvariantFormDual" );
DeclareGlobalFunction( "ClassicalForms_InvariantFormFrobenius" );

DeclareGlobalFunction( "ClassicalForms_PossibleScalarsSesquilinear" );
DeclareGlobalFunction( "ClassicalForms_GeneratorsWithBetterScalarsSesquilinear" );
DeclareGlobalFunction( "ClassicalForms_InvariantForms" );

DeclareGlobalFunction( "ClassicalForms_ScalarMultipleFrobenius" );
DeclareGlobalFunction( "ClassicalForms_GeneratorsWithoutScalarsFrobenius" );
DeclareGlobalFunction( "ClassicalForms_ScalarMultipleDual" );
DeclareGlobalFunction( "ClassicalForms_GeneratorsWithoutScalarsDual" );
DeclareGlobalFunction( "ClassicalForms_Signum2" );
DeclareGlobalFunction( "ClassicalForms_Signum" );
DeclareGlobalFunction( "ClassicalForms_QuadraticForm2" );
DeclareGlobalFunction( "ClassicalForms_QuadraticForm" );
DeclareGlobalFunction( "PossibleClassicalForms" );

DeclareGlobalFunction( "ScalarsOfPreservedForm" );

DeclareGlobalFunction( "TestPreservedSesquilinearForms" );

#helper function for classic.gi
DeclareGlobalFunction( "TransposedFrobeniusMat" );


#############################################################################
# Methods (to be used by the user):
#############################################################################

DeclareOperation( "ScalarOfSimilarity", [ IsMatrix, IsSesquilinearForm ]);
DeclareOperation( "PreservedFormsOp", [ IsMatrixGroup ] ); #jdb 19/09/2018: was PreservedForms.

DeclareOperation( "PreservedForms", [ IsMatrixGroup ] );
DeclareOperation( "PreservedQuadraticForms", [ IsMatrixGroup ] );
DeclareOperation( "PreservedSesquilinearForms", [ IsMatrixGroup ] );

InfoForms := NewInfoClass("InfoForms");;

