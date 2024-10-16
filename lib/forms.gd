#############################################################################
##
##  forms.gd              'Forms' package
##                                                              John Bamberg
##                                                              Jan De Beule
##
##  Copyright 2024, Vrije Universiteit Brussel
##  Copyright 2024, The University of Western Austalia
##
##  Declaration matter for quadratic and sesquilinear forms
##
#############################################################################

#############################################################################
# Categories and representations:
#############################################################################

BindGlobal( "FormFamily", NewFamily( "FormFamily" ) );
DeclareCategory( "IsForm", IsAttributeStoringRep );

# jdb 18/09/18 dealing with the evaluation of subspaces of a projective space under forms,
# it seems logical to add a "vectorspace" field to a Forms object upon creation. This
# will make life much easier when incorporating the necessary checks.
DeclareRepresentation( "IsFormRep", IsForm, [ "matrix", "basefield", "type", "vectorspace" ] );

DeclareCategory( "IsQuadraticForm", IsForm );
DeclareCategoryCollections( "IsQuadraticForm" );

DeclareCategory( "IsSesquilinearForm", IsForm );
DeclareCategoryCollections( "IsSesquilinearForm" );

DeclareCategory( "IsBilinearForm", IsSesquilinearForm );
DeclareCategoryCollections( "IsBilinearForm" );

DeclareCategory( "IsHermitianForm", IsSesquilinearForm );
DeclareCategoryCollections( "IsHermitianForm" );

DeclareCategory( "IsTrivialForm", IsForm );
DeclareCategoryCollections( "IsTrivialForm" );

BindGlobal( "QuadraticFormFamily",
            NewFamily( "QuadraticFormFamily", IsObject, IsQuadraticForm ) );
BindGlobal( "QuadraticFormCollFamily", CollectionsFamily(QuadraticFormFamily) );
BindGlobal( "QuadraticFormType",
            NewType( QuadraticFormFamily, IsQuadraticForm and IsFormRep) );

BindGlobal( "SesquilinearFormFamily",
            NewFamily( "SesquilinearFormFamily", IsObject, IsSesquilinearForm ) );
BindGlobal( "SesquilinearFormCollFamily", CollectionsFamily(SesquilinearFormFamily) );
BindGlobal( "SesquilinearFormType",
            NewType( SesquilinearFormFamily, IsSesquilinearForm and IsFormRep) );

BindGlobal( "BilinearFormFamily",
            NewFamily( "BilinearFormFamily", IsObject, IsBilinearForm ) );
BindGlobal( "BilinearFormCollFamily", CollectionsFamily(BilinearFormFamily) );
BindGlobal( "BilinearFormType",
            NewType( BilinearFormFamily, IsBilinearForm and IsFormRep) );

BindGlobal( "HermitianFormFamily",
            NewFamily( "HermitianFormFamily", IsObject, IsHermitianForm ) );
BindGlobal( "HermitianFormCollFamily", CollectionsFamily(HermitianFormFamily) );
BindGlobal( "HermitianFormType",
            NewType( HermitianFormFamily, IsHermitianForm and IsFormRep) );

BindGlobal( "TrivialFormFamily",
            NewFamily( "TrivialFormFamily", IsObject, IsTrivialForm ) );
BindGlobal( "TrivialFormCollFamily", CollectionsFamily(TrivialFormFamily) );
BindGlobal( "TrivialFormType",
            NewType( TrivialFormFamily, IsTrivialForm and IsFormRep) );

#############################################################################
# Constructor operations:
#############################################################################

## the user probably won't use this one. it will not be documented.
DeclareOperation( "FormByMatrix", [IsMatrix and IsFFECollColl, IsField, IsString] );
DeclareOperation( "BilinearFormByMatrixOp", [IsMatrix and IsFFECollColl, IsField] );
DeclareOperation( "QuadraticFormByMatrixOp", [IsMatrix and IsFFECollColl, IsField]);

## For the users...
DeclareOperation( "BilinearFormByMatrix", [IsMatrix and IsFFECollColl, IsField] );
DeclareOperation( "BilinearFormByMatrix", [IsMatrix and IsFFECollColl] );
DeclareOperation( "HermitianFormByMatrix", [IsMatrix and IsFFECollColl, IsField] );
DeclareOperation( "QuadraticFormByMatrix", [IsMatrix and IsFFECollColl, IsField] );
DeclareOperation( "QuadraticFormByMatrix", [IsMatrix and IsFFECollColl] );

## For the users...
DeclareOperation( "BilinearFormByPolynomial",
                      [ IsPolynomial, IsFiniteFieldPolynomialRing, IsInt ] );
DeclareOperation( "QuadraticFormByPolynomial",
                      [ IsPolynomial, IsFiniteFieldPolynomialRing, IsInt ] );
DeclareOperation( "HermitianFormByPolynomial",
                      [ IsPolynomial, IsFiniteFieldPolynomialRing, IsInt ] );

## if no dimension is specified, then we use a natural default (i.e.,
## the dimension is just the number of indeterminates)
## see forms.gi for more information.
DeclareOperation( "BilinearFormByPolynomial",
                      [ IsPolynomial, IsFiniteFieldPolynomialRing] );
DeclareOperation( "QuadraticFormByPolynomial",
                      [ IsPolynomial, IsFiniteFieldPolynomialRing] );
DeclareOperation( "HermitianFormByPolynomial",
                      [ IsPolynomial, IsFiniteFieldPolynomialRing] );

## exploring the nice connections between bilinear and quadratic forms...
DeclareOperation( "BilinearFormByQuadraticForm",
                      [ IsQuadraticForm ] );

DeclareOperation( "QuadraticFormByBilinearForm",
                      [ IsBilinearForm ] );

## helper operations, also not meant for users -> not documented.
## computing the Gram matrix from a polynomial...
DeclareOperation( "UpperTriangleMatrixByPolynomialForForm",
                      [IsPolynomial, IsField, IsInt, IsList] );
DeclareOperation( "GramMatrixByPolynomialForHermitianForm",
                      [IsPolynomial, IsField, IsInt, IsList] );

## computing base changes. See package documentation for information.
DeclareOperation( "BaseChangeOrthogonalBilinear",
                      [IsMatrix and IsFFECollColl, IsField and IsFinite] );
DeclareOperation( "BaseChangeOrthogonalQuadratic",
                      [IsMatrix and IsFFECollColl, IsField and IsFinite] );
DeclareOperation( "BaseChangeSymplectic",
                      [IsMatrix and IsFFECollColl, IsField and IsFinite] );
DeclareOperation( "BaseChangeHermitian",
                      [IsMatrix and IsFFECollColl, IsField and IsFinite] );

#############################################################################
# Functions to support base change operations. (not to be used by the user):
# 29/8/2011: old names are in the commented out part. We changed this on the occasion of the nearby release of gap4r5
#############################################################################

DeclareGlobalFunction( "Forms_SWR" ); #DeclareGlobalFunction( "SWR" );
DeclareGlobalFunction( "Forms_SUM_OF_SQUARES" ); #DeclareGlobalFunction( "SUM_OF_SQUARES" );
DeclareGlobalFunction( "Forms_REDUCE2" ); #DeclareGlobalFunction( "REDUCE2" );
DeclareGlobalFunction( "Forms_REDUCE4" ); #DeclareGlobalFunction( "REDUCE4" );
DeclareGlobalFunction( "Forms_DIFF_2_S" ); #DeclareGlobalFunction( "DIFF_2_S" );
DeclareGlobalFunction( "Forms_RESET" ); #DeclareGlobalFunction( "RESET" );
DeclareGlobalFunction( "Forms_SQRT2" ); #DeclareGlobalFunction( "SQRT2" );
DeclareGlobalFunction( "Forms_PERM_VAR" ); #DeclareGlobalFunction( "PERM_VAR" );
DeclareGlobalFunction( "Forms_C1" ); #DeclareGlobalFunction( "C1" );
DeclareGlobalFunction( "Forms_QUAD_EQ" ); #DeclareGlobalFunction( "QUAD_EQ");
DeclareGlobalFunction( "Forms_HERM_CONJ" ); #DeclareGlobalFunction( "HERM_CONJ" );
DeclareGlobalFunction( "BaseChangeSymplectic_blockchange" ); #a name clash is inprobable, so we leave these unchanged.
DeclareGlobalFunction( "BaseChangeSymplectic_cleanup" );

#############################################################################
# Operations to check input (not for the user):
#############################################################################

DeclareOperation( "FORMS_IsSymplecticMatrix",[IsFFECollColl, IsField] );
DeclareOperation( "FORMS_IsSymmetricMatrix", [IsFFECollColl] );
DeclareOperation( "FORMS_IsHermitianMatrix", [IsFFECollColl, IsField] );

#############################################################################
# User operations:
#############################################################################

## Operations:
DeclareOperation("RadicalOfFormBaseMat", [IsForm]); #not documented

## Attributes
DeclareAttribute( "GramMatrix", IsForm );
DeclareAttribute( "CompanionAutomorphism", IsSesquilinearForm );
DeclareAttribute( "RadicalOfForm", IsForm );
DeclareAttribute( "BaseChangeToCanonical", IsForm );
DeclareAttribute( "WittIndex", IsForm );
DeclareAttribute( "IsometricCanonicalForm", IsForm );
DeclareAttribute( "DiscriminantOfForm", IsForm );
DeclareAttribute( "PolynomialOfForm", IsForm );
DeclareAttribute( "AssociatedBilinearForm", IsQuadraticForm );
# new in forms 1.2.3
DeclareAttribute( "TypeOfForm", IsBilinearForm);
DeclareAttribute( "TypeOfForm", IsHermitianForm);
DeclareAttribute( "TypeOfForm", IsQuadraticForm);

## Properties
DeclareProperty( "IsReflexiveForm", IsForm );
DeclareProperty( "IsAlternatingForm", IsForm );
DeclareProperty( "IsSymmetricForm", IsForm );

DeclareProperty( "IsDegenerateForm", IsForm );
DeclareProperty( "IsSingularForm", IsQuadraticForm );
DeclareProperty( "IsSingularForm", IsTrivialForm );  #new in 1.2.1

DeclareProperty( "IsOrthogonalForm", IsForm );
DeclareProperty( "IsPseudoForm", IsForm );
DeclareProperty( "IsSymplecticForm", IsForm );

DeclareProperty( "IsEllipticForm", IsForm );
DeclareProperty( "IsParabolicForm", IsForm );
DeclareProperty( "IsHyperbolicForm", IsForm );

## More Operations
DeclareOperation( "BaseChangeHomomorphism", [ IsMatrix and IsFFECollColl, IsField ] );

#changed all declarations below on 23/2/9.
DeclareOperation( "EvaluateForm", [ IsSesquilinearForm,
       IsVector and IsFFECollection, IsVector and IsFFECollection]);
DeclareOperation( "EvaluateForm", [ IsTrivialForm,
       IsVector and IsFFECollection, IsVector and IsFFECollection]);
DeclareOperation( "EvaluateForm", [ IsSesquilinearForm, IsFFECollColl, IsFFECollColl]);
DeclareOperation( "EvaluateForm", [ IsTrivialForm, IsFFECollColl, IsFFECollColl]);
DeclareOperation( "EvaluateForm", [ IsQuadraticForm, IsVector and IsFFECollection]);
DeclareOperation( "EvaluateForm", [ IsQuadraticForm, IsFFECollColl]);
DeclareOperation( "EvaluateForm", [ IsTrivialForm, IsVector and IsFFECollection]);

DeclareOperation("OrthogonalSubspaceMat", [IsForm, IsVector and IsFFECollection]);
DeclareOperation("OrthogonalSubspaceMat", [IsForm, IsMatrix]);

DeclareOperation("OrthogonalSubspace", [IsForm, IsVector and IsFFECollection]);
DeclareOperation("OrthogonalSubspace", [IsForm, IsMatrix]);

DeclareOperation("IsIsotropicVector", [IsForm, IsVector and IsFFECollection]);
DeclareOperation("IsTotallyIsotropicSubspace", [IsForm, IsMatrix]);

DeclareOperation("IsSingularVector", [IsQuadraticForm, IsVector and IsFFECollection]);
DeclareOperation("IsTotallySingularSubspace", [IsQuadraticForm, IsMatrix]);


