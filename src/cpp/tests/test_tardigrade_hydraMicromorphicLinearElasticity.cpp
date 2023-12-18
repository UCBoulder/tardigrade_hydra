/**
  * \file test_tardigrade_hydraMicromorphicLinearElasticity.cpp
  *
  * Tests for tardigrade_hydraMicromorphicLinearElasticity
  */

#include<tardigrade_hydraMicromorphicLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphicLinearElasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::micromorphicLinearElasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::micromorphicLinearElasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::micromorphicLinearElasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

typedef tardigradeHydra::micromorphicLinearElasticity::variableType variableType; //!< Redefinition of the variable type
typedef tardigradeHydra::micromorphicLinearElasticity::variableVector variableVector; //!< Redefinition of the vector of variable types
typedef tardigradeHydra::micromorphicLinearElasticity::variableMatrix variableMatrix; //!< Redefinition of the matrix of variable types

typedef tardigradeHydra::micromorphicLinearElasticity::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicLinearElasticity::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicLinearElasticity::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicLinearElasticity::constantVector constantVector; //!< Redefinition of the vector of constants
typedef tardigradeHydra::micromorphicLinearElasticity::constantMatrix constantMatrix; //!< Redefinition of the matrix of constants

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( testFormIsotropicA ){
    /*!
     * Test the formation of the isotropic A stiffness tensor.
     *
     */

    parameterType lambda = 4;
    parameterType mu = 7;

    parameterVector answer = { 18.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,  7.,  0.,  7.,
                                0.,  0.,  0.,  0.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  7.,  0.,
                                0.,  0.,  7.,  0.,  7.,  0.,  0.,  0.,  0.,  0.,  4.,  0.,  0.,
                                0., 18.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  0.,  0.,  7.,  0.,
                                7.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.,
                                0.,  0.,  0.,  7.,  0.,  7.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,
                                0.,  0., 18. };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicA( lambda, mu, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testFormIsotropicB ){
    /*!
     * Test the formation of the isotropic B stiffness tensor.
     *
     */

    parameterType eta = 3;
    parameterType tau = 5;
    parameterType kappa = 6;
    parameterType nu = 8;
    parameterType sigma = 4;

    parameterVector answer = {  4.,  0.,  0.,  0., -2.,  0.,  0.,  0., -2.,  0.,  2.,  0.,  4.,
                                0.,  0.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  4.,  0.,
                                0.,  0.,  4.,  0.,  2.,  0.,  0.,  0.,  0.,  0., -2.,  0.,  0.,
                                0.,  4.,  0.,  0.,  0., -2.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,
                                4.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  0.,
                                0.,  0.,  0.,  4.,  0.,  2.,  0., -2.,  0.,  0.,  0., -2.,  0.,
                                0.,  0.,  4. };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicB( eta, tau, kappa, nu, sigma, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testFormIsotropicC ){
    /*!
     * Test the formation of the isotropic C stiffness tensor.
     *
     */

    parameterVector taus = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    parameterVector answer = {
        97.,  0.,  0.,  0., 13.,  0.,  0.,  0., 13.,  0., 16.,  0.,  9.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  9.,  0.,
        0.,  0., 23.,  0., 22.,  0.,  0.,  0.,  0.,  0., 23.,  0.,  0.,
        0.,  9.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,
        4.,  0.,  0.,  0., 23.,  0.,  0.,  0., 22.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  4.,  0.,  3.,  0., 23.,  0.,  0.,  0.,  2.,  0.,
        0.,  0.,  9.,  0., 22.,  0., 27.,  0.,  0.,  0.,  0.,  0., 26.,
        0.,  0.,  0., 16.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,
        7.,  0.,  3.,  0., 13.,  0.,  0.,  0., 23.,  0.,  0.,  0.,  5.,
        0., 26.,  0., 23.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,
        0.,  0.,  2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  8.,  0., 10.,
        0.,  0.,  0., 11.,  0.,  0.,  0.,  9.,  0.,  0.,  0.,  9.,  0.,
       12.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 22.,  0.,  0.,  0., 27.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,  7.,  0., 26.,  0.,
        0.,  0.,  6.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  0.,  0., 10.,
        0.,  8.,  0.,  0.,  0.,  9.,  0.,  0.,  0., 12.,  0.,  0.,  0.,
       11.,  0.,  9.,  0.,  0.,  0.,  0.,  0., 13.,  0.,  0.,  0.,  5.,
        0.,  0.,  0., 23.,  0.,  6.,  0.,  2.,  0.,  0.,  0.,  0.,  0.,
        0.,  0., 26.,  0.,  0.,  0., 23.,  0.,  0.,  0., 23.,  0., 26.,
        0.,  0.,  0.,  0.,  0., 23.,  0.,  0.,  0., 13.,  0.,  0.,  0.,
        5.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,  2.,  0., 16.,  0.,  0.,
        0., 26.,  0.,  0.,  0.,  6.,  0., 27.,  0., 22.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,
        0.,  0.,  0., 11.,  0.,  9.,  0.,  0.,  0.,  8.,  0.,  0.,  0.,
       10.,  0.,  0.,  0., 12.,  0.,  9.,  0.,  0.,  0.,  0.,  0.,  9.,
        0.,  0.,  0., 23.,  0.,  0.,  0.,  2.,  0., 22.,  0., 23.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  4.,  0.,  0.,
        0.,  9.,  0., 16.,  0.,  0.,  0.,  0.,  0., 13.,  0.,  0.,  0.,
       97.,  0.,  0.,  0., 13.,  0.,  0.,  0.,  0.,  0., 16.,  0.,  9.,
        0.,  0.,  0.,  4.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,
        0.,  0., 23.,  0., 22.,  0.,  2.,  0.,  0.,  0., 23.,  0.,  0.,
        0.,  9.,  0.,  0.,  0.,  0.,  0.,  9.,  0., 12.,  0.,  0.,  0.,
       10.,  0.,  0.,  0.,  8.,  0.,  0.,  0.,  9.,  0., 11.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,
        0.,  0.,  0.,  0., 22.,  0., 27.,  0.,  6.,  0.,  0.,  0., 26.,
        0.,  0.,  0., 16.,  0.,  2.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,
        5.,  0.,  0.,  0., 13.,  0.,  0.,  0., 23.,  0.,  0.,  0.,  0.,
        0., 26.,  0., 23.,  0.,  0.,  0., 23.,  0.,  0.,  0., 26.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,  6.,  0., 23.,  0.,  0.,
        0.,  5.,  0.,  0.,  0., 13.,  0.,  0.,  0.,  0.,  0.,  9.,  0.,
       11.,  0.,  0.,  0., 12.,  0.,  0.,  0.,  9.,  0.,  0.,  0.,  8.,
        0., 10.,  0.,  0.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  6.,  0.,
        0.,  0., 26.,  0.,  7.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  0.,
        0., 27.,  0.,  0.,  0., 22.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
       12.,  0.,  9.,  0.,  0.,  0.,  9.,  0.,  0.,  0., 11.,  0.,  0.,
        0., 10.,  0.,  8.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,
        0.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 23.,  0., 26.,
        0.,  5.,  0.,  0.,  0., 23.,  0.,  0.,  0., 13.,  0.,  3.,  0.,
        7.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,  0.,  0., 16.,  0.,  0.,
        0., 26.,  0.,  0.,  0.,  0.,  0., 27.,  0., 22.,  0.,  9.,  0.,
        0.,  0.,  2.,  0.,  0.,  0., 23.,  0.,  3.,  0.,  4.,  0.,  0.,
        0.,  0.,  0.,  0.,  0., 22.,  0.,  0.,  0., 23.,  0.,  0.,  0.,
        4.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  9.,
        0.,  0.,  0., 23.,  0.,  0.,  0.,  0.,  0., 22.,  0., 23.,  0.,
        0.,  0.,  9.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  9.,  0., 16.,  0., 13.,  0.,  0.,  0., 13.,  0.,  0.,  0.,
       97.
    };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicC( taus, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testFormIsotropicD ){
    /*!
     * Test the formation of the isotropic D stiffness tensor.
     *
     */

    parameterType tau = 5;
    parameterType sigma = 4;

    parameterVector answer = { 13.,  0.,  0.,  0.,  5.,  0.,  0.,  0.,  5.,  0.,  4.,  0.,  4.,
                                0.,  0.,  0.,  0.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,
                                0.,  0.,  4.,  0.,  4.,  0.,  0.,  0.,  0.,  0.,  5.,  0.,  0.,
                                0., 13.,  0.,  0.,  0.,  5.,  0.,  0.,  0.,  0.,  0.,  4.,  0.,
                                4.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  0.,
                                0.,  0.,  0.,  4.,  0.,  4.,  0.,  5.,  0.,  0.,  0.,  5.,  0.,
                                0.,  0., 13. };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicD( tau, sigma, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testAssembleFundamentalDeformationMeasures ){
    /*!
     * Assemble the fundamental deformation measures from the degrees of freedom.
     *
     */

    double grad_u[ 3 ][ 3 ] = { { 1, 2, 3 },
                                { 4, 5, 6 },
                                { 7, 8, 9 } };

    double phi[ 9 ] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    double grad_phi[ 9 ][ 3 ] = { {  1,  2,  3 },
                                  {  4,  5,  6 },
                                  {  7,  8,  9 },
                                  { 10, 11, 12 },
                                  { 13, 14, 15 },
                                  { 16, 17, 18 },
                                  { 19, 20, 21 },
                                  { 22, 23, 24 },
                                  { 25, 26, 27 } };

    variableVector answerDeformationGradient = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerMicroDeformation = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerGradientMicroDeformation = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                     10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                     19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector resultF, resultChi, resultGradChi;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                                                            resultF, resultChi, resultGradChi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultF, answerDeformationGradient ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultChi, answerMicroDeformation ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultGradChi, answerGradientMicroDeformation ) );

    //Test the Jacobians
    variableVector resultFJ, resultChiJ, resultGradChiJ;
    variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

    error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                                                   resultFJ, resultChiJ, resultGradChiJ,
                                                                                                   dFdGradU, dChidPhi, dGradChidGradPhi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultFJ, answerDeformationGradient ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultChiJ, answerMicroDeformation ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultGradChiJ, answerGradientMicroDeformation ) );

    //Test the jacobians w.r.t. the gradient of the displacement
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta( 3, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_u[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 3 ][ 3 ] =
        {
                { grad_u[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] + delta[ 2 ][ 2 ] }
        };

        double negative_perturb[ 3 ][ 3 ] =
        {
                { grad_u[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] - delta[ 2 ][ 2 ] }
        };

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( positive_perturb, phi, grad_phi,
                                                                                                       FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( negative_perturb, phi, grad_phi,
                                                                                                       FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dFdGradU[ j ][ i ] ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 9; i++ ){
        constantVector delta( 9, 0 );

        delta[ i ] = eps * fabs( phi[ i ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ] = { phi[ 0 ] + delta[ 0 ], phi[ 1 ] + delta[ 1 ], phi[ 2 ] + delta[ 2 ],
                                         phi[ 3 ] + delta[ 3 ], phi[ 4 ] + delta[ 4 ], phi[ 5 ] + delta[ 5 ],
                                         phi[ 6 ] + delta[ 6 ], phi[ 7 ] + delta[ 7 ], phi[ 8 ] + delta[ 8 ] };

        double negative_perturb[ 9 ] = { phi[ 0 ] - delta[ 0 ], phi[ 1 ] - delta[ 1 ], phi[ 2 ] - delta[ 2 ],
                                         phi[ 3 ] - delta[ 3 ], phi[ 4 ] - delta[ 4 ], phi[ 5 ] - delta[ 5 ],
                                         phi[ 6 ] - delta[ 6 ], phi[ 7 ] - delta[ 7 ], phi[ 8 ] - delta[ 8 ] };

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, positive_perturb, grad_phi,
                                                                                                       FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, negative_perturb, grad_phi,
                                                                                                       FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dChidPhi[ j ][ i ] ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_u[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ][ 3 ] =
        {
                { grad_phi[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] + delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] + delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] + delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] + delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] + delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] + delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] + delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] + delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] + delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] + delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] + delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] + delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] + delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] + delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] + delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] + delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] + delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] + delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] + delta[ 8 ][ 2 ] }
        };

        double negative_perturb[ 9 ][ 3 ] =
        {
                { grad_phi[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] - delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] - delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] - delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] - delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] - delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] - delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] - delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] - delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] - delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] - delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] - delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] - delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] - delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] - delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] - delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] - delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] - delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] - delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] - delta[ 8 ][ 2 ] }
        };

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, positive_perturb,
                                                                                      FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, negative_perturb,
                                                                                                       FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dGradChidGradPhi[ j ][ i ] ) );
        }
    }
}
