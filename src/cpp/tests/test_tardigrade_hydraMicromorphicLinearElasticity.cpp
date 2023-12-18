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

typedef tardigradeHydra::micromorphicLinearElasticity::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicLinearElasticity::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicLinearElasticity::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicLinearElasticity::constantVector constantVector; //!< Redefinition of the vector of constants

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
