/**
  * \file test_tardigrade_hydraMicromorphicDruckerPragerPlasticity.cpp
  *
  * Tests for tardigrade_hydraMicromorphicDruckerPragerPlasticity
  */

#include<tardigrade_hydraMicromorphicDruckerPragerPlasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphicDruckerPragerPlasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::variableType variableType; //!< Redefinition of the variable type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::variableVector variableVector; //!< Redefinition of the vector of variable types
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::variableMatrix variableMatrix; //!< Redefinition of the matrix of variable types

typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::constantVector constantVector; //!< Redefinition of the vector of constants
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::constantMatrix constantMatrix; //!< Redefinition of the matrix of constants

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

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

}

BOOST_AUTO_TEST_CASE( testSayHello ){
    /*!
     * Test message printed to stdout in sayHello function
     */

    //Setup redirect variables for stdout
    std::stringbuf buffer;
    cout_redirect rd(&buffer);
    boost::test_tools::output_test_stream result;

    //Initialize test variables
    std::string message;
    std::string answer;
    errorOut error = NULL;

    cout_redirect guard( result.rdbuf() );

    //Check normal operation
    message = "World!";
    answer = "Hello World!\n";
    error = tardigradeHydra::sayHello(message);
    BOOST_CHECK( ! error );
    BOOST_CHECK( result.is_equal( answer ) );

    //Reset error code between tests
    error = NULL;

    //Check for "George" error
    message = "George";
    error = tardigradeHydra::sayHello(message);
    BOOST_CHECK( error );

}
