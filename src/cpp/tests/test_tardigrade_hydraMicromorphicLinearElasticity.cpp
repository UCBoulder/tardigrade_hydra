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
typedef tardigradeHydra::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

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

BOOST_AUTO_TEST_CASE( testSayHello ){
    /*!
     * Test message printed to stdout in sayHello function
     */

}
