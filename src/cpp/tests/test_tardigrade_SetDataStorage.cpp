/**
  * \file test_tardigrade_SetDataStorage.cpp
  *
  * Tests for tardigrade_SetDataStorage
  */

#include<tardigrade_SetDataStorage.h>

#define BOOST_TEST_MODULE test_tardigrade_SetDataStorage
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node

bool tolerantCheck( const std::vector< double > &v1, const std::vector< double > &v2, double eps = 1e-6, double tol = 1e-9 ){

    if ( v1.size( ) != v2.size( ) ){

        return false;

    }

    BOOST_CHECK( v1.size( ) == v2.size( ) );

    const unsigned int len = v1.size( );

    for ( unsigned int i = 0; i < len; i++ ){

        if ( std::fabs( v1[ i ] ) < tol ){

            if ( std::fabs( v1[ i ] - v2[ i ] ) > eps ){

                return false;

            }

            BOOST_CHECK( std::fabs( v1[ i ] - v2[ i ] ) <= eps );

        }
        else{

            if ( ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) > eps ) ||
                 ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) > eps ) ){

                return false;

            }

            BOOST_TEST( v1[ i ] == v2[ i ] );

        }

    }

    return true;

}

BOOST_AUTO_TEST_CASE( test_SetDataStorageBase, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class ResidualMock{

        public:

            ResidualMock( ){ }

            tardigradeHydra::DataStorage< double > myScalarData;

            tardigradeHydra::DataStorage< std::vector< double > > myVectorData;

    };

    ResidualMock residual;

    BOOST_TEST( !residual.myScalarData.first );

    BOOST_TEST( !residual.myVectorData.first );

    std::vector< double > vectorAnswer = { 123.4, 234.5, 345.6 };

    {

        tardigradeHydra::SetDataStorageBase< double > setFloatType( &residual.myScalarData );

        *setFloatType.value = 123.4;

        tardigradeHydra::SetDataStorageBase< std::vector< double > > setVectorType( &residual.myVectorData );

        *setVectorType.value = { 123.4, 234.5, 345.6 };

    }

    BOOST_TEST( residual.myScalarData.first );

    BOOST_TEST( residual.myScalarData.second == 123.4 );

    BOOST_TEST( residual.myVectorData.first );

    BOOST_TEST( residual.myVectorData.second == vectorAnswer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SetDataStorageBase2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class ResidualMock{

        public:

            ResidualMock( ){ }

            tardigradeHydra::DataStorage< double > _myScalarData;

            tardigradeHydra::DataStorage< std::vector< double > > _myVectorData;

            TARDIGRADE_HYDRA_DECLARE_SETDATASTORAGE_GETTER(myScalarData,tardigradeHydra::SetDataStorageBase,double);

            TARDIGRADE_HYDRA_DECLARE_SETDATASTORAGE_GETTER(myVectorData,tardigradeHydra::SetDataStorageBase,std::vector< double >);

    };

    ResidualMock residual;

    BOOST_TEST( !residual._myScalarData.first );

    BOOST_TEST( !residual._myVectorData.first );

    std::vector< double > vectorAnswer = { 123.4, 234.5, 345.6 };

    {

        tardigradeHydra::SetDataStorageBase< double > setFloatType = residual.get_SetDataStorage_myScalarData( );

        *setFloatType.value = 123.4;

        tardigradeHydra::SetDataStorageBase< std::vector< double > > setVectorType = residual.get_SetDataStorage_myVectorData( );

        *setVectorType.value = { 123.4, 234.5, 345.6 };

    }

    BOOST_TEST( residual._myScalarData.first );

    BOOST_TEST( residual._myScalarData.second == 123.4 );

    BOOST_TEST( residual._myVectorData.first );

    BOOST_TEST( residual._myVectorData.second == vectorAnswer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SetDataStorageIteration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class ResidualMock{

        public:

            ResidualMock( ){ }

            tardigradeHydra::DataStorage< double > myScalarData;

            tardigradeHydra::DataStorage< std::vector< double > > myVectorData;

            std::vector< tardigradeHydra::dataBase* > iteration_data;

            void resetIterationData( ){

                for ( auto v = std::begin( iteration_data ); v != std::end( iteration_data ); ++v ){

                    (*v)->clear( );

                }

                iteration_data.clear( );

            }

            void addIterationData( tardigradeHydra::dataBase* v ){

                iteration_data.push_back( v );

            }

    };

    ResidualMock residual;

    BOOST_TEST( !residual.myScalarData.first );

    BOOST_TEST( !residual.myVectorData.first );

    std::vector< double > vectorAnswer = { 123.4, 234.5, 345.6 };

    {

        tardigradeHydra::SetDataStorageIterationBase< ResidualMock, double > setFloatType( &residual.myScalarData, &residual );

        *setFloatType.value = 123.4;

        tardigradeHydra::SetDataStorageIterationBase< ResidualMock, std::vector< double > > setVectorType( &residual.myVectorData, &residual );

        *setVectorType.value = { 123.4, 234.5, 345.6 };

    }

    BOOST_TEST( residual.myScalarData.first );

    BOOST_TEST( residual.myScalarData.second == 123.4 );

    BOOST_TEST( residual.myVectorData.first );

    BOOST_TEST( residual.myVectorData.second == vectorAnswer, CHECK_PER_ELEMENT );

    BOOST_TEST( residual.iteration_data.size( ) == 2 );

    residual.resetIterationData( );

    BOOST_TEST( residual.iteration_data.size( ) == 0 );

    BOOST_TEST( !residual.myScalarData.first );

    BOOST_TEST( !residual.myVectorData.first );

}

namespace tardigradeHydra{

    namespace test_SetDataStorageIteration2{

        class ResidualMock{

            public:

                ResidualMock( ){ }

                tardigradeHydra::DataStorage< double > _myScalarData;

                tardigradeHydra::DataStorage< std::vector< double > > _myVectorData;

                template<typename T>
                class SetDataStorageIteration : public tardigradeHydra::SetDataStorageIterationBase<ResidualMock,T>{
                    using tardigradeHydra::SetDataStorageIterationBase<ResidualMock,T>::SetDataStorageIterationBase;
                };

                TARDIGRADE_HYDRA_DECLARE_CONTROLLED_SETDATASTORAGE_GETTER(myScalarData,SetDataStorageIteration,double,this);

                TARDIGRADE_HYDRA_DECLARE_CONTROLLED_SETDATASTORAGE_GETTER(myVectorData,SetDataStorageIteration,std::vector<double>,this);

                std::vector< tardigradeHydra::dataBase* > iteration_data = { };

                void resetIterationData( ){

                    for ( auto v = std::begin( iteration_data ); v != std::end( iteration_data ); ++v ){

                        (*v)->clear( );

                    }

                    iteration_data.clear( );

                }

                void addIterationData( tardigradeHydra::dataBase* v ){

                    iteration_data.push_back( v );

                }

        };

    }

}

BOOST_AUTO_TEST_CASE( test_SetDataStorageIteration2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::test_SetDataStorageIteration2::ResidualMock residual;

    BOOST_TEST( !residual._myScalarData.first );

    BOOST_TEST( !residual._myVectorData.first );

    std::vector< double > vectorAnswer = { 123.4, 234.5, 345.6 };

    {

        auto setFloatType = residual.get_SetDataStorage_myScalarData( );

        *setFloatType.value = 123.4;

        auto setVectorType = residual.get_SetDataStorage_myVectorData( );

        *setVectorType.value = { 123.4, 234.5, 345.6 };

    }

    BOOST_TEST( residual._myScalarData.first );

    BOOST_TEST( residual._myScalarData.second == 123.4 );

    BOOST_TEST( residual._myVectorData.first );

    BOOST_TEST( residual._myVectorData.second == vectorAnswer, CHECK_PER_ELEMENT );

    BOOST_TEST( residual.iteration_data.size( ) == 2 );

    residual.resetIterationData( );

    BOOST_TEST( residual.iteration_data.size( ) == 0 );

    BOOST_TEST( !residual._myScalarData.first );

    BOOST_TEST( !residual._myVectorData.first );

}
