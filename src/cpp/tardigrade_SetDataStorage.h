/**
  ******************************************************************************
  * \file tardigrade_SetDataStorage.h
  ******************************************************************************
  * A C++ library for constructing caching data storage objects
  ******************************************************************************
  */

#ifndef TARDIGRADE_SETDATASTORAGE
#define TARDIGRADE_SETDATASTORAGE

#include"tardigrade_error_tools.h"
#include"tardigrade_CoreDefinitions.h"
#include"Eigen/Dense"

/*!
 * \brief Declares a named setDataStorage getter function
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param classtype: The type of the class to be defined
 * \param vartype: The ctype of the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_SETDATASTORAGE_GETTER(getname,varname,classtype,vartype,...) \
    classtype< vartype > getname( ){                                                          \
        /*!                                                                                         \
         * Get an object derived from setDataStorageBase to set values                              \
         */                                                                                         \
        return classtype< vartype >( &_##varname, ##__VA_ARGS__ );                                  \
    }

/*!
 * \brief Declares a setDataStorage getter function
 * \param varname: The name of the variable
 * \param classtype: The type of the class to be defined
 * \param vartype: The ctype of the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_SETDATASTORAGE_GETTER(varname,classtype,vartype,...) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_SETDATASTORAGE_GETTER(get_setDataStorage_##varname,varname,classtype,vartype, ##__VA_ARGS__)

/*!
 * \brief Declares a named getter function
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param uncall:  The function that is called if the variable is undefined.
 *     This should set the variable or throw an error.
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(getname,varname,vartype,uncall) \
    const vartype *getname( ){                                                \
        /*!                                                                   \
         * Get the value of varname                                           \
         */                                                                   \
        if(!_##varname.first){                                                \
            TARDIGRADE_ERROR_TOOLS_CATCH( uncall( ) )                         \
        }                                                                     \
        return &_##varname.second;                                            \
    }

/*!
 * \brief Declares a getter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param uncall:  The function that is called if the variable is undefined.
 *     This should set the variable or throw an error
 */
#define TARDIGRADE_HYDRA_DECLARE_GETTER(varname,vartype,uncall)                \
    TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(get_##varname,varname,vartype,uncall)

/*!
 * \brief Declares a named setter function
 * \param setname: The name of the setter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER(setname,varname,vartype,setfun) \
    void setname(const vartype &varname ){                                    \
    /*!                                                                       \
     * Set the storage variable varname of type vartype                       \
     *                                                                        \
     * \param &varname: The value of varname                                  \
     */                                                                       \
        TARDIGRADE_ERROR_TOOLS_CATCH( setfun( varname, _##varname ) )         \
    }

/*!
 * \brief Declares a setter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_SETTER(varname,vartype,setfun)                      \
    TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER( set_##varname, varname, vartype, setfun )

/*!
 * \brief Declare a dataStorage variable and the associated named setters and getters
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter function
 * \param getname: The name of the getter function
 * \param getsetdatastoragename: The name of the getter function for the setDataStorage object
 * \param varname: The name of the variable
 * \param classtype: The type of the class for the setDataStorage object
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,getsetdatastoragename,varname,classtype,vartype,setfun,uncall,...) \
    private: dataStorage<vartype> _##varname;                       \
    public: TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER(setname,varname,vartype,setfun) \
    public: TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(getname,varname,vartype,uncall) \
    public: TARDIGRADE_HYDRA_DECLARE_NAMED_SETDATASTORAGE_GETTER(getsetdatastoragename,varname,classtype,vartype, ##__VA_ARGS__) \
    context:

/*!
 * \brief Declare a dataStorage variable and the associated setters and getters
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_STORAGE(context,varname,vartype,setfun,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,set_##varname,get_##varname,get_setDataStorage_##varname,varname,tardigradeHydra::setDataStorageBase,vartype,setfun,uncall)

/*!
 * \brief Declare a named dataStorage variable that uses setIterationData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(context,setname,getname,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,get_setDataStorage_##varname,varname,setDataStorageIteration,vartype,setIterationData,uncall,this)

/*!
 * \brief Declare a dataStorage variable that uses setIterationData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,set_##varname,get_##varname,get_setDataStorage_##varname,varname,setDataStorageIteration,vartype,setIterationData,uncall,this)

/*!
 * \brief Declare a dataStorage variable that uses setNLStepData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NLSTEP_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,set_##varname,get_##varname,get_setDataStorage_##varname,varname,setDataStorageNLStep,vartype,setNLStepData,uncall,this)

/*!
 * \brief Declare a named dataStorage variable that uses setPreviousData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(context,setname,getname,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,get_setDataStorage_##varname,varname,setDataStoragePrevious,vartype,setPreviousData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setPreviousData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,set_##varname,get_##varname,get_setDataStorage_##varname,varname,setDataStoragePrevious,vartype,setPreviousData,uncall)

/*!
 * \brief Declare a named dataStorage variable that uses setConstantData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_CONSTANT_STORAGE(context,setname,getname,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,get_setDataSTorage_##varname,varname,setDataStorageConstant,vartype,setConstantData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setContantData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,set_##varname,get_##varname,get_setDataStorage_##varname,varname,setDataStorageConstant,vartype,setConstantData,uncall)

namespace tardigradeHydra{

//    namespace core{

        /*!
         * Base class for data objects which defines the clear command
         */
        class dataBase{

            public:

                virtual void clear( ){
                    /*!
                     * The function to erase the current values stored
                     */

                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "clear not implemented!" ) );

                }

        };

        /*!
         * Custom data storage object that allows for smart storage of objects
         */
        template < typename T >
        class dataStorage : public dataBase{

            public:

                bool first = false; //!< The flag for whether the data has been stored

                T second; //!< The stored data

                dataStorage( ){ };

                /*!
                 * Constructor for a data-storage object setting first and second
                 * 
                 * \param &_first: The flag for whether the data storage has been set
                 * \param &_second: The data contained by the object
                 */
                dataStorage( const bool &_first, const T &_second ) : first( _first ), second( _second ) { }

                virtual void clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and clearing second
                     */

                    first = false;

                    second.clear( );

                }

                virtual void zero( ){
                    /*!
                     * The function to set the current values to zero
                     */

                    std::fill( std::begin( second ), std::end( second ), 0 );

                }

                virtual void zero( const unsigned int size ){
                    /*!
                     * The function to resize second and set the current values to zero
                     *
                     * \param size: The size to resize the vector to
                     */

                    second.resize( size );

                    zero( );

                }

        };

        template <>
        inline void dataStorage< int >::clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and second to zero
                     */

            first = false;

            second = 0;

        }

        template <>
        inline void dataStorage< int >::zero( ){
                    /*!
                     * The function to set the value to zero
                     */

            second = 0;

        }

        template <>
        inline void dataStorage< int >::zero( const unsigned int size ){
                    /*!
                     * The function to set the value to zero
                     */

            throw std::runtime_error( "A scalar value cannot have a size!");

        }

        template <>
        inline void dataStorage< unsigned int >::clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and second to zero
                     */

            first = false;

            second = 0;

        }

        template <>
        inline void dataStorage< unsigned int >::zero( ){
                    /*!
                     * The function to set the value to zero
                     */

            second = 0;

        }

        template <>
        inline void dataStorage< unsigned int >::zero( const unsigned int size ){
                    /*!
                     * The function to set the value to zero
                     */

            throw std::runtime_error( "A scalar value cannot have a size!");

        }

        /*!
         * A custom object that handles setting dataStorage objects in place
         * When the destructor is run the dataStorage is assumed to be set.
         */
        template< typename T >
        class setDataStorageBase {

            public:

                //! Default constructor which points to no data
                setDataStorageBase( ) : _ds( NULL ){ }

                //! Constructor that reads in a data storage object and sets the value
                setDataStorageBase( dataStorage< T > *ds ) : _ds( ds ){

                    if ( _ds ){

                        value = &_ds->second;

                    }

                }

                //! Destructor which indicates that the data storage object has been set
                ~setDataStorageBase( ){

                    if ( _ds ){

                        _ds->first = true;

                    }

                }

                T * value; //!< A pointer to the value of the data storage object

                //! Get the starting iterator of the data storage object
                auto begin( ){ return std::begin( *value ); }

                //! Get the stopping iterator of the data storage object
                auto end( ){ return std::end( *value ); }

                //! Initialize the data storage object to zero
                void zero( ){ _ds->zero( ); }

                void zero( const unsigned int size ){
                    /*!
                     * Resize (if possible) the data storage object and set to zero
                     * 
                     * \param size: The size of the data storage object
                     */

                    _ds->zero( size );

                }

                template< typename X, unsigned int R, unsigned int C >
                Eigen::Map< Eigen::Matrix< X, R, C, Eigen::RowMajor > > zeroMap( ){
                    /*!
                     * Create a zeroed Eigen::Map of the quantity
                     */
                    zero( R * C );
                    return Eigen::Map< Eigen::Matrix< X, R, C, Eigen::RowMajor > > ( value->data( ), R, C );
                }

                template< typename X, unsigned int R >
                Eigen::Map< Eigen::Matrix< X, R, -1, Eigen::RowMajor > > zeroMap( const unsigned int C ){
                    /*!
                     * Create a zeroed Eigen::Map of the quantity with dynamic columns
                     * 
                     * \param C: The number of columns in the matrix
                     */
                    zero( R * C );
                    return Eigen::Map< Eigen::Matrix< X, R, -1, Eigen::RowMajor > > ( value->data( ), R, C );
                }

                template< typename X, unsigned int R >
                Eigen::Map< Eigen::Vector< X, R > > zeroMap( ){
                    /*!
                     * Create a zeroed Eigen::Map of the quantity with dynamic columns as a vector
                     */
                    zero( R );
                    return Eigen::Map< Eigen::Vector< X, R > > ( value->data( ), R );
                }

                template< typename X >
                Eigen::Map< Eigen::Vector< X, -1 > > zeroMap( const unsigned int R ){
                    /*!
                     * Create a zeroed Eigen::Map of the quantity with dynamic columns as a vector
                     * 
                     * \param R: The number of rows in the vector
                     */
                    zero( R );
                    return Eigen::Map< Eigen::Vector< X, -1 > > ( value->data( ), R );
                }

                template< typename X >
                Eigen::Map< Eigen::Matrix< X, -1, -1, Eigen::RowMajor > > zeroMap( const unsigned int R, const unsigned int C ){
                    /*!
                     * Create a zeroed Eigen::Map of the quantity with dynamic columns
                     * 
                     * \param R: The number of rows in the matrix
                     * \param C: The number of columns in the matrix
                     */
                    zero( R * C );
                    return Eigen::Map< Eigen::Matrix< X, -1, -1, Eigen::RowMajor > > ( value->data( ), R, C );
                }

          protected:

              dataStorage< T > *_ds; //!< Pointer to the data for the data storage object

        };

        template< typename container, typename T>
        class setDataStorageIterationBase : public setDataStorageBase< T >{
          public:

              setDataStorageIterationBase( dataStorage< T > *ds, container * rp ) : setDataStorageBase< T >( ds ), _rp( rp ){
                  /*!
                   * Create a data storage object that will be reset at each new iteration
                   * 
                   * \param *ds: The data storage object
                   * \param *rp: The class that contains the data storage object. Must have a function addIterationData
                   *     wich accepts a pointer to the data storage object
                   */ }

              //! Destructor object that adds the data storage object to the iteration data list
              ~setDataStorageIterationBase( ){

                  _rp->addIterationData( this->_ds );

              }

          protected:

              container *_rp; //!< The containing hydraBase class

        };

        // CoreDefinitions dataStorage specifications
        template <>
        inline void dataStorage< floatType >::clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and second to zero
                     */

            first = false;

            second = 0;

        }

        template <>
        inline void dataStorage< floatType >::zero( ){
                    /*!
                     * The function to set the value to zero
                     */

            second = 0;

        }

        template <>
        inline void dataStorage< floatType >::zero( const unsigned int size ){
                    /*!
                     * The function to set the value to zero
                     */

            throw std::runtime_error( "A scalar value cannot have a size!");

        }

//    }

}

#include "tardigrade_SetDataStorage.cpp"

#endif
