/**
  ******************************************************************************
  * \file tardigrade_hydra.h
  ******************************************************************************
  * A C++ library for constructing finite deformation constitutive models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_H
#define TARDIGRADE_HYDRA_H

#include<sstream>
#include<functional>

#include<tardigrade_error_tools.h>
//!We will use the functions that depend on Eigen
#define USE_EIGEN
#include<tardigrade_vector_tools.h>

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
#include<libxsmm.h>
#endif

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

    // forward class definitions
    class hydraBase;

    class residualBase;

    template<typename T>
    class setDataStorageBase;

    template<typename T>
    class setDataStorageIteration;

    template<typename T>
    class setDataStoragePrevious;

    template<typename T>
    class setDataStorageConstant;

    namespace unit_test{
        class hydraBaseTester;
    }

    constexpr const char* str_end(const char *str) {
        /*! Recursively search string for last character
         * \param *str: pointer to string START of UNIX path like string
         * \return *str: pointer to last character in string
         */
        return *str ? str_end(str + 1) : str;
    }
    constexpr bool str_slant(const char *str) {
        /*! Recursively search string for leftmost UNIX path separator from the left
         * \param *str: pointer to string START of UNIX path like string
         * \return bool: True if string contains UNIX path separator. Else false.
         */
        return *str == '/' ? true : (*str ? str_slant(str + 1) : false);
    }
    constexpr const char* r_slant(const char* str) {
        /*! Recursively search string for rightmost UNIX path separator from the right
         * \param *str: pointer to string END of UNIX path like string
         * \return *str: pointer to start of base name
         */
        return *str == '/' ? (str + 1) : r_slant(str - 1);
    }
    constexpr const char* file_name(const char* str) {
        /*! Return the current file name with extension at compile time
         * \param *str: pointer to string START of UNIX path like string
         * \return str: file base name
         */
        return str_slant(str) ? r_slant(str_end(str)) : str;
    }
    //Return filename for constructing debugging messages
    //https://stackoverflow.com/questions/31050113/how-to-extract-the-source-filename-without-path-and-suffix-at-compile-time
    const std::string __BASENAME__ = file_name(__FILE__); //!< The base filename which will be parsed
    const std::string __FILENAME__ = __BASENAME__.substr(0, __BASENAME__.find_last_of(".")); //!< The parsed filename for error handling

    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    //Define tensors of known size
    typedef std::vector< floatType > dimVector; //!< Dimension vector
    typedef std::vector< floatType > secondOrderTensor; //!< Second order tensors
    typedef std::vector< floatType > thirdOrderTensor; //!< Third order tensors
    typedef std::vector< floatType > fourthOrderTensor; //!< Fourth order tensors

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
    typedef libxsmm_mmfunction<floatType> kernel_type; //!< The libxsmm kernel type
#endif

    typedef void ( hydraBase::*hydraBaseFxn )( ); //!< Typedef for passing pointers to hydraBase functions

    template<typename T, int R, int C>
    Eigen::Map< Eigen::Matrix< T, R, C, Eigen::RowMajor > > getFixedSizeMatrixMap( T *p ){
        /*!
         * Get a matrix of type T with a fixed size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */
        return Eigen::Map< Eigen::Matrix< T, R, C, Eigen::RowMajor > >( p, R, C );
    }

    template<typename T, int R>
    Eigen::Map< Eigen::Matrix< T, R, -1, Eigen::RowMajor > > getDynamicColumnSizeMatrixMap( T *p, int C ){
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param C: The number of columns in the map
         */
        return Eigen::Map< Eigen::Matrix< T, R, -1, Eigen::RowMajor > >( p, R, C );
    }

    template<typename T>
    Eigen::Map< Eigen::Matrix< T, -1, -1, Eigen::RowMajor > > getDynamicSizeMatrixMap( T *p, int R, int C ){
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         * \param C: The number of columns in the map
         */
        return Eigen::Map< Eigen::Matrix< T, -1, -1, Eigen::RowMajor > >( p, R, C );
    }

    template<typename T, int R, int C>
    Eigen::Map< const Eigen::Matrix< T, R, C, Eigen::RowMajor > > getFixedSizeMatrixMap( const T *p ){
        /*!
         * Get a matrix of type T with a fixed size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */
        return Eigen::Map< const Eigen::Matrix< T, R, C, Eigen::RowMajor > >( p, R, C );
    }

    template<typename T, int R>
    Eigen::Map< const Eigen::Matrix< T, R, -1, Eigen::RowMajor > > getDynamicColumnSizeMatrixMap( const T *p, int C ){
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param C: The number of columns in the map
         */
        return Eigen::Map< const Eigen::Matrix< T, R, -1, Eigen::RowMajor > >( p, R, C );
    }

    template<typename T>
    Eigen::Map< const Eigen::Matrix< T, -1, -1, Eigen::RowMajor > > getDynamicSizeMatrixMap( const T *p, int R, int C ){
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         * \param C: The number of columns in the map
         */
        return Eigen::Map< const Eigen::Matrix< T, -1, -1, Eigen::RowMajor > >( p, R, C );
    }

    template<typename T, int R>
    Eigen::Map< Eigen::Vector< T, R > > getFixedSizeVectorMap( T *p ){
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */

        return Eigen::Map< Eigen::Vector< T, R > >( p, R );
    }

    template<typename T, int R>
    Eigen::Map< const Eigen::Vector< T, R > > getFixedSizeVectorMap( const T *p ){
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */

        return Eigen::Map< const Eigen::Vector< T, R > >( p, R );
    }

    template<typename T>
    Eigen::Map< Eigen::Vector< T, -1 > > getDynamicSizeVectorMap( T *p, int R ){
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         */

        return Eigen::Map< Eigen::Vector< T, -1 > >( p, R );
    }

    template<typename T>
    Eigen::Map< const Eigen::Vector< T, -1 > > getDynamicSizeVectorMap( const T *p, int R ){
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         */

        return Eigen::Map< const Eigen::Vector< T, -1 > >( p, R );
    }

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
    inline void dataStorage< std::vector< residualBase* > >::zero( ){
                /*!
                 * The function to set the value to zero
                 */

        throw std::runtime_error( "Zeroing the residualBase pointer vector is not allowed" );

    }

    template <>
    inline void dataStorage< std::vector< residualBase* > >::zero( const unsigned int size ){
                /*!
                 * The function to set the value to zero
                 */

        throw std::runtime_error( "Zeroing the residualBase pointer vector is not allowed" );

    }

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

    /*!
     * A custom error for use with failures in convergence of the solver.
     */
    class convergence_error : public std::exception{

        private:
            std::string message_; //!< The output message

        public:

            //! Constructor
            explicit convergence_error( const std::string& message ) : message_( message ) { }

            const char *what( ) const noexcept override {
                /*!
                 * Output the message
                 */

                return message_.c_str( );

            }

    };

    /*!
     * A class to contain the residual computations associated with some part of a non-linear solve
     */
    class residualBase{

        public:

            /*!
             * Default residual
             */
            residualBase( ) : hydra( NULL ), _numEquations( 0 ){ };

            /*!
             * Main utilization constructor
             * 
             * \param *_hydra: A pointer to a hydraBase object
             * \param &_numEquations: The number of equations defined by the residual
             */
            residualBase( hydraBase *_hydra, const unsigned int &_numEquations ) : hydra( _hydra ), _numEquations( _numEquations ){ }

            /*!
             * Copy constructor
             * 
             * \param &r: The residual to be copied
             */
            residualBase( residualBase &r ) : hydra( r.hydra ), _numEquations( *r.getNumEquations( ) ), _numConstraints( *r.getNumConstraints( ) ){ }

            hydraBase* hydra; //!< The hydra class which owns the residualBase object

            // User defined setter functions

            virtual void setResidual( ){
                /*!
                 * The user-defined residual equation. Must have a size of numEquations
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The residual is not implemented" ) );

            }

            virtual void setJacobian( ){
                /*!
                 * The user-defined jacobian equation. Must have a size of numEquations x numUnknowns
                 * 
                 * The order of the unknowns are the cauchy stress, the configurations in order (minus the first one),
                 * and the state variables solved for in the non-linear solve.
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The jacobian is not implemented" ) );

            }

            virtual void setdRdF( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. the deformation gradient.
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The derivative of the residual w.r.t. the deformation gradient is not implemented" ) );
 
            }

            virtual void setdRdT( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. the temperature
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The derivative of the residual w.r.t. the temperature is not implemented" ) );
 
            }

            virtual void setdRdAdditionalDOF( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. the additional DOF
                 */

            }

            virtual void setAdditionalDerivatives( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. additional values
                 */

            }

            virtual void setStress( ){
                /*!
                 * Compute the current stress
                 * 
                 * Only needs to be defined for the first residual
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the stress is not implemented" ) );

            }

            virtual void setPreviousStress( ){
                /*!
                 * Compute the previous stress
                 * 
                 * Only needs to be defined for the first residual
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the previous stress is not implemented" ) );

            }

            virtual void setConstraints( ){
                /*!
                 * Compute the contraints
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the constraints is not implemented" ) );

            }

            virtual void setConstraintJacobians( ){
                /*!
                 * Compute the contraint Jacobians
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the constraint jacobians is not implemented" ) );

            }

            virtual void setCurrentAdditionalStateVariables( ){
                /*!
                 * Set the current additional state variables
                 * 
                 * Doesn't need to be defined
                 */

                setCurrentAdditionalStateVariables( floatVector( 0, 0 ) );

                return;

            }

            virtual void suggestInitialIterateValues( std::vector< unsigned int >   &indices,
                                                      std::vector< floatType > &values ){

                /*!
                 * Function which is called which allows the residual to suggest initial values for given
                 * configurations. This is called when the unknown vector is being initialized. If more than
                 * one residual attempts to set the initial vector the last residual will override all of the others.
                 *
                 * After the initial iterate has been suggested, the iteration data is cleared so that the residual
                 * starts the iteration in a clean state.
                 * 
                 * \param &indices: The indices of the unknown vector to set
                 * \param &values:  The values to be set in the unknown vector
                 */

                indices.clear( );
                values.clear( );

            }

            virtual void projectSuggestedX( std::vector< floatType > &trialX,
                                            const std::vector< floatType > &Xp ){
                /*!
                 * Project the suggested unknown vector to the allowable space
                 * 
                 * Called whenever hydra calls updateUnknownVector. It is assumed that the
                 * initial value as suggested by `residual::suggestInitialIterationValues` is
                 * in the allowable space.
                 * 
                 * \param &trialX: The trial value of X
                 * \param &Xp: The previously accepted value of X
                 */

            }

            void setUseProjection( const bool &value ){
                /*!
                 * Set whether to use the projection or not
                 * 
                 * \param &value: The value of the parameter
                 */

                _useProjection = value;

            }

            virtual bool checkRelaxedConvergence( ){
                /*!
                 * When performing a relaxed solve the residuals must return if they are converged or not.
                 * This function returns true or false (default true)
                 */

                return true;

            }

            virtual void modifyGlobalResidual( ){
                /*!
                 * Function that is called to modify the global residual.
                 *
                 * Called after all of the residuals are agglomerated into the whole.
                 *
                 * A mutable version of the global residual is accessable with hydra->getMutableResidual( )
                 */
            }

            virtual void modifyGlobalJacobian( ){
                /*!
                 * Function that is called to modify the global jacobian.
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global jacobian is accessable with hydra->getMutableJacobian( )
                 */
            }

            virtual void modifyGlobaldRdT( ){
                /*!
                 * Function that is called to modify the global derivative of the residual w.r.t. the temperature
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global dRdT is accessable with hydra->getMutabledRdT( )
                 */
            }

            virtual void modifyGlobaldRdF( ){
                /*!
                 * Function that is called to modify the global derivative of the residual w.r.t. the deformation
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global dRdF is accessable with hydra->getMutabledRdF( )
                 */
            }

            virtual void modifyGlobaldRdAdditionalDOF( ){
                /*!
                 * Function that is called to modify the global derivative of the residual w.r.t. the additional DOF
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global dRdAdditionalDOF is accessable with hydra->getMutabledRdAdditionalDOF( )
                 */
            }

            virtual void preNLSolve( ){
                /*!
                 * Function that is called prior to a nonlinear solve
                 */
            };

            virtual void postNLSolve( ){
                /*!
                 * Function that is called after a nonlinear solve
                 */
            };

            virtual void successfulNLStep( ){
                /*!
                 * Function that is called whenever a successful nonlinear step is taken
                 */
            };

            virtual void setupRelaxedStep( const unsigned int &relaxedStep );

            //! Get the flag for whether to use the projection or not
            const bool *getUseProjection( ){ return &_useProjection; }

            // Getter functions

            //! Get the number of equations the residual defined
            const unsigned int* getNumEquations( ){ return &_numEquations; }

            //! Get the number of constraints the residual defined
            const unsigned int* getNumConstraints( ){ return &_numConstraints; }

            void addIterationData( dataBase *data );

            void addNLStepData( dataBase *data );

            template<class T>
            void setIterationData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding iteration data
                 *
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

                addIterationData( &storage );

            }

            template<class T>
            void setNLStepData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding nonlinear step data
                 *
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

                addNLStepData( &storage );

            }

            template<class T>
            void setPreviousData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding previous data
                 * 
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

            }

            template<class T>
            void setConstantData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding constant data
                 * 
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

            }

            virtual void addParameterizationInfo( std::string &parameterization_info ){
                /*!
                 * Add parameterization information to the provided docstring
                 * 
                 * Information on the current parameters and their values should be added to the
                 * incoming string.
                 * 
                 * \param &parameterization_info: The parameterization information string. Append
                 *     this class' information
                 */

                parameterization_info += "NO INFORMATION DEFINED FOR CLASS\n";

            }

            virtual void updateAdditionalStateVariables( floatVector &additionalStateVariables ){
                /*!
                 * Update the additional state variable vector
                 *
                 * \param &additionalStateVariables: The additional state variable vector
                 */

            }

        protected:

            void setNumConstraints( const unsigned int numConstraints ){
                /*!
                 * Set the number of constraints for the solve
                 * 
                 * \param numConstraints: The number of constraints
                 */

                _numConstraints = numConstraints;

            }

            void setPenaltyIndices( const std::vector< unsigned int > &indices ){
                /*!
                 * Set the indices where the penalties are defined
                 * 
                 * \param &indices: The indices where the penalty is defined
                 */
                _penalty_indices = indices;
            }

            void unexpectedError( ){
                /*!
                 * Function to throw for an unexpected error. A user should never get here!
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "You shouldn't have gotten here. If you aren't developing the code then contact a developer with the stack trace." ) )

            }

            //! Class which defines data storage objects which are reset at each iteration
            template< typename T >
            class setDataStorageIteration : public setDataStorageBase< T > {

              public:

                  setDataStorageIteration( dataStorage< T > *ds, residualBase * rp ) : setDataStorageBase< T >( ds ), _rp( rp ){
                      /*!
                       * Create a data storage object that will be reset at each new iteration
                       * 
                       * \param *ds: The data storage object
                       * \param *rp: The residual class that contains the data storage object
                       */
                  }

                  ~setDataStorageIteration( ){
                      /*!
                       * Destructor that says the data storage object has been set
                       */

                      _rp->addIterationData( this->_ds );

                  }

              protected:

                  residualBase *_rp; //!< The containing residual class

            };

            //! Class which defines data storage objects which are reset at each nonlinear step
            template< typename T >
            class setDataStorageNLStep : public setDataStorageBase< T > {

              public:

                  setDataStorageNLStep( dataStorage< T > *ds, residualBase * rp ) : setDataStorageBase< T >( ds ), _rp( rp ){
                      /*!
                       * Create a data storage object that will be reset at each nonlinear step
                       * 
                       * \param *ds: The data storage object
                       * \param *rp: The residual class that contains the data storage object
                       */
                  }

                  ~setDataStorageNLStep( ){
                      /*!
                       * Destructor that says the data storage object has been set
                       */

                      _rp->addNLStepData( this->_ds );

                  }

              protected:

                  residualBase *_rp; //!< The containing residual class

            };

            //! Class which defines data storage objects for values defined at the previous timestep
            template< typename T >
            class setDataStoragePrevious : public setDataStorageBase< T > {

                public:

                    setDataStoragePrevious( dataStorage< T > *ds ) : setDataStorageBase< T >( ds ){
                        /*!
                         * Constructor for data storage objects for temporally previous objects
                         * 
                         * \param *ds: The data storage object to modify
                         */
                    }

            };

            /*!
             * Class that is a constant data storage object
             */
            template< typename T >
            class setDataStorageConstant : public setDataStorageBase< T > {

                public:

                    setDataStorageConstant( dataStorage< T > *ds ) : setDataStorageBase< T >( ds ){
                        /*!
                         * Constructor for constant data storage objects
                         * 
                         * \param *ds: The data storage object
                         */

                    }

            };

        private:

            unsigned int _numEquations; //!< The number of residual equations

            unsigned int _numConstraints = 0; //!< The number of constraint equations

            bool _useProjection = false; //!< Flag for whether to use the projection or not

            std::vector< unsigned int > _penalty_indices; //!< The indices of the variables which should be penalized for negative values

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setResidual,                         getResidual,                        residual,                        floatVector, setResidual )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setJacobian,                         getJacobian,                        jacobian,                        floatVector, setJacobian )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdF,                             getdRdF,                            dRdF,                            floatVector, setdRdF )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdT,                             getdRdT,                            dRdT,                            floatVector, setdRdT )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdAdditionalDOF,                 getdRdAdditionalDOF,                dRdAdditionalDOF,                floatVector, setdRdAdditionalDOF )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setAdditionalDerivatives,            getAdditionalDerivatives,           additionalDerivatives,           floatVector, setAdditionalDerivatives )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setStress,                           getStress,                          stress,                          floatVector, setStress )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setPreviousStress,                   getPreviousStress,                  previousStress,                  floatVector, setPreviousStress )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraints,                      getConstraints,                     constraints,                     floatVector, setConstraints )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraintJacobians,              getConstraintJacobians,             constraintJacobians,             floatVector, setConstraintJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setCurrentAdditionalStateVariables,  getCurrentAdditionalStateVariables, currentAdditionalStateVariables, floatVector, setCurrentAdditionalStateVariables )

    };

    /*!
     * hydraBase: A base class which can be used to construct finite deformation material models.
     * 
     * The hydra class seeks to provide utilities for the construction of finite deformation constitutive models
     * more rapidly than would be possible previously. The user can define as many different configurations as desired
     * and provide a calculation of the Cauchy stress.
     * 
     * A non-linear problem which is of the size ( dimension**2 * num_configurations + num_ISVs ) will be solved.
     */
    class hydraBase{

        public:

            // Constructors
            //! Default constructor for hydraBase
            hydraBase( ) : _configuration_unknown_count( 0 ){ }

            //! Main constructor for objects of type hydraBase. Sets all quantities required for most solves.
            hydraBase( const floatType &time, const floatType &deltaTime,
                       const floatType &temperature, const floatType &previousTemperature,
                       const secondOrderTensor &deformationGradient, const secondOrderTensor &previousDeformationGradient,
                       const floatVector &additionalDOF, const floatVector &previousAdditionalDOF,
                       const floatVector &previousStateVariables, const floatVector &parameters,
                       const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                       const unsigned int dimension=3, const unsigned int configuration_unknown_count=9,
                       const floatType tolr=1e-9, const floatType tola=1e-9,
                       const unsigned int maxIterations=20, const unsigned int maxLSIterations=5, const floatType lsAlpha=1e-4,
                       const bool use_preconditioner=false, const unsigned int preconditioner_type=0 );

            // User defined functions

            // Setter functions

            // Getter functions
            //! Get a reference to the number of unknowns in each configuration
            const unsigned int* getConfigurationUnknownCount( ){ return &_configuration_unknown_count; }

            //! Get a reference to the number of components of the stress
            const unsigned int* getStressSize( ){ return &_stress_size; }

            //! Get a reference to the current time
            const floatType* getTime( ){ return getScaledTime( ); }

            //! Get a reference to the change in time
            const floatType* getDeltaTime( ){ return getScaledDeltaTime( ); }

            //! Get a reference to the current temperature
            const floatType* getTemperature( ){ return getScaledTemperature( ); };

            //! Get a reference to the previous temperature
            const floatType* getPreviousTemperature( ){ return &_previousTemperature; };

            //! Get a reference to the deformation gradient
            const secondOrderTensor* getDeformationGradient( ){ return getScaledDeformationGradient( ); }

            //! Get a reference to the previous deformation gradient
            const secondOrderTensor* getPreviousDeformationGradient( ){ return &_previousDeformationGradient; }

            //! Get a reference to the additional degrees of freedom
            const floatVector* getAdditionalDOF( ){ return getScaledAdditionalDOF( ); }

            //! Get a reference to the previous additional degrees of freedom
            const floatVector* getPreviousAdditionalDOF( ){ return &_previousAdditionalDOF; }

            //! Get a reference to the previous values of the state variables
            const floatVector* getPreviousStateVariables( ){ return &_previousStateVariables; }

            //! Get a reference to the model parameters
            const floatVector* getParameters( ){ return &_parameters; }

            //! Get a reference to the number of configurations
            const unsigned int* getNumConfigurations( ){ return &_numConfigurations; }

            //! Get a reference to the number of state variables involved in the non-linear solve
            const unsigned int* getNumNonLinearSolveStateVariables( ){ return &_numNonLinearSolveStateVariables; }

            //! Get the number of terms in the unknown vector
            virtual const unsigned int getNumUnknowns( ){ return ( *getNumConfigurations( ) ) * ( *getConfigurationUnknownCount( ) ) + *getNumNonLinearSolveStateVariables( ); }

            //! Get the number of additional degrees of freedom
            virtual const unsigned int getNumAdditionalDOF( ){ return getAdditionalDOF( )->size( ); }

            //! Get the value of the number of constraint equations
            virtual const unsigned int getNumConstraints( ){

                unsigned int value = 0;

                for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->end( ); v++ ){

                    value += *( *v )->getNumConstraints( );

                }

                return value;

            }

            //! Get the current residual index
            const unsigned int getCurrentResidualIndex( ){

                TARDIGRADE_ERROR_TOOLS_CHECK( currentResidualIndexMeaningful( ), "The current residual index isn't meaningful" );
                return _current_residual_index;

            }

            const unsigned int getCurrentResidualOffset( ){
                /*!
                 * Get the offset of the current residual
                 */
                unsigned int offset = 0;
                for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->begin( ) + getCurrentResidualIndex( ); ++v ){
                    offset += *( *v )->getNumEquations( );
                }
                return offset;
            }

            //! Get a reference to the dimension
            constexpr unsigned int getDimension( ){ return _dimension; }

            //! Get a reference to a second order tensor's dimension
            constexpr unsigned int getSOTDimension( ){ return _dimension * _dimension; }

            //! Get a reference to a third order tensor's dimension
            constexpr unsigned int getTOTDimension( ){ return _dimension * _dimension * _dimension; }

            //! Get a reference to a fourth order tensor's dimension
            constexpr unsigned int getFOTDimension( ){ return _dimension * _dimension * _dimension * _dimension; }

            //! Get a reference to the relative tolerance
            const floatType* getRelativeTolerance( ){ return &_tolr; }

            //! Get a reference to the absolute tolerance
            const floatType* getAbsoluteTolerance( ){ return &_tola; }

            //! Get a reference to the line-search alpha
            const floatType* getLSAlpha( ){ return &_lsAlpha; }

            //! Get a reference to whether to use a preconditioner
            const bool* getUsePreconditioner( ){ return &_use_preconditioner; }

            //! Get a reference to the preconditioner type
            const unsigned int* getPreconditionerType( ){ return &_preconditioner_type; }

            //! Get whether the preconditioner is diagonal or not
            const bool* getPreconditionerIsDiagonal( ){ return &_preconditioner_is_diagonal; }

            //!< Get a reference to the gradient descent sigma parameter
            const floatType* getGradientSigma( ){ return &_gradientSigma; }

            //!< Get a reference to the gradient descent beta parameter
            const floatType* getGradientBeta( ){ return &_gradientBeta; }

            //!< Get a reference to the max allowable number of gradient iterations
            const unsigned int* getMaxGradientIterations( ){ return &_maxGradientIterations; }

            //!< Get a reference to the gradient descent rho parameter
            const floatType* getGradientRho( ){ return &_gradientRho; }

            //!< Get a reference to the gradient descent p parameter
            const floatType* getGradientP( ){ return &_gradientP; }

            //!< Get a reference to the Levenberg-Marquardt mu parameter
            const floatType* getLMMu( ){ return &_lm_mu; }

            //!< Get a reference to the current value of mu_k
            const floatType* getMuk( ){ return &_mu_k; }

            //!< Get a reference to whether the Newton step should be a LevenbergMarquardt step
            const bool* getUseLevenbergMarquardt( ){ return &_use_LM_step; }

            //!< Get a reference to whether the Newton step should be a relaxed solve
            const bool* getUseRelaxedSolve( ){ return &_use_relaxed_solve; }

            //!< Get a reference to whether Gradient descent is allowed
            const bool* getUseGradientDescent( ){ return &_use_gradient_descent; }

            //!< Get a reference to the flag for whether to throw an error if the LHS matrix is rank-deficient
            const bool* getRankDeficientError( ){ return &_rank_deficient_error; }

            //!< Set the gradient descent sigma parameter
            void setGradientSigma( const floatType &value ){
               /*!
                * Set the value of the sigma parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */

                _gradientSigma = value;

            }

            //!< Set the gradient descent beta parameter
            void setGradientBeta( const floatType &value ){
               /*!
                * Set the value of the beta parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */
 
                _gradientBeta = value;

            }

            //!< Set the max allowable number of gradient iterations
            void setMaxGradientIterations( const unsigned int &value ){
               /*!
                * Set the value of the maximum number of iterations for gradient descent steps
                *
                * \param &value: The value of the parameter
                */

                _maxGradientIterations = value;

            }

            //!< Set the gradient descent rho parameter
            void setGradientRho( const floatType &value ){
               /*!
                * Set the value of the rho parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */
 
                _gradientRho = value;

            }

            //!< Set the gradient descent p parameter
            void setGradientP( const floatType &value ){
               /*!
                * Set the value of the p parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */
 
                _gradientP = value;

            }

            //!< Set the Levenberg-Marquardt mu parameter
            void setLMMu( const floatType &value ){
               /*!
                * Set the value of the mu parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _lm_mu = value;

            }

            //!< Set the Levenberg-Marquardt mu_k
            void setMuk( const floatType &value ){
               /*!
                * Set the value of the mu_k parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _mu_k = value;

            }

            //!< Set whether to attempt a Levenberg-Marquardt step
            void setUseLevenbergMarquardt( const bool &value ){
                /*!
                 * Set whether to attempt a Levenberg-Marquardt step
                 * 
                 * \param &value: The value of the parameter
                 */

                setUseGradientDescent( value );

                _use_LM_step = value;

            }

            //!< Set whether to attempt a Relaxed-solve
            void setUseRelaxedSolve( const bool &value ){
                /*!
                 * Set whether to attempt a relaxed solve
                 * 
                 * \param &value: The value of the parameter
                 */

                _use_relaxed_solve = value;

            }

            //!< Set whether to use gradient descent
            void setUseGradientDescent( const bool &value ){
                /*!
                 * Set whether to attempt a gradient descent step
                 * 
                 * \param &value: The value of the parameter
                 */

                _use_gradient_descent = value;

            }

            //!< Set whether rank deficiency is a reason to throw an error
            void setRankDeficientError( const bool &value ){
                /*!
                 * Set whether a rank-deficient LHS will cause an error
                 * 
                 * \param &value: The value of the parameter
                 */

                _rank_deficient_error = value;

            }

            secondOrderTensor getSubConfiguration( const floatVector &configurations, const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getSubConfigurationJacobian( const floatVector &configurations, const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getPrecedingConfiguration( const unsigned int &index );

            secondOrderTensor getFollowingConfiguration( const unsigned int &index );

            secondOrderTensor getConfiguration( const unsigned int &index );

            secondOrderTensor getPreviousSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getPreviousPrecedingConfiguration( const unsigned int &index );

            secondOrderTensor getPreviousFollowingConfiguration( const unsigned int &index );

            secondOrderTensor getPreviousConfiguration( const unsigned int &index );

            floatVector getSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPrecedingConfigurationJacobian( const unsigned int &index );

            floatVector getFollowingConfigurationJacobian( const unsigned int &index );

            floatVector getPreviousSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPreviousPrecedingConfigurationJacobian( const unsigned int &index );

            floatVector getPreviousFollowingConfigurationJacobian( const unsigned int &index );

            const floatType* getLSResidualNorm( );

            virtual void setResidualClasses( );

            void setResidualClasses( std::vector< residualBase* > &residualClasses );

            std::vector< residualBase* >* getResidualClasses( );

            const floatVector* getResidual( );

            const floatVector* getFlatJacobian( );

            virtual const floatVector* getFlatNonlinearLHS( );

            virtual const floatVector* getNonlinearRHS( );

            const floatVector* getFlatPreconditioner( );

            floatMatrix getJacobian( );

            const floatVector* getFlatdRdF( );

            floatMatrix getdRdF( );

            const floatVector* getdRdT( );

            const floatVector* getFlatdRdAdditionalDOF( );

            floatMatrix getdRdAdditionalDOF( );

            const floatVector* getFlatAdditionalDerivatives( );

            floatMatrix getAdditionalDerivatives( );

            const floatVector* getUnknownVector( );

            const floatVector* getTolerance( );

            virtual bool checkConvergence( );

            virtual bool checkLSConvergence( );

            virtual bool checkGradientConvergence( const floatVector &X0 );

            virtual bool checkDescentDirection( const floatVector &dx );

            const floatVector* getStress( );

            const floatVector* getPreviousStress( );

            virtual void evaluate( const bool &use_subcycler = false );

            virtual void computeTangents( );

            virtual void computedXdAdditionalDOF( );

            const floatVector *getFlatdXdF( );

            const floatVector *getFlatdXdT( );

            const floatVector *getFlatdXdAdditionalDOF( );

            //! Return the flag which indicates whether hydra should initialize the unknown vector
            const bool *getInitializeUnknownVector( ){ return &_initializeUnknownVector; }

            //! Add data to the vector of values which will be cleared after each iteration
            void addIterationData( dataBase *data ){ _iterationData.push_back( data ); }

            //! Add data to the vector of values which will be cleared after each non-linear step
            void addNLStepData( dataBase *data ){ _nlStepData.push_back( data ); }

            unsigned int getNumNewton( ){ /*! Get the number of newton steps performed */  return _NUM_NEWTON; }

            unsigned int getNumLS( ){ /*! Get the number of line search steps performed */ return _NUM_LS; }

            unsigned int getNumGrad( ){ /*! Get the number of gradient descent steps performed */  return _NUM_GRAD; }

            const bool *getUseSQPSolver( ){ /*! Return a flag for whether to use the SQP solver */ return &_useSQPSolver; }

            const void setMaxRelaxedIterations( const unsigned int &value ){
                /*!
                 * Set the maximum allowable number of relaxed iterations
                 * 
                 * \param &value: The number of relaxed iterations
                 */

                _maxRelaxedIterations = value;

            }

            const unsigned int *getMaxRelaxedIterations( ){
                /*!
                 * Get the maximum number of relaxed iterations
                 */

                return &_maxRelaxedIterations;

            }

            void setFailureVerbosityLevel( const unsigned int &value ){
                /*!
                 * Set the verbosity level for failures
                 * 
                 * \param &value: The verbosity level of the failure (defaults to zero)
                 */

                _failure_verbosity_level = value;

            }

            //! Set the failure output string to use scientific notation
            void setFailureOutputScientific( ){ _failure_output << std::scientific; }

            //! Get the verbosity level for failure outputs
            const unsigned int *getFailureVerbosityLevel( ){ return &_failure_verbosity_level; }

            //! Add a string to the failure output string
            void addToFailureOutput( const std::string &additional ){ _failure_output << additional; }

            //! Add a floating point vector to the output string
            void addToFailureOutput( const floatVector &value ){ for ( auto v = value.begin( ); v != value.end( ); v++ ){ _failure_output << *v << ", "; } _failure_output << "\n"; }

            //! Add a floating point value to the output string
            void addToFailureOutput( const floatType &value ){ _failure_output << value; _failure_output << "\n"; }

            //! Get the failure output string
            const std::string getFailureOutput( ){ return _failure_output.str( ); }

            //! Get a scale factor for the deformation
            const floatType *getScaleFactor( ){ return &_scale_factor; }

            //! Set the value of the scale factor. Will automatically re-calculate the deformation and trial stresses
            void setScaleFactor( const floatType &value );

            const floatType   *getScaledTime( ){ /*! Get the value of the scaled current time */ return &_scaled_time; }

            const floatType   *getScaledDeltaTime( ){ /*! Get the value of the scaled changed in time */ return &_scaled_deltaTime; }

            const floatType   *getScaledTemperature( ){ /*! Get the value of the scaled current temperature */ return &_scaled_temperature; }

            const floatVector *getScaledDeformationGradient( ){ /*! Get the value of the scaled current deformation gradient */ return &_scaled_deformationGradient; }

            const floatVector *getScaledAdditionalDOF( ){ /*! Get the value of the scaled current additional DOF */ return &_scaled_additionalDOF; }

            const floatType *getCutbackFactor( ){ /*! Get the value of the cutback factor */ return &_cutback_factor; }

            const unsigned int *getNumGoodControl( ){ /*! Get the number of good iterations we need to have before increasing the timestep */ return &_num_good_control; }

            const floatType *getGrowthFactor( ){ /*! Get the growth factor for the timestep increase */ return &_growth_factor; }

            const floatType *getMinDS( ){ /*! Get the minimum allowable ratio of the total timestep to the cutback timestep */ return &_minDS; }

            void setCutbackFactor( const floatType &value ){ /*! Get the current value of the cutback factor. \param &value: The value of the cutback */  _cutback_factor = value; }

            void setNumGoodControl( const unsigned int &value ){ /*! Set the number of good iterations that need to happen before the timestep increases. \param &value: The value of the number of good iterations prior to increasing the relative timestep */  _num_good_control = value; }

            void setGrowthFactor( const floatType &value ){ /*! Set the relative growth factor for the local timestep increase \param &value: The new value */ _growth_factor = value; }

            void setMinDS( const floatType &value ){ /*! Set the minimum value of the relative cutback timestep \param &value: The new value */  _minDS = value; }

            const bool allowStepGrowth( const unsigned int &num_good );

            floatVector *getMutableResidual( ){
                /*! Get a reference to the full residual that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobalResidual methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_residual ){
                    return &_residual.second;
                }

                return NULL;

            }

            floatVector *getMutableJacobian( ){
                /*! Get a reference to the full jacobian that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobalJacobian methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_jacobian ){
                    return &_jacobian.second;
                }

                return NULL;

            }

            floatVector *getMutabledRdT( ){
                /*! Get a reference to the full dRdT that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobaldRdT methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_dRdT ){
                    return &_dRdT.second;
                }

                return NULL;

            }

            floatVector *getMutabledRdF( ){
                /*! Get a reference to the full dRdF that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobaldRdF methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_dRdF ){
                    return &_dRdF.second;
                }

                return NULL;

            }

            floatVector *getMutabledRdAdditionalDOF( ){
                /*! Get a reference to the full dRdAdditionalDOF that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobaldRdAdditionalDOF methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_dRdAdditionalDOF ){
                    return &_dRdAdditionalDOF.second;
                }

                return NULL;

            }

            const bool currentResidualIndexMeaningful( ){
                /*!
                 * Return if the current residual index is meaningful or not
                 */
                return _current_residual_index_set;
            }

            std::string getResidualParameterizationInfo( ){
                /*!
                 * Get the parameterization information of the residual classes
                 */

                std::string parameterization_info = "########################################\n# RESIDUAL PARAMETERIZATION INFORMATION#\n########################################\n\n";

                for ( auto v = std::begin( *getResidualClasses( ) ); v != std::end( *getResidualClasses( ) ); ++v ){

                    parameterization_info += "RESIDUAL CLASS:";
                    parameterization_info += " " + std::to_string( ( unsigned int )( v - std::begin( *getResidualClasses( ) ) ) ) + "\n\n";

                    ( *v )->addParameterizationInfo( parameterization_info );

                    parameterization_info += "\n";

                }

                return parameterization_info;

            }

            void updateAdditionalStateVariables( ){
                /*!
                 * Update the additional state variable vector
                 */

                for ( auto v = std::begin( *getResidualClasses( ) ); v != std::end( *getResidualClasses( ) ); ++v ){

                    ( *v )->updateAdditionalStateVariables( _additionalStateVariables.second );

                }

            }

        protected:

            // Setters that the user may need to access but not override

            const void setCurrentResidualIndexMeaningful( const bool &value ){
                /*!
                 * Set if the current residual index is meaningful
                 * 
                 * \param &value: Set if the current residual index is meaningful or not
                 */

                _current_residual_index_set = value;

            }

            const void setCurrentResidualIndex( const unsigned int value ){
                /*!
                 * Set if the current residual index is meaningful
                 * 
                 * \param value: Set the value of the current residual index
                 */

                _current_residual_index = value;

            }

            void setStress( const floatVector &stress );

            // Utility functions
            virtual void computeConfigurations( const floatVector *data_vector, const unsigned int start_index,
                                                const floatVector &total_transformation,
                                                floatVector &configurations, floatVector &inverseConfigurations,
                                                const bool add_eye=false );

            virtual void extractStress( );

            virtual void updateConfigurationsFromUnknownVector( );

            virtual void decomposeUnknownVector( );

            virtual void decomposeStateVariableVector( );

            virtual void formNonLinearResidual( );

            virtual void formNonLinearDerivatives( );

            virtual void formPreconditioner( );

            virtual void solveNonLinearProblem( );

            virtual void formMaxRowPreconditioner( );

            virtual void initializeUnknownVector( );

            virtual void setTolerance( );

            virtual void initializePreconditioner( );

            virtual void evaluateInternal( );

            //! Update the line-search lambda parameter
            virtual void updateLambda( ){ _lambda *= 0.5; }

            virtual void updateUnknownVector( const floatVector &newUnknownVector );

            virtual void calculateFirstConfigurationJacobians( const floatVector &configurations, fourthOrderTensor &dC1dC, floatVector &dC1dCn );

            virtual void performArmijoTypeLineSearch( const floatVector &X0, const floatVector &deltaX );

            virtual void performGradientStep( const floatVector &X0 );

            //! Update the scaled quantities
            virtual void setScaledQuantities( );

            const floatType *get_baseResidualNorm( );

            const floatVector *get_basedResidualNormdX( );

            dataStorage< floatType > _baseResidualNorm; //!< The base value of the norm of the residual

            dataStorage< floatVector > _basedResidualNormdX; //!< The base value of the derivative of the norm of the residual w.r.t. the unknown vector

            void set_baseResidualNorm( const floatType &value ){ /*! Set the base value of the residual norm \param &value: The new value */  setNLStepData( value, _baseResidualNorm ); }

            void set_basedResidualNormdX( const floatVector &value ){ /*! Set the base derivative of the residual norm w.r.t. the unknown vector \param &value: The new value */  setNLStepData( value, _basedResidualNormdX ); }

            virtual void setBaseQuantities( );

            void resetNumNewton( ){ /*! Reset the number of newton steps */ _NUM_NEWTON = 0; }

            void resetNumLS( ){ /*! Reset the number of line search steps */ _NUM_LS = 0; }

            void resetNumGrad( ){ /*! Reset the number of gradient descent steps */ _NUM_GRAD = 0; }

            void incrementNumNewton( ){ /*! Reset the number of newton steps */ _NUM_NEWTON++; }

            void incrementNumLS( ){ /*! Reset the number of line search steps */ _NUM_LS++; }

            void incrementNumGrad( ){ /*! Reset the number of gradient descent steps */ _NUM_GRAD++; }

            void resetIterations( ){ /*! Reset the number of iterations */ _iteration = 0; }

            virtual void callResidualSuccessfulNLStep( );

            virtual void callResidualPreNLSolve( );

            virtual void callResidualPostNLSolve( );

            template<class T>
            void setIterationData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding iteration data. These values are cleared
                 * every time the unknown vector is updated.
                 *
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

                addIterationData( &storage );

            }

            template<class T>
            void setNLStepData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding nonlinear step data. These values are cleared
                 * every time the nonlinear step is advanced.
                 *
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

                addNLStepData( &storage );

            }

            template<class T>
            void setPreviousData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding previous data
                 * 
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

            }

            template<class T>
            void setConstantData( const T &data, dataStorage<T> &storage ){
                /*!
                 * Template function for adding constant data
                 * 
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

            }

            void unexpectedError( ){
                /*!
                 * Function to throw for an unexpected error. A user should never get here!
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "You shouldn't have gotten here. If you aren't developing the code then contact a developer with the stack trace." ) )

            }

            //!A pass through function that does nothing
            void passThrough( ){ }

            void setX( const floatVector &X ){
                /*!
                 * Set the value of the unknown vector
                 * 
                 * \param &X: The unknown vector
                 */

                _X.second = X;

                _X.first = true;

            }

            virtual void setResidualNorm( );

            virtual void setdResidualNormdX( );

            std::string build_upper_index_out_of_range_error_string( const unsigned int upperIndex, const unsigned int num_configurations );
            std::string build_lower_index_out_of_range_error_string( const unsigned int lowerIndex, const unsigned int upperIndex );

            floatVector _initialX; //!< The initial value of the unknown vector

            //! A data storage class that resets at every iteration
            template< typename T >
            class setDataStorageIteration : public setDataStorageBase< T > {

              public:

                  setDataStorageIteration( dataStorage< T > *ds, hydraBase * rp ) : setDataStorageBase< T >( ds ), _rp( rp ){
                      /*!
                       * Create a data storage object that will be reset at each new iteration
                       * 
                       * \param *ds: The data storage object
                       * \param *rp: The base hydra class that contains the data storage object
                       */ }

                  //! Destructor object that adds the data storage object to the iteration data list
                  ~setDataStorageIteration( ){

                      _rp->addIterationData( this->_ds );

                  }

              protected:

                  hydraBase *_rp; //!< The containing hydraBase class

            };

            /*!
             * Class which defines setting values defined at the previous timestep
             */
            template< typename T >
            class setDataStoragePrevious : public setDataStorageBase< T > {

                public:

                    //! Create a data storage object that will be reset whenever the previous value gets reset
                    setDataStoragePrevious( dataStorage< T > *ds ) : setDataStorageBase< T >( ds ){
                        /*!
                         * Constructor for data storage objects for temporally previous objects
                         * 
                         * \param *ds: The data storage object to modify
                         */
                    }

            };

            /*!
             * Class which defines setting constant values regardless of the timestep
             */
            template< typename T >
            class setDataStorageConstant : public setDataStorageBase< T > {

                public:

                    setDataStorageConstant( dataStorage< T > *ds ) : setDataStorageBase< T >( ds ){
                        /*!
                         * Constructor for constant data storage objects
                         * 
                         * \param *ds: The data storage object
                         */
                    }

            };

            virtual tardigradeHydra::hydraBase::setDataStorageConstant<floatVector> get_setDataStorage_tolerance( );

            virtual tardigradeHydra::hydraBase::setDataStorageIteration<secondOrderTensor> get_setDataStorage_stress( );

            virtual void assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            virtual void updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            virtual void assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints );

            virtual void solveConstrainedQP( floatVector &dx, const unsigned int kmax=100 );

            virtual void setConstraints( );

            virtual void setConstraintJacobians( );

            virtual void initializeActiveConstraints( std::vector< bool > &active_constraints );

            void setUseSQPSolver( const unsigned int &value ){ /*! Set whether to use the SQP solver \param &value: The updated value */ _useSQPSolver = value; }

            void setInitializeUnknownVector( const bool &value ){
                /*!
                 * Set the initialize unknown vector flag
                 * 
                 * \param &value: The value of the flag
                 */

                _initializeUnknownVector = value;

            }

            virtual void performRelaxedSolve( );

            void setAllowModifyGlobalResidual( const bool value){
                /*!
                 * Set a flag for if the global residual can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_residual = value;
            }

            void setAllowModifyGlobalJacobian( const bool value){
                /*!
                 * Set a flag for if the global jacobian can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_jacobian = value;
            }

            void setAllowModifyGlobaldRdT( const bool value){
                /*!
                 * Set a flag for if the global dRdT can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_dRdT = value;
            }

            void setAllowModifyGlobaldRdF( const bool value){
                /*!
                 * Set a flag for if the global dRdF can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_dRdF = value;

            }

            void setAllowModifyGlobaldRdAdditionalDOF( const bool value){
                /*!
                 * Set a flag for if the global dRdAdditionalDOF can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_dRdAdditionalDOF = value;

            }

        private:

            // Friend classes
            friend class unit_test::hydraBaseTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

            unsigned int _dimension; //!< The spatial dimension of the problem

            unsigned int _configuration_unknown_count; //!< The number of unknowns required for a configuration. Used to ensure that the unknown and state variable vectors are the right size. Must be set by all inheriting classes. For 3D classical continuum this will be 9, for higher order theories this will change.

            unsigned int _stress_size; //!< The number of terms in the stress measures. For 3D classical continuum this will be 9, for higher order theories this will change.

            floatType _time; //!< The current time

            floatType _deltaTime; //!< The change in time

            floatType _temperature; //!< The current temperature

            floatType _previousTemperature; //!< The previous temperature

            secondOrderTensor _deformationGradient; //!< The current deformation gradient

            secondOrderTensor _previousDeformationGradient; //!< The previous deformation gradient

            floatVector _additionalDOF; //!< The current additional degrees of freedom

            floatVector _previousAdditionalDOF; //!< The previous additional degrees of freedom

            floatVector _previousStateVariables; //!< The previous state variables

            floatVector _parameters; //!< The model parameters

            floatType _scaled_time; //!< The current time scaled by the scaling factor

            floatType _scaled_deltaTime; //!< The change in time scaled by the scaling factor

            floatType _scaled_temperature; //!< The current temperature scaled by the scaling factor

            secondOrderTensor _scaled_deformationGradient; //!< The current deformation gradient scaled by the scaling factor

            floatVector _scaled_additionalDOF; //!< The current additional degrees of freedom scaled by the scaling factor

            unsigned int _numConfigurations; //!< The number of configurations

            unsigned int _numNonLinearSolveStateVariables; //!< The number of state variables which will be solved in the Newton-Raphson loop

            floatType _tolr; //!< The relative tolerance

            floatType _tola; //!< The absolute tolerance

            unsigned int _maxIterations; //!< The maximum number of allowable iterations

            unsigned int _maxLSIterations; //!< The maximum number of line-search iterations

            unsigned int _maxGradientIterations = 10; //!< The maximum number of gradient iterations

            unsigned int _NUM_NEWTON = 0; //!< The number of Newton steps performed

            unsigned int _NUM_LS = 0; //!< The number of line search steps performed

            unsigned int _NUM_GRAD = 0; //!< The number of gradient descent steps performed

            bool _use_LM_step = false; //!< Flag for whether to attempt a Levenberg-Marquardt step

            bool _use_relaxed_solve = true; //!< Flag for whether to attempt a relaxed solve in case of failure

            bool _use_gradient_descent = false; //!< Flag for whether to attempt a gradient descent step

            bool _rank_deficient_error = true; //!< Flag for whether a rank-deficient LHS will throw a convergence error

            bool _initializeUnknownVector = true; //!< Flag for whether to initialize the unknown vector in the non-linear solve

            unsigned int _maxRelaxedIterations = 5; //!< The number of allowed relaxed iterations

            floatType _mu_k = -1; //!< The Levenberg-Marquardt scaling parameter

            floatType _lsAlpha; //!< The line-search alpha value i.e., the term by which it is judged that the line-search is converging

            std::vector< dataBase* > _iterationData; //!< A vector of pointers to data which should be cleared at each iteration

            std::vector< dataBase* > _nlStepData; //!< A vector of pointers to data which should be cleared after each nonlinear step

            dataStorage< std::vector< residualBase* > > _residualClasses; //!< A vector of classes which compute the terms in the residual equation

            dataStorage< floatVector > _residual; //!< The residual vector for the global solve

            dataStorage< floatVector > _jacobian; //!< The jacobian matrix in row-major form for the global solve

            dataStorage< floatVector > _nonlinearRHS; //!< The right hand side vector for the Newton solve

            dataStorage< floatVector > _flatNonlinearLHS; //!< The left hand side vector for the Newton solve

            dataStorage< floatVector > _preconditioner; //!< The pre-conditioner matrix in row-major form for the global solve

            bool _use_preconditioner; //!< Flag for whether to pre-condition the Jacobian or not

            unsigned int _preconditioner_type; //<! The type of preconditioner to use

            bool _preconditioner_is_diagonal; //!< Flag for if the pre-conditioner only stores the diagonal elements

            floatType _gradientSigma = 1e-4; //!< The sigma parameter for the gradient descent step

            floatType _gradientBeta  = 0.9; //!< The beta parameter for the gradient descent step

            floatType _gradientRho   = 1e-8; //!< The rho parameter for the gradient descent step

            floatType _gradientP     = 2.1; //!< The p parameter for the gradient descent step

            floatType _lm_mu         = 1e-8; //!< The mu parameter for Levenberg-Marquardt iterations

            dataStorage< floatVector > _dRdF; //!< The gradient of the residual w.r.t. the deformation gradient in row-major form for the global solve

            dataStorage< floatVector > _dRdT; //!< The gradient of the residual w.r.t. the temperature for the global solve

            dataStorage< floatVector > _dRdAdditionalDOF; //!< The derivatives of the residual w.r.t. the additional degrees of freedom

            dataStorage< floatVector > _additionalDerivatives; //!< Additional derivatives of the residual

            dataStorage< floatVector > _X; //!< The unknown vector { stress, F1, ..., Fn, xi1, ..., xim }

            dataStorage< floatVector > _tolerance; //!< The tolerance vector for the non-linear solve

            dataStorage< floatType > _lsResidualNorm; //!< The reference residual norm for the line-search convergence criteria

            dataStorage< floatVector > _stress; //!< The stress in the current configuration as determined from the current state

            dataStorage< floatVector > _previousStress; //!< The previous value of the stress in the current configuration as determined from the previous state

            dataStorage< floatVector > _flatdXdF; //!< The total derivative of the unknown vector w.r.t. the deformation in row-major form

            dataStorage< floatVector > _flatdXdT; //!< The total derivative of the unknown vector w.r.t. the temperature

            dataStorage< floatVector > _flatdXdAdditionalDOF; //!< The total derivative of the unknown vector w.r.t. the additional DOF

            unsigned int _iteration = 0; //!< The current iteration of the non-linear problem

            unsigned int _LSIteration = 0; //!< The current line search iteration of the non-linear problem

            unsigned int _gradientIteration = 0; //!< The current gradient iteration of the non-linear problem

            floatType _lambda = 1;

            bool _useSQPSolver = false;

            void setFirstConfigurationJacobians( );

            void setPreviousFirstConfigurationJacobians( );

            void performPreconditionedSolve( floatVector &deltaX_tr );

            void solveNewtonUpdate( floatVector &deltaX_tr );

            void setTolerance( const floatVector &tolerance );

            void incrementIteration( ){ _iteration++; resetLSIteration( ); }

            void incrementLSIteration( ){ _LSIteration++; }

            void incrementGradientIteration( ){ _gradientIteration++; }

            void resetLSIteration( ){ _LSIteration = 0; _lambda = 1.0;
                                      _lsResidualNorm.second = tardigradeVectorTools::l2norm( *getResidual( ) );
                                      _lsResidualNorm.first = true; }

            void resetGradientIteration( ){ _gradientIteration = 0; }

            const floatType* getLambda( ){ return &_lambda; }

            bool checkIteration( ){ return _iteration < _maxIterations; }

            bool checkLSIteration( ){ return _LSIteration < _maxLSIterations; }

            bool checkGradientIteration( ){ return _gradientIteration < _maxGradientIterations; }

            void resetIterationData( );

            void resetNLStepData( );

            unsigned int _failure_verbosity_level = 0; //!< The verbosity level for failure.

            std::stringstream _failure_output; //!< Additional failure output information

            floatType _scale_factor = 1.0; //!< A scale factor applied to the incoming loading (deformation, temperature, etc.)

            floatType _cutback_factor = 0.5; //!< The factor by which the pseudo-time will be scaled if a solve fails

            floatType _growth_factor  = 1.2; //!< The factor by which the pseudo-time will be scaled if we can grow the pseudo-timestep

            unsigned int _num_good_control = 2; //!< The number of good iterations we need to have before we try and increase the timestep

            floatType _minDS = 1e-2; //!< The minimum allowable pseudo-timestep

            bool _allow_modify_global_residual = false; //!< Flag for if the global residual can be modified

            bool _allow_modify_global_jacobian = false; //!< Flag for if the global jacobian can be modified

            bool _allow_modify_global_dRdT = false; //!< Flag for if the global dRdT can be modified

            bool _allow_modify_global_dRdF = false; //!< Flag for if the global dRdF can be modified

            bool _allow_modify_global_dRdAdditionalDOF = false; //!< Flag for if the global dRdAdditionalDOF can be modified

            bool _current_residual_index_set = false; //!< Flag for whether the current residual index has been set

            int _current_residual_index = 0; //!< The current residual index

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, configurations,                       floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousConfigurations,               floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, inverseConfigurations,                floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousInverseConfigurations,        floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, nonLinearSolveStateVariables,         floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousNonLinearSolveStateVariables, floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, additionalStateVariables,             floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousAdditionalStateVariables,     floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, residualNorm,       floatType,          setResidualNorm )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dResidualNormdX,    floatVector,        setdResidualNormdX )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraints,         getConstraints,         constraints,         floatVector,       setConstraints )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraintJacobians, getConstraintJacobians, constraintJacobians, floatVector,       setConstraintJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, set_dF1dF,              get_dF1dF,              dF1dF,               secondOrderTensor, setFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, set_dF1dFn,             get_dF1dFn,             dF1dFn,              floatVector,       setFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(  private, set_previousdF1dF,      get_previousdF1dF,      previousdF1dF,       secondOrderTensor, setPreviousFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(  private, set_previousdF1dFn,     get_previousdF1dFn,     previousdF1dFn,      floatVector,       setPreviousFirstConfigurationJacobians )

    };

    //! A data storage class that updates at every iteration
    template< typename T >
    class setDataStorageIteration : public setDataStorageBase< T > {

      public:

          setDataStorageIteration( dataStorage< T > *ds, residualBase * rp ) : setDataStorageBase< T >( ds ), _rp( rp ){
                      /*!
                       * Create a data storage object that will be reset at each new iteration
                       * 
                       * \param *ds: The data storage object
                       * \param *rp: The residual class that contains the data storage object
                       */
          }

          //! The destructor that says that the data storage object has been set
          ~setDataStorageIteration( ){

              _rp->addIterationData( this->_ds );

          }

      protected:

          residualBase *_rp; //!< The containing residual base class

    };

    //! A data storage class that updates whenever the previous values change
    template< typename T >
    class setDataStoragePrevious : public setDataStorageBase< T > {

        public:

            setDataStoragePrevious( dataStorage< T > *ds ) : setDataStorageBase< T >( ds ){
                /*!
                 * Constructor for data storage objects for temporally previous objects
                 * 
                 * \param *ds: The data storage object to modify
                 */
            }

    };

    //! A data storage class that is constant
    template< typename T >
    class setDataStorageConstant : public setDataStorageBase< T > {

        public:

            setDataStorageConstant( dataStorage< T > *ds ) : setDataStorageBase< T >( ds ){
                /*!
                 * Constructor for constant data storage objects
                 * 
                 * \param *ds: The data storage object
                 */
            }

    };

    /// Say hello
    /// @param message The message to print
    void sayHello(std::string message);

    void abaqusInterface( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                          double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                          const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                          const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                          const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                          const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                          const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                          const int *JSTEP,     const int &KINC );

    void dummyMaterialModel( floatVector &stress,             floatVector &statev,        floatMatrix &ddsdde,       floatType &SSE,            floatType &SPD,
                             floatType &SCD,                  floatType &RPL,             floatVector &ddsddt,       floatVector &drplde,       floatType &DRPLDT,
                             const floatVector &strain,       const floatVector &dstrain, const floatVector &time,   const floatType &DTIME,    const floatType &TEMP,
                             const floatType &DTEMP,          const floatVector &predef,  const floatVector &dpred,  const std::string &cmname, const int &NDI,
                             const int &NSHR,                 const int &NTENS,           const int &NSTATV,         const floatVector &props,  const int &NPROPS,
                             const floatVector &coords,       const floatMatrix &drot,    floatType &PNEWDT,         const floatType &CELENT,   const floatMatrix &dfgrd0,
                             const floatMatrix &dfgrd1,       const int &NOEL,            const int &NPT,            const int &LAYER,          const int &KSPT,
                             const std::vector< int > &jstep, const int &KINC );

}

#endif
