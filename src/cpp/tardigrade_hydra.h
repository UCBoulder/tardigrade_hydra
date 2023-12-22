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
#include<tardigrade_abaqus_tools.h>

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
        setfun( varname, _##varname );                                        \
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
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setfun,uncall) \
    private: dataStorage<vartype> _##varname;                       \
    public: TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER(setname,varname,vartype,setfun) \
    public: TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(getname,varname,vartype,uncall) \
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
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,set_##varname,get_##varname,varname,vartype,setfun,uncall)

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
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setIterationData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setIterationData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_STORAGE(context,varname,vartype,setIterationData,uncall)

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
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setPreviousData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setPreviousData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_STORAGE(context,varname,vartype,setPreviousData,uncall)

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
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setConstantData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setContantData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_HYDRA_DECLARE_STORAGE(context,varname,vartype,setConstantData,uncall)

namespace tardigradeHydra{

    // forward class definitions
    class hydraBase;

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

    typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
    typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    typedef void ( hydraBase::*hydraBaseFxn )( ); //!< Typedef for passing pointers to hydraBase functions

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
    inline void dataStorage< unsigned int >::clear( ){
                /*!
                 * The function to erase the current values stored by setting first to false and second to zero
                 */

        first = false;

        second = 0;

    }

    template <>
    inline void dataStorage< floatType >::clear( ){
                /*!
                 * The function to erase the current values stored by setting first to false and second to zero
                 */

        first = false;

        second = 0;

    }

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
            residualBase( residualBase &r ) : hydra( r.hydra ), _numEquations( *r.getNumEquations( ) ){ }

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

            virtual void setCurrentAdditionalStateVariables( ){
                /*!
                 * Set the current additional state variables
                 * 
                 * Doesn't need to be defined
                 */

                setCurrentAdditionalStateVariables( floatVector( 0, 0 ) );

                return;

            }

            // Getter functions

            //! Get the number of equations the residual defined
            const unsigned int* getNumEquations( ){ return &_numEquations; }

            void addIterationData( dataBase *data );

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

        protected:

            void unexpectedError( ){
                /*!
                 * Function to throw for an unexpected error. A user should never get here!
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "You shouldn't have gotten here. If you aren't developing the code then contact a developer with the stack trace." ) )

            }

        private:

            unsigned int _numEquations; //!< The number of residual equations

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setResidual,                         getResidual,                        residual,                        floatVector, setResidual )
                                                                                                                                              
            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setJacobian,                         getJacobian,                        jacobian,                        floatMatrix, setJacobian )
                                                                                                                                              
            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdF,                             getdRdF,                            dRdF,                            floatMatrix, setdRdF )
                                                                                                                                              
            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdT,                             getdRdT,                            dRdT,                            floatVector, setdRdT )
                                                                                                                                              
            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setAdditionalDerivatives,            getAdditionalDerivatives,           additionalDerivatives,           floatMatrix, setAdditionalDerivatives )
                                                                                                                                              
            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setStress,                           getStress,                          stress,                          floatVector, setStress )
                                                                                                                                              
            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setPreviousStress,                   getPreviousStress,                  previousStress,                  floatVector, setPreviousStress )
                                                                                                                                              
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
                       const floatVector &deformationGradient, const floatVector &previousDeformationGradient,
                       const floatVector &previousStateVariables, const floatVector &parameters,
                       const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                       const unsigned int dimension=3, const unsigned int configuration_unknown_count=9,
                       const floatType tolr=1e-9, const floatType tola=1e-9,
                       const unsigned int maxIterations=20, const unsigned int maxLSIterations=5, const floatType lsAlpha=1e-4 );

            // User defined functions

            // Setter functions

            // Getter functions
            //! Get a reference to the number of unknowns in each configuration
            const unsigned int* getConfigurationUnknownCount( ){ return &_configuration_unknown_count; }

            //! Get a reference to the current time
            const floatType* getTime( ){ return &_time; }

            //! Get a reference to the change in time
            const floatType* getDeltaTime( ){ return &_deltaTime; }

            //! Get a reference to the current temperature
            const floatType* getTemperature( ){ return &_temperature; };

            //! Get a reference to the previous temperature
            const floatType* getPreviousTemperature( ){ return &_previousTemperature; };

            //! Get a reference to the deformation gradient
            const floatVector* getDeformationGradient( ){ return &_deformationGradient; }

            //! Get a reference to the previous deformation gradient
            const floatVector* getPreviousDeformationGradient( ){ return &_previousDeformationGradient; }

            //! Get a reference to the previous values of the state variables
            const floatVector* getPreviousStateVariables( ){ return &_previousStateVariables; }

            //! Get a reference to the model parameters
            const floatVector* getParameters( ){ return &_parameters; }

            //! Get a reference to the number of configurations
            const unsigned int* getNumConfigurations( ){ return &_numConfigurations; }

            //! Get a reference to the number of state variables involved in the non-linear solve
            const unsigned int* getNumNonLinearSolveStateVariables( ){ return &_numNonLinearSolveStateVariables; }

            //! Get a reference to the dimension
            const unsigned int* getDimension( ){ return &_dimension; }

            //! Get a reference to the relative tolerance
            const floatType* getRelativeTolerance( ){ return &_tolr; }

            //! Get a reference to the absolute tolerance
            const floatType* getAbsoluteTolerance( ){ return &_tola; }

            //! Get a reference to the line-search alpha
            const floatType* getLSAlpha( ){ return &_lsAlpha; }

            floatVector getSubConfiguration( const floatMatrix &configurations, const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatMatrix getSubConfigurationJacobian( const floatMatrix &configurations, const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPrecedingConfiguration( const unsigned int &index );

            floatVector getFollowingConfiguration( const unsigned int &index );

            floatVector getConfiguration( const unsigned int &index );

            floatVector getPreviousSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPreviousPrecedingConfiguration( const unsigned int &index );

            floatVector getPreviousFollowingConfiguration( const unsigned int &index );

            floatVector getPreviousConfiguration( const unsigned int &index );

            floatMatrix getSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatMatrix getPrecedingConfigurationJacobian( const unsigned int &index );

            floatMatrix getFollowingConfigurationJacobian( const unsigned int &index );

            floatMatrix getPreviousSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatMatrix getPreviousPrecedingConfigurationJacobian( const unsigned int &index );

            floatMatrix getPreviousFollowingConfigurationJacobian( const unsigned int &index );

            const floatType* getLSResidualNorm( );

            virtual void setResidualClasses( );

            void setResidualClasses( std::vector< residualBase* > &residualClasses );

            std::vector< residualBase* >* getResidualClasses( );

            const floatVector* getResidual( );

            const floatVector* getFlatJacobian( );

            floatMatrix getJacobian( );

            const floatVector* getFlatdRdF( );

            floatMatrix getdRdF( );

            const floatVector* getdRdT( );

            const floatVector* getFlatAdditionalDerivatives( );

            floatMatrix getAdditionalDerivatives( );

            const floatVector* getUnknownVector( );

            const floatVector* getTolerance( );

            virtual bool checkConvergence( );

            virtual bool checkLSConvergence( );

            const floatVector* getStress( );

            const floatVector* getPreviousStress( );

            virtual void evaluate( );

            //! Add data to the vector of values which will be cleared after each iteration
            void addIterationData( dataBase *data ){ _iterationData.push_back( data ); }

        protected:

            // Setters that the user may need to access but not override

            void setStress( const floatVector &stress );

            // Utility functions
            virtual void computeConfigurations( const floatVector *data_vector, const unsigned int start_index,
                                                const floatVector &total_transformation,
                                                floatMatrix &configurations, floatMatrix &inverseConfigurations,
                                                const bool add_eye=false );

            virtual void extractStress( );

            virtual void decomposeUnknownVector( );

            virtual void decomposeStateVariableVector( );

            virtual void formNonLinearProblem( );

            virtual void initializeUnknownVector( );

            virtual void setTolerance( );

            //! Update the line-search lambda parameter
            virtual void updateLambda( ){ _lambda *= 0.5; }

            virtual void updateUnknownVector( const floatVector &newUnknownVector );

            virtual void calculateFirstConfigurationJacobians( const floatMatrix &configurations, floatMatrix &dC1dC, floatMatrix &dC1dCn );

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

        private:

            // Friend classes
            friend class unit_test::hydraBaseTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

            unsigned int _dimension; //!< The spatial dimension of the problem

            const unsigned int _configuration_unknown_count; //!< The number of unknowns required for a configuration. Used to ensure that the unknown and state variable vectors are the right size. Must be set by all inheriting classes. For 3D classical continuum this will be 9, for higher order theories this will change.

            floatType _time; //!< The current time

            floatType _deltaTime; //!< The change in time

            floatType _temperature; //!< The current temperature

            floatType _previousTemperature; //!< The previous temperature

            floatVector _deformationGradient; //!< The current deformation gradient

            floatVector _previousDeformationGradient; //!< The previous deformation gradient

            floatVector _previousStateVariables; //!< The previous state variables

            floatVector _parameters; //!< The model parameters

            unsigned int _numConfigurations; //!< The number of configurations

            unsigned int _numNonLinearSolveStateVariables; //!< The number of state variables which will be solved in the Newton-Raphson loop

            floatType _tolr; //!< The relative tolerance

            floatType _tola; //!< The absolute tolerance

            unsigned int _maxIterations; //!< The maximum number of allowable iterations

            unsigned int _maxLSIterations; //!< The maximum number of line-search iterations

            floatType _lsAlpha; //!< The line-search alpha value i.e., the term by which it is judged that the line-search is converging

            std::vector< dataBase* > _iterationData; //!< A vector of pointers to data which should be cleared at each iteration

            dataStorage< std::vector< residualBase* > > _residualClasses; //!< A vector of classes which compute the terms in the residual equation

            dataStorage< floatVector > _residual; //!< The residual vector for the global solve

            dataStorage< floatVector > _jacobian; //!< The jacobian matrix in row-major form for the global solve

            dataStorage< floatVector > _dRdF; //!< The gradient of the residual w.r.t. the deformation gradient in row-major form for the global solve

            dataStorage< floatVector > _dRdT; //!< The gradient of the residual w.r.t. the temperature for the global solve

            dataStorage< floatVector > _additionalDerivatives; //!< Additional derivatives of the residual

            dataStorage< floatVector > _X; //!< The unknown vector { stress, F1, ..., Fn, xi1, ..., xim }

            dataStorage< floatVector > _tolerance; //!< The tolerance vector for the non-linear solve

            dataStorage< floatType > _lsResidualNorm; //!< The reference residual norm for the line-search convergence criteria

            dataStorage< floatVector > _stress; //!< The stress in the current configuration as determined from the current state

            dataStorage< floatVector > _previousStress; //!< The previous value of the stress in the current configuration as determined from the previous state

            unsigned int _iteration = 0; //!< The current iteration of the non-linear problem

            unsigned int _LSIteration = 0; //!< The current line search iteration of the non-linear problem

            floatType _lambda = 1;

            void setFirstConfigurationJacobians( );

            void setPreviousFirstConfigurationJacobians( );

            void solveNonLinearProblem( );

            void setTolerance( const floatVector &tolerance );

            void incrementIteration( ){ _iteration++; }

            void incrementLSIteration( ){ _LSIteration++; }

            void resetLSIteration( ){ _LSIteration = 0; _lambda = 1.0;
                                      _lsResidualNorm.second = tardigradeVectorTools::l2norm( *getResidual( ) );
                                      _lsResidualNorm.first = true; }

            const floatType* getLambda( ){ return &_lambda; }

            bool checkIteration( ){ return _iteration < _maxIterations; }

            bool checkLSIteration( ){ return _LSIteration < _maxLSIterations; }

            void resetIterationData( );

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, configurations,                       floatMatrix, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousConfigurations,               floatMatrix, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, inverseConfigurations,                floatMatrix, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousInverseConfigurations,        floatMatrix, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, nonLinearSolveStateVariables,         floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousNonLinearSolveStateVariables, floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, additionalStateVariables,             floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousAdditionalStateVariables,     floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, set_dF1dF,          get_dF1dF,          dF1dF,          floatMatrix, setFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, set_dF1dFn,         get_dF1dFn,         dF1dFn,         floatMatrix, setFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(  private, set_previousdF1dF,  get_previousdF1dF,  previousdF1dF,  floatMatrix, setPreviousFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(  private, set_previousdF1dFn, get_previousdF1dFn, previousdF1dFn, floatMatrix, setPreviousFirstConfigurationJacobians )

    };

    /// Say hello
    /// @param message The message to print
    errorOut sayHello(std::string message);

    void abaqusInterface( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                          double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                          const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                          const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                          const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                          const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                          const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                          const int *JSTEP,     const int &KINC );

    errorOut dummyMaterialModel( floatVector &stress,             floatVector &statev,        floatMatrix &ddsdde,       floatType &SSE,            floatType &SPD,
                             floatType &SCD,                  floatType &RPL,             floatVector &ddsddt,       floatVector &drplde,       floatType &DRPLDT,
                             const floatVector &strain,       const floatVector &dstrain, const floatVector &time,   const floatType &DTIME,    const floatType &TEMP,
                             const floatType &DTEMP,          const floatVector &predef,  const floatVector &dpred,  const std::string &cmname, const int &NDI,
                             const int &NSHR,                 const int &NTENS,           const int &NSTATV,         const floatVector &props,  const int &NPROPS,
                             const floatVector &coords,       const floatMatrix &drot,    floatType &PNEWDT,         const floatType &CELENT,   const floatMatrix &dfgrd0,
                             const floatMatrix &dfgrd1,       const int &NOEL,            const int &NPT,            const int &LAYER,          const int &KSPT,
                             const std::vector< int > &jstep, const int &KINC );

}

#endif
