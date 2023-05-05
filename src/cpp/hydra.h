/**
  ******************************************************************************
  * \file hydra.h
  ******************************************************************************
  * A C++ library for printing messages to stdout. Used as a stub repo example.
  ******************************************************************************
  */

#ifndef HYDRA_H
#define HYDRA_H

#include<sstream>
#include<functional>

#include<error_tools.h>
#define USE_EIGEN
#include<vector_tools.h>
#include<abaqus_tools.h>

namespace hydra{

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
    const std::string __BASENAME__ = file_name(__FILE__);
    const std::string __FILENAME__ = __BASENAME__.substr(0, __BASENAME__.find_last_of("."));

    typedef errorTools::Node errorNode; //!< Redefinition for the error node
    typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    typedef void ( hydraBase::*hydraBaseFxn )( ); //!< Typedef for passing pointers to hydraBase functions

    class dataBase{

        public:

            virtual void clear( ){
                /*!
                 * The function to erase the current values stored
                 */

                ERROR_TOOLS_CATCH( throw std::runtime_error( "clear not implemented!" ) );

            }

    };

    template < typename T >
    class dataStorage : public dataBase{

        public:

            bool first = false; //!The flag for whether the data has been stored

            T second; //!The stored data

            dataStorage( ){ };

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

    class hydraBase{
        /*!
         * hydraBase: A base class which can be used to construct finite deformation material models.
         * 
         * The hydra class seeks to provide utilities for the construction of finite deformation constitutive models
         * more rapidly than would be possible previously. The user can define as many different configurations as desired
         * and provide a calculation of the Cauchy stress.
         * 
         * A non-linear problem which is of the size ( dimension**2 * num_configurations + num_ISVs ) will be solved.
         */

        public:

            // Constructors
            //! Default constructor for hydraBase
            hydraBase( ){ }

            //! Main constructor for objects of type hydraBase. Sets all quantities required for most solves.
            hydraBase( const floatType &time, const floatType &deltaTime,
                       const floatType &temperature, const floatType &previousTemperature,
                       const floatVector &deformationGradient, const floatVector &previousDeformationGradient,
                       const floatVector &previousStateVariables, const floatVector &parameters,
                       const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                       const unsigned int dimension=3 ) : _time( time ), _deltaTime( deltaTime ),
                                                          _temperature( temperature ), _previousTemperature( previousTemperature ),
                                                          _deformationGradient( deformationGradient ),
                                                          _previousDeformationGradient( previousDeformationGradient ),
                                                          _previousStateVariables( previousStateVariables ),
                                                          _parameters( parameters ),
                                                          _numConfigurations( numConfigurations ),
                                                          _numNonLinearSolveStateVariables( numNonLinearSolveStateVariables ),
                                                          _dimension( dimension ){ }

            // User defined functions

            // Setter functions

            // Getter functions
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

            //! Get a reference to the configurations
            const floatMatrix* getConfigurations( ){ return &_configurations.second; }

            //! Get a reference to the previous configurations
            const floatMatrix* getPreviousConfigurations( ){ return &_previousConfigurations.second; }

            //! Get a reference to the inverse configurations
            const floatMatrix* getInverseConfigurations( ){ return &_inverseConfigurations.second; }

            //! Get a reference to the previous inverse configurations
            const floatMatrix* getPreviousInverseConfigurations( ){ return &_previousInverseConfigurations.second; }

            //! Get a reference to the state variables used in the unknown vector of the non-linear solve
            const floatVector* getNonLinearSolveStateVariables( ){ return &_nonLinearSolveStateVariables.second; }

            //! Get a reference to the previous values of the state variables used in the unknown vector of the non-linear solve
            const floatVector* getPreviousNonLinearSolveStateVariables( ){ return &_previousNonLinearSolveStateVariables.second; }

            //! Get a reference to the state variables not used in the unknown vector for the non-linear solve
            const floatVector* getAdditionalStateVariables( ){ return &_additionalStateVariables.second; }

            //! Get a reference to the previous values of the state variables not used in the unknown vector for the non-linear solve
            const floatVector* getPreviousAdditionalStateVariables( ){ return &_previousAdditionalStateVariables.second; }

        private:

            // Friend classes
            friend class unit_test::hydraBaseTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

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
            unsigned int _dimension; //!< The spatial dimension of the problem

            dataStorage< floatMatrix > _configurations; //!< The current values of the configurations

            dataStorage< floatMatrix > _previousConfigurations; //!< The previous values of the configurations

            dataStorage< floatMatrix > _inverseConfigurations; //!< The inverses of the configurations

            dataStorage< floatMatrix > _previousInverseConfigurations; //!< The inverses of the previous configurations

            dataStorage< floatVector > _nonLinearSolveStateVariables; //!< The current values of the state variables involved in the non-linear solve
            dataStorage< floatVector > _previousNonLinearSolveStateVariables; //!< The previous values of the state variables involved in the non-linear solve

            dataStorage< floatVector > _additionalStateVariables; //!< The current values of the additional state variables

            dataStorage< floatVector > _previousAdditionalStateVariables; //!< The previous values of the additional state variables

            virtual void decomposeStateVariableVector( );

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
