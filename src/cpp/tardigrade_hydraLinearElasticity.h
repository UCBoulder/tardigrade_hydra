/**
  ******************************************************************************
  * \file tardigrade_hydraLinearElasticity.h
  ******************************************************************************
  * An implementation of linear elasticity using the hydra framework. Used as an
  * example and as the basis for more complex models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_LINEAR_ELASTICITY_H
#define TARDIGRADE_HYDRA_LINEAR_ELASTICITY_H

#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    namespace linearElasticity{

        // forward class definitions
        namespace unit_test{
            class residualTester;
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
        const std::string __BASENAME__ = file_name(__FILE__);  //!< The base filename which will be parsed
    const std::string __FILENAME__ = __BASENAME__.substr(0, __BASENAME__.find_last_of(".")); //!< The parsed filename for error handling
    
        typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
        typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
        typedef double floatType; //!< Define the float values type.
        typedef std::vector< floatType > floatVector; //!< Define a vector of floats
        typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats 

        /*!
         * A residual class for a linear-elastic material model where the stress is computed
         * in the reference configuration and pushed forward to the current configuration.
         */
        class residual : public tardigradeHydra::residualBase {
        
            public:

                /*!
                 * The main initialization constructor for the linear elastic residual
                 * 
                 * \param *hydra: A pointer to the containing hydra class
                 * \param &numEquations: The number of equations the residual defines
                 * \param &parameters: The parameter vector
                 */
                residual( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations, const floatVector &parameters ) : tardigradeHydra::residualBase( hydra, numEquations ){
    
                    TARDIGRADE_ERROR_TOOLS_CATCH( decomposeParameterVector( parameters ) );
    
                }
                //! Get a pointer to the value of the lambda Lame parameter
                const floatType* getLambda( ){ return &_lambda; }

                //! Get a pointer to the value of the mu Lame parameter
                const floatType* getMu( ){ return &_mu; }
        
                /*!
                 * Set the value of the lambda Lame parameter
                 * 
                 * \param &lambda: The lambda Lame parameter
                 */ 
                void setLambda( const floatType &lambda ){ _lambda = lambda; }
   
                /*!
                 * Set the value of the mu Lame parameter
                 * 
                 * \param &mu: The mu Lame parameter
                 */ 
                void setMu( const floatType &mu ){ _mu = mu; }

            protected:

                virtual void setFe( const bool isPrevious );

                virtual void setFe( );

                virtual void setPreviousFe( );

                virtual void setFeDerivatives( const bool isPrevious );

                virtual void setdFedF( );

                virtual void setdFedFn( );

                virtual void setPreviousdFedF( );

                virtual void setPreviousdFedFn( );

                virtual void setEe( const bool isPrevious );

                virtual void setEe( );

                virtual void setdEedFe( );
 
                virtual void setPreviousEe( );

                virtual void setPreviousdEedFe( );

                virtual void setPK2Stress( const bool isPrevious );
 
                virtual void setPK2Stress( ); 

                virtual void setdPK2StressdEe( const bool isPrevious );

                virtual void setdPK2StressdEe( );

                virtual void setdPK2StressdFe( );

                virtual void setdPK2StressdPreviousFe( );

                virtual void setPreviousPK2Stress( ); 

                virtual void setPreviousdPK2StressdEe( );

                virtual void setPreviousdPK2StressdFe( );

                virtual void setStress( ) override;

                virtual void setPreviousStress( ) override;

                virtual void setCauchyStress( const bool isPrevious );

                virtual void setCauchyStress( );

                virtual void setPreviousCauchyStress( );

                virtual void setdCauchyStressdPK2Stress( );

                virtual void setdCauchyStressdF( );

                virtual void setdCauchyStressdFn( );

                virtual void setdCauchyStressdPreviousF( );

                virtual void setdCauchyStressdPreviousFn( );

                virtual void setPreviousdCauchyStressdPK2Stress( );

                virtual void setPreviousdCauchyStressdF( );

                virtual void setPreviousdCauchyStressdFn( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

                virtual void decomposeParameterVector( const floatVector &parameters );
     
            private:
    
                // Friend classes
                friend class tardigradeHydra::linearElasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!
        
                using tardigradeHydra::residualBase::residualBase;
        
                using tardigradeHydra::residualBase::setResidual;
        
                using tardigradeHydra::residualBase::setJacobian;
        
                using tardigradeHydra::residualBase::setdRdF;
        
                using tardigradeHydra::residualBase::setdRdT;
        
                using tardigradeHydra::residualBase::setAdditionalDerivatives;
        
                using tardigradeHydra::residualBase::setStress;

                using tardigradeHydra::residualBase::setPreviousStress;
        
                floatType _lambda;
        
                floatType _mu;
        
                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, Fe,                              secondOrderTensor, setFe                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFedF,                           fourthOrderTensor, setdFedF                           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dFedFn,                          floatVector,       setdFedFn                          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousFe,                      secondOrderTensor, setPreviousFe                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdFedF,                   fourthOrderTensor, setPreviousdFedF                   )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdFedFn,                  floatVector,       setPreviousdFedFn                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, Ee,                              secondOrderTensor, setEe                              )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dEedFe,                          fourthOrderTensor, setdEedFe                          )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousEe,                      secondOrderTensor, setPreviousEe                      )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdEedFe,                  fourthOrderTensor, setPreviousdEedFe                  )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, PK2Stress,                       secondOrderTensor, setPK2Stress                       )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2StressdEe,                   fourthOrderTensor, setdPK2StressdEe                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2StressdFe,                   fourthOrderTensor, setdPK2StressdFe                   )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dPK2StressdPreviousFe,           fourthOrderTensor, setdPK2StressdPreviousFe           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousPK2Stress,               secondOrderTensor, setPreviousPK2Stress               )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2StressdEe,           fourthOrderTensor, setPreviousdPK2StressdEe           )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdPK2StressdFe,           fourthOrderTensor, setPreviousdPK2StressdFe           )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, cauchyStress,                    secondOrderTensor, setCauchyStress                    )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdPK2Stress,         fourthOrderTensor, setdCauchyStressdPK2Stress         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdF,                 fourthOrderTensor, setdCauchyStressdF                 )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdFn,                floatVector,       setdCauchyStressdFn                )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdPreviousF,         fourthOrderTensor, setdCauchyStressdPreviousF         )

                TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dCauchyStressdPreviousFn,        floatVector,       setdCauchyStressdPreviousFn        )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousCauchyStress,            secondOrderTensor, setPreviousCauchyStress            )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdPK2Stress, fourthOrderTensor, setPreviousdCauchyStressdPK2Stress )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdF,         fourthOrderTensor, setPreviousdCauchyStressdF         )

                TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdCauchyStressdFn,        floatVector,       setPreviousdCauchyStressdFn        )

        };

    }

}

#endif
