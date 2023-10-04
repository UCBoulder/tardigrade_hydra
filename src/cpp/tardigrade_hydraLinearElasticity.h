/**
  ******************************************************************************
  * \file tardigrade-hydraLinearElasticity.h
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
        
                const floatVector* getEe( );

                const floatMatrix* getdEedFe( );
       
                const floatVector* getPK2Stress( );
    
                const floatMatrix* getdPK2StressdEe( );

                const floatMatrix* getdPK2StressdFe( );
    
                const floatMatrix* getdCauchyStressdPK2Stress( );

                const floatMatrix* getdCauchyStressdF( );
    
                const floatMatrix* getdCauchyStressdFn( );
    
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
     
                void setEe( const floatVector &Ee );
        
                void setdEedFe( const floatMatrix &dEedFe );
        
                void setPK2Stress( const floatVector &PK2Stress );
    
                void setdPK2StressdEe( const floatMatrix &dPK2StressdEe );

                void setdPK2StressdFe( const floatMatrix &dPK2StressdFe );
    
                void setdCauchyStressdPK2Stress( const floatMatrix &dCauchyStressdPK2Stress );

                void setdCauchyStressdF( const floatMatrix &dCauchyStressdF );
    
                void setdCauchyStressdFn( const floatMatrix &dCauchyStressdFn );
    
            private:
    
                // Friend classes
                friend class tardigradeHydra::linearElasticity::unit_test::residualTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!
        
                using tardigradeHydra::residualBase::residualBase;
        
                using tardigradeHydra::residualBase::setResidual;
        
                using tardigradeHydra::residualBase::setJacobian;
        
                using tardigradeHydra::residualBase::setdRdF;
        
                using tardigradeHydra::residualBase::setdRdT;
        
                using tardigradeHydra::residualBase::setAdditionalDerivatives;
        
                using tardigradeHydra::residualBase::setCauchyStress;
        
                floatType _lambda;
        
                floatType _mu;
        
                tardigradeHydra::dataStorage< floatVector > _Ee;

                tardigradeHydra::dataStorage< floatMatrix > _dEedFe;

                tardigradeHydra::dataStorage< floatVector > _PK2Stress;

                tardigradeHydra::dataStorage< floatMatrix > _dPK2StressdEe;

                tardigradeHydra::dataStorage< floatMatrix > _dPK2StressdFe;

                tardigradeHydra::dataStorage< floatMatrix > _dCauchyStressdPK2Stress;

                tardigradeHydra::dataStorage< floatMatrix > _dCauchyStressdF;

                tardigradeHydra::dataStorage< floatMatrix > _dCauchyStressdFn;

                virtual void setEe( );

                virtual void setdEedFe( );
 
                virtual void setPK2Stress( );            

                virtual void setdPK2StressdEe( );

                virtual void setdPK2StressdFe( );

                virtual void setCauchyStress( ) override;

                virtual void setdCauchyStressdPK2Stress( );

                virtual void setdCauchyStressdF( );

                virtual void setdCauchyStressdFn( );

                virtual void setResidual( ) override;

                virtual void setJacobian( ) override;

                virtual void setdRdT( ) override;

                virtual void setdRdF( ) override;

                virtual void decomposeParameterVector( const floatVector &parameters );
    
        };

    }

}

#endif
