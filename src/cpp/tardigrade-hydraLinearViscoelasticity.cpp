/**
  ******************************************************************************
  * \file tardigrade-hydraLinearViscoelasticity.cpp
  ******************************************************************************
  * An implementation of linear elasticity using the hydra framework. Used as an
  * example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade-hydraLinearViscoelasticity.h>
#include<constitutive_tools.h>
#include<stress_tools.h>

namespace tardigradeHydra{

    namespace linearViscoelasticity{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length 2 which defines lambda and mu.
             */
   
            const unsigned int *dim = hydra->getDimension( );
 
            if ( parameters.size( ) < 10 ){
    
                ERROR_TOOLS_CATCH( throw std::runtime_error( "Parameter vector is expected to have a length of at least 10 but has a length of " + std::to_string( parameters.size( ) ) ) );
    
            }

            setNumVolumetricViscousTerms( ( unsigned int )( parameters[ 0 ] + 0.5 ) );

            setNumIsochoricViscousTerms( ( unsigned int )( parameters[ 1 ] + 0.5 ) );

            setNumStateVariables( *getNumVolumetricViscousTerms( ) + ( *dim ) * ( *dim ) * ( *getNumIsochoricViscousTerms( ) ) );

            if ( *getNumStateVariables( ) != ( *getViscoelasticISVUpperIndex( ) - *getViscoelasticISVLowerIndex( ) ) ){

                std::string message = "The number of state variables required by the parameterization is not equal to the number of state variables indicated by the ISV bounds\n";
                message            += "   required # ISVs: " + std::to_string( *getNumStateVariables( ) ) + "\n";
                message            += "   ISV Lower Bound: " + std::to_string( *getViscoelasticISVLowerIndex( ) ) + "\n";
                message            += "   ISV UPper Bound: " + std::to_string( *getViscoelasticISVLowerIndex( ) ) + "\n";

                ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

            }

            setKinf( parameters[ 2 ] );

            setGinf( parameters[ 3 ] );

            setVolumetricTemperatureParameters( floatVector( parameters.begin( ) + 4,
                                                             parameters.begin( ) + 7 ) );

            setIsochoricTemperatureParameters( floatVector( parameters.begin( ) +  7,
                                                            parameters.begin( ) + 10 ) );

            unsigned int parameterCount = 10 + 2 * ( *getNumVolumetricViscousTerms( ) )
                                        + 2 * ( *getNumIsochoricViscousTerms( ) );

            if ( parameters.size( ) != parameterCount ){

                std::string message = "The number of parameters provided is not consistent with the parameter counts\n";
                message            += "  num parameters:      " + std::to_string( parameters.size( ) ) + "\n";
                message            += "  num viscous terms:   " + std::to_string( *getNumVolumetricViscousTerms( ) ) + "\n";
                message            += "  num isochoric terms: " + std::to_string( *getNumIsochoricViscousTerms( ) ) + "\n";
                message            += "The number of parameters is 4 + 2 * ( numVolumetricViscousTerms + numIsochoricViscousTerms )\n";
                message            += "  required parameter count: " + std::to_string( parameterCount ) + "\n";

                ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

            }

            unsigned int lb = 10;
            unsigned int ub = lb + *getNumVolumetricViscousTerms( );

            floatVector Ks( parameters.begin( ) + lb, parameters.begin( ) + ub );

            lb = ub;
            ub = lb + *getNumVolumetricViscousTerms( );

            floatVector Ktaus( parameters.begin( ) + lb, parameters.begin( ) + ub );

            lb = ub;
            ub = lb + *getNumIsochoricViscousTerms( );

            floatVector Gs( parameters.begin( ) + lb, parameters.begin( ) + ub );

            lb = ub;
            ub = lb + *getNumIsochoricViscousTerms( );

            floatVector Gtaus( parameters.begin( ) + lb, parameters.begin( ) + ub );

            setVolumetricModuli( Ks );

            setVolumetricTaus( Ktaus );

            setIsochoricModuli( Gs );

            setIsochoricTaus( Gtaus );

        }

        void residual::setNumVolumetricViscousTerms( const unsigned int &num ){
            /*!
             * Set the number of volumetric prony-series viscous terms
             * 
             * \param &num: The number of volumetric Prony-series terms
             */

            _numVolumetricViscousTerms = num;

        }

        void residual::setNumIsochoricViscousTerms( const unsigned int &num ){
            /*!
             * Set the number of isochoric prony-series viscous terms
             * 
             * \param &num: The number of isochorric Prony-series terms
             */

            _numIsochoricViscousTerms = num;

        }

        void residual::setKinf( const floatType &Kinf ){
            /*!
              * Set the infinite bulk modulus
              * 
              * \param &Kinf: The infinite bulk modulus
              */

            _Kinf = Kinf;

        }

        void residual::setGinf( const floatType &Ginf ){
            /*!
              * Set the infinite shear modulus
              * 
              * \param &Ginf: The infinite shear modulus
              */
        
            _Ginf = Ginf;

        }

        void residual::setVolumetricModuli( const floatVector &Ks ){
            /*!
             * Set the volumetric moduli
             * 
             * \param &Ks: The bulk moduli
             */

            _Ks = Ks;

        }

        void residual::setIsochoricModuli( const floatVector &Gs ){
            /*!
             * Set the isochoric moduli
             * 
             * \param &Gs: The isochoric moduli
             */

            _Gs = Gs;

        }

        void residual::setVolumetricTaus( const floatVector &taus ){
            /*!
             * Set the volumetric time constants
             * 
             * \param &taus: The bulk time constants
             */

            _volumetricTaus = taus;

        }

        void residual::setIsochoricTaus( const floatVector &taus ){
            /*!
             * Set the isochoric time constants
             * 
             * \param &Gs: The isochoric time constants
             */

            _isochoricTaus = taus;

        }

        void residual::decomposeElasticDeformation( ){
            /*!
             * Decompose the elastic deformation into volumetric and isochoric parts
             */

            floatType Je;

            floatVector Fehat;

            ERROR_TOOLS_CATCH( decomposeDeformation( ( *hydra->getConfigurations( ) )[ 0 ], Je, Fehat ) );

            setJe( Je );

            setFehat( Fehat );

        }

        void residual::decomposePreviousElasticDeformation( ){
            /*!
             * Decompose the previous elastic deformation into volumetric and isochoric parts
             */

            floatType previousJe;

            floatVector previousFehat;

            ERROR_TOOLS_CATCH( decomposeDeformation( ( *hydra->getPreviousConfigurations( ) )[ 0 ], previousJe, previousFehat ) );

            setPreviousJe( previousJe );

            setPreviousFehat( previousFehat );

        }

        void residual::decomposeDeformation( const floatVector &F, floatType &J, floatVector &Fhat ){
            /*!
             * Decompose a deformation into volumetric and isochoric parts where
             * 
             * \f$\hat{\bf{F}} = J^{-\frac{1}{3}} \bf{F} \f$
             * 
             * \param &F: The incoming deformation gradient
             * \param &J: The jacobian of deformation
             * \param &Fhat: The isochoric part of the deformation gradient
             */

            const unsigned int *dim = hydra->getDimension( );

            ERROR_TOOLS_CATCH( J = vectorTools::determinant( F, ( *dim ), ( *dim ) ) );

            ERROR_TOOLS_CATCH( Fhat = F / std::pow( J, 1./3 ) );

        }

        void residual::setJe( const floatType &Je ){
            /*!
             * Set the elastic Jacobian of deformation
             * 
             * \param &Je: The elastic Jacobian of deformation
             */

            _Je.second = Je;

            _Je.first = true;

            addIterationData( &_Je );

        }

        void residual::setFehat( const floatVector &Fehat ){
            /*!
             * Set the isochoric part of the elastic deformation gradient
             * 
             * \param &Fehat: The isochoric part of the elastic deformation gradient
             */

            _Fehat.second = Fehat;

            _Fehat.first = true;

            addIterationData( &_Fehat );

        }

        void residual::setPreviousJe( const floatType &previousJe ){
            /*!
             * Set the previous elastic Jacobian of deformation
             * 
             * \param &previousJe: The previous elastic Jacobian of deformation
             */

            _previousJe.second = previousJe;

            _previousJe.first = true;

        }

        void residual::setPreviousFehat( const floatVector &previousFehat ){
            /*!
             * Set the previous isochoric part of the elastic deformation gradient
             * 
             * \param &previousFehat: The previous isochoric part of the elastic deformation gradient
             */

            _previousFehat.second = previousFehat;

            _previousFehat.first = true;

            addIterationData( &_previousFehat );

        }

        const floatType* residual::getJe( ){
            /*!
             * Get the elastic Jacobian of deformation
             */

            if ( !_Je.first ){

                ERROR_TOOLS_CATCH( decomposeElasticDeformation( ) );

            }

            return &_Je.second;

        }

        const floatVector* residual::getFehat( ){
            /*!
             * Get the isochoric part of the elastic deformation gradient
             */

            if ( !_Fehat.first ){

                ERROR_TOOLS_CATCH( decomposeElasticDeformation( ) );

            }

            return &_Fehat.second;

        }

        const floatType* residual::getPreviousJe( ){
            /*!
             * Get the previous elastic Jacobian of deformation
             */

            if ( !_previousJe.first ){

                ERROR_TOOLS_CATCH( decomposePreviousElasticDeformation( ) );

            }

            return &_previousJe.second;

        }

        const floatVector* residual::getPreviousFehat( ){
            /*!
             * Get the previous isochoric part of the elastic deformation gradient
             */

            if ( !_previousFehat.first ){

                ERROR_TOOLS_CATCH( decomposePreviousElasticDeformation( ) );

            }

            return &_previousFehat.second;

        }

        void residual::setdJedFe( ){
            /*!
             * Set the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             */

            const unsigned int* dim = hydra->getDimension( );

            setdJedFe( vectorTools::computeDDetADA( ( *hydra->getConfigurations( ) )[ 0 ], ( *dim ), ( *dim ) ) );

        }

        void residual::setdJedFe( const floatVector &dJedFe ){
            /*!
             * Set the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             * 
             * \param &dJedFe: The derivative of the elastic Jacobian of deformation
             *     w.r.t. the elastic deformation gradient
             */

            _dJedFe.second = dJedFe;

            _dJedFe.first = true;

            addIterationData( &_dJedFe );

        }

        const floatVector* residual::getdJedFe( ){
            /*!
             * Get the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient
             */

            if ( !_dJedFe.first ){

                ERROR_TOOLS_CATCH( setdJedFe( ) );

            }

            return &_dJedFe.second;

        }

        void residual::setdFehatdFe( ){
            /*!
             * Set the derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             */

            const unsigned int* dim = hydra->getDimension( );

            floatMatrix dFehatdFe = vectorTools::eye< floatType >( ( *dim ) * ( *dim ) ) * std::pow( ( *getJe( ) ), -1. / 3 );

            dFehatdFe -= vectorTools::dyadic( ( *hydra->getConfigurations( ) )[ 0 ], *getdJedFe( ) ) * std::pow( ( *getJe( ) ), -4. / 3 ) / 3.;

            setdFehatdFe( dFehatdFe );

        }

        void residual::setdFehatdFe( const floatMatrix &dFehatdFe ){
            /*!
             * Set the derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             * 
             * \param &dFehatdFe: The derivative of the isochoric part of the elastic
             *     deformation gradient the elastic deformation gradient
             */

            _dFehatdFe.second = dFehatdFe;

            _dFehatdFe.first = true;

            addIterationData( &_dFehatdFe );

        }

        const floatMatrix* residual::getdFehatdFe( ){
            /*!
             * Get the derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient
             */

            if ( !_dFehatdFe.first ){

                ERROR_TOOLS_CATCH( setdFehatdFe( ) );

            }

            return &_dFehatdFe.second;

        }

        void residual::decomposeStateVariableVector( floatVector &volumetricStateVariables,
                                                     floatVector &isochoricStateVariables ){
            /*!
             * Decompose the state variable vector into parameters associated with the
             * volumetric and isochoric viscoelasticity
             */

            const unsigned int *dim = hydra->getDimension( );

            unsigned int lb = *getViscoelasticISVLowerIndex( );
            unsigned int ub = lb + *getNumVolumetricViscousTerms( );

            volumetricStateVariables = floatVector( hydra->getAdditionalStateVariables( )->begin( ) + lb,
                                                    hydra->getAdditionalStateVariables( )->begin( ) + ub );

            lb = ub;

            ub = lb + ( *dim ) * ( *dim ) * *getNumIsochoricViscousTerms( );

            isochoricStateVariables = floatVector( hydra->getAdditionalStateVariables( )->begin( ) + lb,
                                                   hydra->getAdditionalStateVariables( )->begin( ) + ub );

        }

        void residual::setNumStateVariables( const unsigned int numStateVariables ){
            /*!
             * Set the number of state variables
             * 
             * \param &numStateVariables: The number of state variables required
             */

            _numStateVariables = numStateVariables;

        }

        floatType residual::computeRateMultiplier( const floatVector &variables,
                                                   const floatVector &parameters ){
            /*!
             * Compute the value of the rate multiplier this implementation uses the
             * WLF function of temperature to compute an increase in evolution rate for
             * high temperatures and a decrease for low temperatures.
             * 
             * \param &variables: The incoming variables { temperature }
             * \param &parameters: The incoming parameters see constitutiveTools::WLF
             */

            if ( variables.size( ) != 1 ){

                ERROR_TOOLS_CATCH( throw std::runtime_error( "The incoming variables must have a size 1" ) );

            }

            floatType invRM;

            ERROR_TOOLS_CATCH_NODE_POINTER( constitutiveTools::WLF( variables[ 0 ], parameters, invRM ) );

            return 1 / invRM;

        }

        floatVector residual::computedRateMultiplierdVariables( const floatVector &variables,
                                                                const floatVector &parameters ){
            /*!
             * Compute the value of the derivative of the rate multiplier with respect to
             * the variables. This implementation uses the WLF function of temperature to
             * compute an increase in evolution rate for high temperatures and a decrease
             * for low temperatures.
             * 
             * \param &variables: The incoming variables { temperature }
             * \param &parameters: The incoming parameters see constitutiveTools::WLF
             */

            if ( variables.size( ) != 1 ){

                ERROR_TOOLS_CATCH( throw std::runtime_error( "The incoming variables must have a size 1" ) );

            }

            floatType invRM;
            floatVector dinvRMdT( 1, 0 );

            ERROR_TOOLS_CATCH_NODE_POINTER( constitutiveTools::WLF( variables[ 0 ], parameters, invRM, dinvRMdT[ 0 ] ) );

            return -1 / ( invRM * invRM ) * dinvRMdT;

        }

        void residual::setVolumetricRateMultiplier( ){
            /*!
             * Set the value of the volumetric rate multiplier
             */

            floatType rateMultiplier;

            ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getTemperature( ) },
                                                                       *getVolumetricTemperatureParameters( ) ) );

            setVolumetricRateMultiplier( rateMultiplier );

        }

        void residual::setPreviousVolumetricRateMultiplier( ){
            /*!
             * Set the previous value of the volumetric rate multiplier
             */

            floatType rateMultiplier;

            ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getPreviousTemperature( ) },
                                                                       *getVolumetricTemperatureParameters( ) ) );

            setPreviousVolumetricRateMultiplier( rateMultiplier );

        }

        void residual::setIsochoricRateMultiplier( ){
            /*!
             * Set the value of the isochoric rate multiplier
             */

            floatType rateMultiplier;

            ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getTemperature( ) },
                                                                       *getIsochoricTemperatureParameters( ) ) );

            setIsochoricRateMultiplier( rateMultiplier );

        }

        void residual::setPreviousIsochoricRateMultiplier( ){
            /*!
             * Set the previous value of the isochoric rate multiplier
             */

            floatType rateMultiplier;

            ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getPreviousTemperature( ) },
                                                                       *getIsochoricTemperatureParameters( ) ) );

            setPreviousIsochoricRateMultiplier( rateMultiplier );

        }

        void residual::setdVolumetricRateMultiplierdT( ){
            /*!
             * Set the value of the derivative of the volumetric rate multiplier
             * with respect to the temperature
             */

            floatType dRateMultiplierdT;

            ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getTemperature( ) },
                                                                                       *getVolumetricTemperatureParameters( ) )[ 0 ] );

            setdVolumetricRateMultiplierdT( dRateMultiplierdT );

        }

        void residual::setdPreviousVolumetricRateMultiplierdPreviousT( ){
            /*!
             * Set the previous value of the derivative of the volumetric rate multiplier
             * with respect to the temperature
             */

            floatType dRateMultiplierdT;

            ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getPreviousTemperature( ) },
                                                                                       *getVolumetricTemperatureParameters( ) )[ 0 ] );

            setdPreviousVolumetricRateMultiplierdPreviousT( dRateMultiplierdT );

        }

        void residual::setdIsochoricRateMultiplierdT( ){
            /*!
             * Set the value of the derivative of the isochoric rate multiplier
             * with respect to the temperature
             */

            floatType dRateMultiplierdT;

            ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getTemperature( ) },
                                                                                       *getIsochoricTemperatureParameters( ) )[ 0 ] );

            setdIsochoricRateMultiplierdT( dRateMultiplierdT );

        }

        void residual::setdPreviousIsochoricRateMultiplierdPreviousT( ){
            /*!
             * Set the previous value of the isochoric rate multiplier
             */

            floatType dRateMultiplierdT;

            ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getPreviousTemperature( ) },
                                                                                       *getIsochoricTemperatureParameters( ) )[ 0 ] );

            setdPreviousIsochoricRateMultiplierdPreviousT( dRateMultiplierdT );

        }

        void residual::setVolumetricRateMultiplier( const floatType &rateMultiplier ){
            /*!
             * Set the volumetric rate multiplier for the viscous term evolution
             * 
             * \param &rateMultiplier: The value of the volumetric rate multiplier
             */

            _volumetricRateMultiplier.second = rateMultiplier;

            _volumetricRateMultiplier.first = true;

            addIterationData( &_volumetricRateMultiplier );

        }

        void residual::setPreviousVolumetricRateMultiplier( const floatType &previousRateMultiplier ){
            /*!
             * Set the previous volumetric rate multiplier for the viscous term evolution
             * 
             * \param &previousRateMultiplier: The value of the previous volumetric rate multiplier
             */

            _previousVolumetricRateMultiplier.second = previousRateMultiplier;

            _previousVolumetricRateMultiplier.first = true;

        }

        void residual::setIsochoricRateMultiplier( const floatType &rateMultiplier ){
            /*!
             * Set the isochoric rate multiplier for the viscous term evolution
             * 
             * \param &rateMultiplier: The value of the isochoric rate multiplier
             */

            _isochoricRateMultiplier.second = rateMultiplier;

            _isochoricRateMultiplier.first = true;

            addIterationData( &_isochoricRateMultiplier );

        }

        void residual::setPreviousIsochoricRateMultiplier( const floatType &previousRateMultiplier ){
            /*!
             * Set the previous isochoric rate multiplier for the viscous term evolution
             * 
             * \param &previousRateMultiplier: The value of the previous isochoric rate multiplier
             */

            _previousIsochoricRateMultiplier.second = previousRateMultiplier;

            _previousIsochoricRateMultiplier.first = true;

        }

        void residual::setdVolumetricRateMultiplierdT( const floatType &dRateMultiplierdT ){
            /*!
             * Set the derivative of the volumetric rate multiplier for the viscous
             * term evolution w.r.t. the temperature
             * 
             * \param &dRateMultiplierdT: The value of the derivative of the
             *     volumetric rate multiplier w.r.t. T
             */

            _dVolumetricRateMultiplierdT.second = dRateMultiplierdT;

            _dVolumetricRateMultiplierdT.first = true;

            addIterationData( &_dVolumetricRateMultiplierdT );

        }

        void residual::setdPreviousVolumetricRateMultiplierdPreviousT( const floatType &dPreviousRateMultiplierdPreviousT ){
            /*!
             * Set the derivative of the previous volumetric rate multiplier for the
             * viscous term evolution w.r.t. the previous temperature.
             * 
             * \param &dPreviousRateMultiplierdPreviousT: The value of the derivative of the
             *     previous volumetric rate multiplier w.r.t. the previous temperature
             */

            _dPreviousVolumetricRateMultiplierdPreviousT.second = dPreviousRateMultiplierdPreviousT;

            _dPreviousVolumetricRateMultiplierdPreviousT.first = true;

        }

        void residual::setdIsochoricRateMultiplierdT( const floatType &dRateMultiplierdT ){
            /*!
             * Set the derivative of the isochoric rate multiplier for the viscous
             * term evolution w.r.t. the temperature
             * 
             * \param &dRateMultiplierdT: The value of the derivative of the isochoric
             *     rate multiplier w.r.t. the temperature
             */

            _dIsochoricRateMultiplierdT.second = dRateMultiplierdT;

            _dIsochoricRateMultiplierdT.first = true;

            addIterationData( &_dIsochoricRateMultiplierdT );

        }

        void residual::setdPreviousIsochoricRateMultiplierdPreviousT( const floatType &dPreviousRateMultiplierdPreviousT ){
            /*!
             * Set the derivative of the previous isochoric rate multiplier for the
             * viscous term evolution w.r.t. the previous temperature
             * 
             * \param &dPreviousRateMultiplierdPreviousT: The value of the derivative
             *     of the previous isochoric rate multiplier w.r.t. the previous temperature
             */

            _dPreviousIsochoricRateMultiplierdPreviousT.second = dPreviousRateMultiplierdPreviousT;

            _dPreviousIsochoricRateMultiplierdPreviousT.first = true;

        }

        const floatType* residual::getVolumetricRateMultiplier( ){
            /*!
             * Get the current value of the volumetric rate multiplier
             */

            if ( !_volumetricRateMultiplier.first ){

                ERROR_TOOLS_CATCH( setVolumetricRateMultiplier( ) );

            }

            return &_volumetricRateMultiplier.second;

        }

        const floatType* residual::getPreviousVolumetricRateMultiplier( ){
            /*!
             * Get the current value of the previous volumetric rate multiplier
             */

            if ( !_previousVolumetricRateMultiplier.first ){

                ERROR_TOOLS_CATCH( setPreviousVolumetricRateMultiplier( ) );

            }

            return &_previousVolumetricRateMultiplier.second;

        }

        const floatType* residual::getIsochoricRateMultiplier( ){
            /*!
             * Get the current value of the isochoric rate multiplier
             */

            if ( !_isochoricRateMultiplier.first ){

                ERROR_TOOLS_CATCH( setIsochoricRateMultiplier( ) );

            }

            return &_isochoricRateMultiplier.second;

        }

        const floatType* residual::getPreviousIsochoricRateMultiplier( ){
            /*!
             * Get the current value of the previous isochoric rate multiplier
             */

            if ( !_previousIsochoricRateMultiplier.first ){

                ERROR_TOOLS_CATCH( setPreviousIsochoricRateMultiplier( ) );

            }

            return &_previousIsochoricRateMultiplier.second;

        }

        const floatType* residual::getdVolumetricRateMultiplierdT( ){
            /*!
             * Get the current value of the derivative of the 
             * volumetric rate multiplier w.r.t. temperature
             */

            if ( !_dVolumetricRateMultiplierdT.first ){

                ERROR_TOOLS_CATCH( setdVolumetricRateMultiplierdT( ) );

            }

            return &_dVolumetricRateMultiplierdT.second;

        }

        const floatType* residual::getdPreviousVolumetricRateMultiplierdPreviousT( ){
            /*!
             * Get the current value of the derivative of the previous volumetric
             * rate multiplier w.r.t. the previous temperature
             */

            if ( !_dPreviousVolumetricRateMultiplierdPreviousT.first ){

                ERROR_TOOLS_CATCH( setdPreviousVolumetricRateMultiplierdPreviousT( ) );

            }

            return &_dPreviousVolumetricRateMultiplierdPreviousT.second;

        }

        const floatType* residual::getdIsochoricRateMultiplierdT( ){
            /*!
             * Get the current value of the derivative of the isochoric rate multiplier
             * w.r.t. the temperature
             */

            if ( !_dIsochoricRateMultiplierdT.first ){

                ERROR_TOOLS_CATCH( setdIsochoricRateMultiplierdT( ) );

            }

            return &_dIsochoricRateMultiplierdT.second;

        }

        const floatType* residual::getdPreviousIsochoricRateMultiplierdPreviousT( ){
            /*!
             * Get the current value of the previous isochoric rate multiplier
             */

            if ( !_dPreviousIsochoricRateMultiplierdPreviousT.first ){

                ERROR_TOOLS_CATCH( setdPreviousIsochoricRateMultiplierdPreviousT( ) );

            }

            return &_dPreviousIsochoricRateMultiplierdPreviousT.second;

        }

        void residual::setVolumetricTemperatureParameters( const floatVector &parameters ){
            /*!
             * Set the volumetric temperature parameters
             * 
             * \param &parameters: The parameters for the volumetric temperature dependence
             */

            _volumetricTemperatureParameters = parameters;

        }

        void residual::setIsochoricTemperatureParameters( const floatVector &parameters ){
            /*!
             * Set the isochoric temperature parameters
             * 
             * \param &parameters: The parameters for the isochoric temperature dependence
             */

            _isochoricTemperatureParameters = parameters;

        }

        floatVector residual::getVolumetricViscoelasticParameters( ){
            /*!
             * Get the volumetric viscoelastic parameters prepared for stressTools::linearViscoelasticity
             */

            floatVector parameters( 1 + 2 * *getNumVolumetricViscousTerms( ), 0 );

            parameters[ 0 ] = *getKinf( );

            for ( unsigned int i = 1; i <= *getNumVolumetricViscousTerms( ); i++ ){

                parameters[ i ] = ( *getVolumetricTaus( ) )[ i - 1 ];
                parameters[ i + *getNumVolumetricViscousTerms( ) ] = ( *getVolumetricModuli( ) )[ i - 1 ];

            }

            return parameters;

        }

        floatVector residual::getIsochoricViscoelasticParameters( ){
            /*!
             * Get the isochoric viscoelastic parameters prepared for stressTools::linearViscoelasticity
             */

            floatVector parameters( 1 + 2 * *getNumIsochoricViscousTerms( ), 0 );

            parameters[ 0 ] = 2 * ( *getGinf( ) );

            for ( unsigned int i = 1; i <= *getNumIsochoricViscousTerms( ); i++ ){

                parameters[ i ] = ( *getIsochoricTaus( ) )[ i - 1 ];
                parameters[ i + *getNumIsochoricViscousTerms( ) ] = 2 * ( *getIsochoricModuli( ) )[ i - 1 ];

            }

            return parameters;

        }

        void residual::setPK2MeanStress( ){
            /*!
             * Set the mean stress for the Second Piola-Kirchhoff stress
             */

            // Compute the strain measures

            floatVector volumetricStrain = { ( *getJe( ) - 1 ) };

            floatVector previousVolumetricStrain = { ( *getPreviousJe( ) - 1 ) };

            // Get the previous state variable values

            floatVector previousVolumetricStateVariables;

            floatVector previousIsochoricStateVariables;

            ERROR_TOOLS_CATCH( decomposeStateVariableVector( previousVolumetricStateVariables,
                                                             previousIsochoricStateVariables ) );

            floatVector PK2MeanStress;

            floatVector deltaPK2MeanStress;

            floatMatrix dPK2MeanStressdJe;

            floatVector dPK2MeanStressdRateModifier;

            floatVector currentVolumetricStateVariables;

            // Compute the viscous mean stress

            ERROR_TOOLS_CATCH_NODE_POINTER( stressTools::linearViscoelasticity( *hydra->getTime( ), volumetricStrain,
                                                                                *hydra->getTime( ) - *hydra->getDeltaTime( ), previousVolumetricStrain,
                                                                                *getVolumetricRateMultiplier( ),
                                                                                *getPreviousVolumetricRateMultiplier( ),
                                                                                previousVolumetricStateVariables,
                                                                                getVolumetricViscoelasticParameters( ),
                                                                                *getIntegrationAlpha( ), deltaPK2MeanStress,
                                                                                PK2MeanStress, currentVolumetricStateVariables, dPK2MeanStressdJe, dPK2MeanStressdRateModifier ) );

            setPK2MeanStress( PK2MeanStress[ 0 ] );

            setUpdatedVolumetricViscoelasticStateVariables( currentVolumetricStateVariables );

            setdPK2MeanStressdFe( dPK2MeanStressdJe[ 0 ][ 0 ] * ( *getdJedFe( ) ) );

            setdPK2MeanStressdT( dPK2MeanStressdRateModifier[ 0 ] * ( *getdVolumetricRateMultiplierdT( ) ) );

        }

        void residual::setdPK2MeanStressdT( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the temperature
             */

            getPK2MeanStress( );

        }

        void residual::setdPK2MeanStressdFe( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the elastic deformation gradient
             */

            getPK2MeanStress( );

        }

        void residual::setPK2MeanStress( const floatType &PK2MeanStress ){
            /*!
             * Set the PK2 mean stress
             * 
             * \param &PK2MeanStress: The second Piola-Kirchhoff mean stress
             */

            _PK2MeanStress.second = PK2MeanStress;

            _PK2MeanStress.first = true;

            addIterationData( &_PK2MeanStress );

        }

        void residual::setdPK2MeanStressdT( const floatType &dPK2MeanStressdT ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the temperature
             * 
             * \param &dPK2MeanStressdT: The derivative of the second Piola-Kirchhoff
             *     mean stress w.r.t. the temperature
             */

            _dPK2MeanStressdT.second = dPK2MeanStressdT;

            _dPK2MeanStressdT.first = true;

            addIterationData( &_dPK2MeanStressdT );

        }

        void residual::setdPK2MeanStressdFe( const floatVector &dPK2MeanStressdFe ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the elastic deformation
             * gradient
             * 
             * \param &dPK2MeanStressdFe: The gradient of the second Piola-Kirchhoff
             *     mean stress w.r.t. the temperature
             */

            _dPK2MeanStressdFe.second = dPK2MeanStressdFe;

            _dPK2MeanStressdFe.first = true;

            addIterationData( &_dPK2MeanStressdFe );

        }

        const floatType* residual::getPK2MeanStress( ){
            /*!
             * Get the PK2 mean stress
             */

            if ( !_PK2MeanStress.first ){

                ERROR_TOOLS_CATCH( setPK2MeanStress( ) );

            }

            return &_PK2MeanStress.second;

        }

        const floatType* residual::getdPK2MeanStressdT( ){
            /*!
             * Get the derivative of the PK2 mean stress w.r.t. the temperature
             */

            if ( !_dPK2MeanStressdT.first ){

                ERROR_TOOLS_CATCH( setdPK2MeanStressdT( ) );

            }

            return &_dPK2MeanStressdT.second;

        }

        const floatVector* residual::getdPK2MeanStressdFe( ){
            /*!
             * Get the derivative of the PK2 mean stress w.r.t. the elastic
             * deformation gradient
             */

            if ( !_dPK2MeanStressdFe.first ){

                ERROR_TOOLS_CATCH( setdPK2MeanStressdFe( ) );

            }

            return &_dPK2MeanStressdFe.second;

        }

        void residual::setPK2IsochoricStress( ){
            /*!
             * Set the isochoric second Piola-Kirchhoff stress
             */

            // Compute the strain measures

            floatVector isochoricStrain, previousIsochoricStrain;
            floatMatrix dEehatdFehat;
            ERROR_TOOLS_CATCH_NODE_POINTER( constitutiveTools::computeGreenLagrangeStrain( *getFehat( ), isochoricStrain, dEehatdFehat ) );

            ERROR_TOOLS_CATCH_NODE_POINTER( constitutiveTools::computeGreenLagrangeStrain( *getPreviousFehat( ), previousIsochoricStrain ) );

            floatMatrix dEehatdFe = vectorTools::dot( dEehatdFehat, *getdFehatdFe( ) );

            // Get the previous state variable values

            floatVector previousVolumetricStateVariables;

            floatVector previousIsochoricStateVariables;

            ERROR_TOOLS_CATCH( decomposeStateVariableVector( previousVolumetricStateVariables,
                                                             previousIsochoricStateVariables ) );

            floatVector PK2IsochoricStress;

            floatVector deltaPK2IsochoricStress;

            floatVector currentIsochoricStateVariables;

            floatMatrix dPK2IsochoricStressdEe;

            floatVector dPK2IsochoricStressdRateMultiplier;

            // Compute the viscous isochoric stress
            ERROR_TOOLS_CATCH_NODE_POINTER( stressTools::linearViscoelasticity( *hydra->getTime( ), isochoricStrain,
                                                                                *hydra->getTime( ) - *hydra->getDeltaTime( ), previousIsochoricStrain,
                                                                                *getIsochoricRateMultiplier( ),
                                                                                *getPreviousIsochoricRateMultiplier( ),
                                                                                previousIsochoricStateVariables,
                                                                                getIsochoricViscoelasticParameters( ),
                                                                                *getIntegrationAlpha( ), deltaPK2IsochoricStress,
                                                                                PK2IsochoricStress, currentIsochoricStateVariables, dPK2IsochoricStressdEe, dPK2IsochoricStressdRateMultiplier ) );

            setPK2IsochoricStress( PK2IsochoricStress );

            setUpdatedIsochoricViscoelasticStateVariables( currentIsochoricStateVariables );

            setdPK2IsochoricStressdFe( vectorTools::dot( dPK2IsochoricStressdEe, dEehatdFe ) );

            setdPK2IsochoricStressdT( dPK2IsochoricStressdRateMultiplier * ( *getdIsochoricRateMultiplierdT( ) ) );

        }

        void residual::setdPK2IsochoricStressdFe( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the elastic
             * deformation gradient
             */

            getPK2IsochoricStress( );

        }

        void residual::setdPK2IsochoricStressdT( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the temperature
             */

            getPK2IsochoricStress( );

        }

        void residual::setdPK2IsochoricStressdFe( const floatMatrix &dPK2IsochoricStressdFe ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the elastic
             * deformation gradient
             * 
             * \param &dPK2IsochoricStressdFe: The derivative of the PK2 stress w.r.t.
             * the elastic deformation gradient
             */

            _dPK2IsochoricStressdFe.second = dPK2IsochoricStressdFe;

            _dPK2IsochoricStressdFe.first = true;

            addIterationData( &_dPK2IsochoricStressdFe );

        }

        void residual::setdPK2IsochoricStressdT( const floatVector &dPK2IsochoricStressdT ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the temperature
             * 
             * \param &dPK2IsochoricStressdFe: The derivative of the PK2 stress w.r.t.
             * the temperature
             */

            _dPK2IsochoricStressdT.second = dPK2IsochoricStressdT;

            _dPK2IsochoricStressdT.first = true;

            addIterationData( &_dPK2IsochoricStressdT );

        }

        const floatMatrix* residual::getdPK2IsochoricStressdFe( ){
            /*!
             * Get the derivative of the isochoric PK2 stress w.r.t. the elastic
             * deformation gradient
             */

            if ( !_dPK2IsochoricStressdFe.first ){

                ERROR_TOOLS_CATCH( setdPK2IsochoricStressdFe( ) );

            }

            return &_dPK2IsochoricStressdFe.second;

        }

        const floatVector* residual::getdPK2IsochoricStressdT( ){
            /*!
             * Get the derivative of the isochoric PK2 stress w.r.t. the temperature
             */

            if ( !_dPK2IsochoricStressdT.first ){

                ERROR_TOOLS_CATCH( setdPK2IsochoricStressdT( ) );

            }

            return &_dPK2IsochoricStressdT.second;

        }

        void residual::setUpdatedVolumetricViscoelasticStateVariables( ){
            /*!
             * Set the updated values of the volumetric viscoelastic state variables
             */

           getPK2MeanStress( );

        }

        void residual::setUpdatedVolumetricViscoelasticStateVariables( floatVector &volumetricViscoelasticStateVariables ){
            /*!
             * Set the updated values of the volumetric viscoelastic state variables
             * 
             * \param &volumetricViscoelasticStateVariables: The updated volumetric viscoelastic state variables
             */

            _volumetricViscoelasticStateVariables.second = volumetricViscoelasticStateVariables;

            _volumetricViscoelasticStateVariables.first = true;

            addIterationData( &_volumetricViscoelasticStateVariables );

        }

        const floatVector* residual::getUpdatedVolumetricViscoelasticStateVariables( ){
            /*!
             * Get the updated values of the volumetric viscoelastic state variables
             */

            if ( !_volumetricViscoelasticStateVariables.first ){

                ERROR_TOOLS_CATCH( setUpdatedVolumetricViscoelasticStateVariables( ) );

            }

            return &_volumetricViscoelasticStateVariables.second;

        }

        void residual::setUpdatedIsochoricViscoelasticStateVariables( ){
            /*!
             * Set the updated values of the isochoric viscoelastic state variables
             */

           getPK2IsochoricStress( );

        }

        void residual::setUpdatedIsochoricViscoelasticStateVariables( floatVector &isochoricViscoelasticStateVariables ){
            /*!
             * Set the updated values of the isochoric viscoelastic state variables
             * 
             * \param &isochoricViscoelasticStateVariables: The updated isochoric viscoelastic state variables
             */

            _isochoricViscoelasticStateVariables.second = isochoricViscoelasticStateVariables;

            _isochoricViscoelasticStateVariables.first = true;

            addIterationData( &_isochoricViscoelasticStateVariables );

        }

        const floatVector* residual::getUpdatedIsochoricViscoelasticStateVariables( ){
            /*!
             * Get the updated values of the isochoric viscoelastic state variables
             */

            if ( !_isochoricViscoelasticStateVariables.first ){

                ERROR_TOOLS_CATCH( setUpdatedIsochoricViscoelasticStateVariables( ) );

            }

            return &_isochoricViscoelasticStateVariables.second;

        }

        void residual::setPK2IsochoricStress( const floatVector &PK2IsochoricStress ){
            /*!
             * Set the PK2 isochoric stress
             * 
             * \param &PK2MeanStress: The second Piola-Kirchhoff isochoric stress
             */

            _PK2IsochoricStress.second = PK2IsochoricStress;

            _PK2IsochoricStress.first = true;

            addIterationData( &_PK2IsochoricStress );

        }

        const floatVector* residual::getPK2IsochoricStress( ){
            /*!
             * Get the PK2 isochoric stress
             */

            if ( !_PK2IsochoricStress.first ){

                ERROR_TOOLS_CATCH( setPK2IsochoricStress( ) );

            }

            return &_PK2IsochoricStress.second;

        }

        void residual::setPK2Stress( ){
            /*!
             * Set the PK2 stress
             */

            const unsigned int* dim = hydra->getDimension( );

            floatVector eye( ( *dim ) * ( *dim ), 0 );
            vectorTools::eye( eye );

            floatVector PK2Stress = ( *getPK2IsochoricStress( ) ) + ( *getPK2MeanStress( ) ) * eye;

            setPK2Stress( PK2Stress );

        }

        void residual::setdPK2StressdFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             */

            const unsigned int* dim = hydra->getDimension( );

            floatVector eye( ( *dim ) * ( *dim ), 0 );
            vectorTools::eye( eye );

            floatMatrix dPK2StressdFe = *getdPK2IsochoricStressdFe( ) + vectorTools::dyadic( eye, *getdPK2MeanStressdFe( ) );

            setdPK2StressdFe( dPK2StressdFe );

        }

        void residual::setdPK2StressdT( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the temperature
             */

            const unsigned int* dim = hydra->getDimension( );

            floatVector eye( ( *dim ) * ( *dim ), 0 );
            vectorTools::eye( eye );

            floatVector dPK2StressdT = *getdPK2IsochoricStressdT( ) + *getdPK2MeanStressdT( ) * eye;

            setdPK2StressdT( dPK2StressdT );

        }

        void residual::setdPK2StressdT( const floatVector &dPK2StressdT ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the
             * temperature
             *
             * \param &dPK2StressdT: The derivative of the second Piola-Kirchoff stress
             *     w.r.t. the temperature
             */

            _dPK2StressdT.second = dPK2StressdT;

            _dPK2StressdT.first = true;

            addIterationData( &_dPK2StressdT );

        }

        const floatVector* residual::getdPK2StressdT( ){
            /*!
             * Get the derivative of the second Piola-Kirchhoff stress w.r.t.
             * the temperature.
             */

            if ( !_dPK2StressdT.first ){

                ERROR_TOOLS_CATCH( setdPK2StressdT( ) );

            }

            return &_dPK2StressdT.second;

        }

        void residual::setdCauchyStressdT( ){
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the temperature
             */

            floatVector dCauchyStressdT = vectorTools::dot( *getdCauchyStressdPK2Stress( ), *getdPK2StressdT( ) );

            setdCauchyStressdT( dCauchyStressdT );

        }

        void residual::setdCauchyStressdT( const floatVector &dCauchyStressdT ){
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the
             * temperature
             *
             * \param &dCauchyStressdT: The derivative of the Cauchy stress
             *     w.r.t. the temperature
             */

            _dCauchyStressdT.second = dCauchyStressdT;

            _dCauchyStressdT.first = true;

            addIterationData( &_dCauchyStressdT );

        }

        const floatVector* residual::getdCauchyStressdT( ){
            /*!
             * Get the derivative of the Cauchy stress w.r.t.
             * the temperature.
             */

            if ( !_dCauchyStressdT.first ){

                ERROR_TOOLS_CATCH( setdCauchyStressdT( ) );

            }

            return &_dCauchyStressdT.second;

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature.
             */

            setdRdT( *getdCauchyStressdT( ) );

        }

        void residual::setPK2Stress( const floatVector &PK2Stress ){
            /*!
             * Pass-through function to linearElasticity::setPK2Stress
             *
             * Required because of overloading
             *
             * \param &PK2Stress: The second Piola-Kirchhoff stress
             */

            tardigradeHydra::linearElasticity::residual::setPK2Stress( PK2Stress );
        }

        void residual::setdPK2StressdFe( const floatMatrix &dPK2StressdFe ){
            /*!
             * Pass-through function to linearElasticity::setdPK2StressdFe
             *
             * Required because of overloading
             *
             * \param &dPK2StressdFe: The derivative of the second Piola-Kirchhoff
             *     stress w.r.t. the elastic deformation gradient
             */

            tardigradeHydra::linearElasticity::residual::setdPK2StressdFe( dPK2StressdFe );
        }

        void residual::setdRdT( const floatVector &dRdT ){
            /*!
             * Pass-through function to residualBase::setdRdT
             *
             * Required because of overloading
             *
             * \param &dRdT: The derivative of the residual w.r.t. the temperature
             */

            tardigradeHydra::residualBase::setdRdT( dRdT );
        }

    }

}
