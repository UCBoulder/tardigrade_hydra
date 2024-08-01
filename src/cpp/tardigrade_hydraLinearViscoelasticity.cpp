/**
  ******************************************************************************
  * \file tardigrade_hydraLinearViscoelasticity.cpp
  ******************************************************************************
  * An implementation of linear elasticity using the hydra framework. Used as an
  * example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraLinearViscoelasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_stress_tools.h>

namespace tardigradeHydra{

    namespace linearViscoelasticity{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length 2 which defines lambda and mu.
             */
   
            const unsigned int sot_dim = hydra->getSOTDimension( );
 
            if ( parameters.size( ) < 10 ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Parameter vector is expected to have a length of at least 10 but has a length of " + std::to_string( parameters.size( ) ) ) );
    
            }

            setNumVolumetricViscousTerms( ( unsigned int )( parameters[ 0 ] + 0.5 ) );

            setNumIsochoricViscousTerms( ( unsigned int )( parameters[ 1 ] + 0.5 ) );

            setNumStateVariables( *getNumVolumetricViscousTerms( ) + sot_dim * ( *getNumIsochoricViscousTerms( ) ) );

            if ( *getNumStateVariables( ) != ( *getViscoelasticISVUpperIndex( ) - *getViscoelasticISVLowerIndex( ) ) ){

                std::string message = "The number of state variables required by the parameterization is not equal to the number of state variables indicated by the ISV bounds\n";
                message            += "   required # ISVs: " + std::to_string( *getNumStateVariables( ) ) + "\n";
                message            += "   ISV Lower Bound: " + std::to_string( *getViscoelasticISVLowerIndex( ) ) + "\n";
                message            += "   ISV UPper Bound: " + std::to_string( *getViscoelasticISVLowerIndex( ) ) + "\n";

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

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
             * \param &taus: The isochoric time constants
             */

            _isochoricTaus = taus;

        }

        void residual::decomposeElasticDeformation( ){
            /*!
             * Decompose the elastic deformation into volumetric and isochoric parts
             */

            floatType Je;

            secondOrderTensor Fehat;

            TARDIGRADE_ERROR_TOOLS_CATCH( decomposeDeformation( *get_Fe( ), Je, Fehat ) );

            set_Je( Je );

            set_Fehat( Fehat );

        }

        void residual::decomposePreviousElasticDeformation( ){
            /*!
             * Decompose the previous elastic deformation into volumetric and isochoric parts
             */

            floatType previousJe;

            secondOrderTensor previousFehat;

            TARDIGRADE_ERROR_TOOLS_CATCH( decomposeDeformation( *get_previousFe( ), previousJe, previousFehat ) );

            set_previousJe( previousJe );

            set_previousFehat( previousFehat );

        }

        void residual::decomposeDeformation( const secondOrderTensor &F, floatType &J, secondOrderTensor &Fhat ){
            /*!
             * Decompose a deformation into volumetric and isochoric parts where
             * 
             * \f$\hat{\bf{F}} = J^{-\frac{1}{3}} \bf{F} \f$
             * 
             * \param &F: The incoming deformation gradient
             * \param &J: The jacobian of deformation
             * \param &Fhat: The isochoric part of the deformation gradient
             */

            const unsigned int dim = hydra->getDimension( );

            TARDIGRADE_ERROR_TOOLS_CATCH( J = tardigradeVectorTools::determinant( F, dim, dim ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( Fhat = F / std::pow( J, 1./3 ) );

        }

        void residual::setdJedFe( const bool isPrevious ){
            /*!
             * Set the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             * 
             * \param isPrevious: Flag for if the derivative is of the current (false) or previous (true) value
             */

            const unsigned int dim = hydra->getDimension( );

            const secondOrderTensor *Fe;

            if ( isPrevious ){

                Fe = get_previousFe( );

            }
            else{

                Fe = get_Fe( );

            }

            secondOrderTensor dJedFe = tardigradeVectorTools::computeDDetADA( *Fe, dim, dim );

            if ( isPrevious ){

                set_previousdJedFe( dJedFe );

            }
            else{

                set_dJedFe( dJedFe );

            }

        }

        void residual::setdJedFe( ){
            /*!
             * Set the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             */

            setdJedFe( false );

        }

        void residual::setdFehatdFe( ){
            /*!
             * Set the derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             */

            setdFehatdFe( false );

        }

        void residual::setPreviousdJedFe( ){
            /*!
             * Set the previous derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             */

            setdJedFe( true );

        }

        void residual::setPreviousdFehatdFe( ){
            /*!
             * Set the previous derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             */

            setdFehatdFe( true );

        }

        void residual::setdFehatdFe( const bool isPrevious ){
            /*!
             * Set the derivative of the isochoric part of the elastic deformation gradient
             * w.r.t. the elastic deformation gradient.
             * 
             * \param isPrevious: Flag for if the derivative is of the current (false) or previous (true) value
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const floatType   *Je;

            const secondOrderTensor *Fe;

            const secondOrderTensor *dJedFe;

            if ( isPrevious ){

                Je     = get_previousJe( );

                Fe     = get_previousFe( );

                dJedFe = get_previousdJedFe( );

            }
            else{

                Je     = get_Je( );

                Fe     = get_Fe( );

                dJedFe = get_dJedFe( );

            }

            fourthOrderTensor dFehatdFe( sot_dim * sot_dim, 0 );
            tardigradeVectorTools::eye< floatType >( dFehatdFe );
            dFehatdFe *= std::pow( ( *Je ), -1. / 3 );

            dFehatdFe -= tardigradeVectorTools::matrixMultiply( *Fe, *dJedFe, sot_dim, 1, 1, sot_dim ) * std::pow( ( *Je ), -4. / 3 ) / 3.;

            if ( isPrevious ){

                set_previousdFehatdFe( dFehatdFe );

            }
            else{

                set_dFehatdFe( dFehatdFe );

            }

        }

        void residual::decomposeStateVariableVector( floatVector &volumetricStateVariables,
                                                     floatVector &isochoricStateVariables ){
            /*!
             * Decompose the state variable vector into parameters associated with the
             * volumetric and isochoric viscoelasticity
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            unsigned int lb = *getViscoelasticISVLowerIndex( );
            unsigned int ub = lb + *getNumVolumetricViscousTerms( );

            volumetricStateVariables = floatVector( hydra->get_additionalStateVariables( )->begin( ) + lb,
                                                    hydra->get_additionalStateVariables( )->begin( ) + ub );

            lb = ub;

            ub = lb + sot_dim * *getNumIsochoricViscousTerms( );

            isochoricStateVariables = floatVector( hydra->get_additionalStateVariables( )->begin( ) + lb,
                                                   hydra->get_additionalStateVariables( )->begin( ) + ub );

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
             * \param &parameters: The incoming parameters see tardigradeConstitutiveTools::WLF
             */

            if ( variables.size( ) != 1 ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "The incoming variables must have a size 1" ) );

            }

            floatType invRM;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::WLF( variables[ 0 ], parameters, invRM ) );

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
             * \param &parameters: The incoming parameters see tardigradeConstitutiveTools::WLF
             */

            if ( variables.size( ) != 1 ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "The incoming variables must have a size 1" ) );

            }

            floatType invRM;
            floatVector dinvRMdT( 1, 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::WLF( variables[ 0 ], parameters, invRM, dinvRMdT[ 0 ] ) );

            return -1 / ( invRM * invRM ) * dinvRMdT;

        }

        void residual::setVolumetricRateMultiplier( ){
            /*!
             * Set the value of the volumetric rate multiplier
             */

            floatType rateMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getTemperature( ) },
                                                                                  *getVolumetricTemperatureParameters( ) ) );

            set_volumetricRateMultiplier( rateMultiplier );

        }

        void residual::setPreviousVolumetricRateMultiplier( ){
            /*!
             * Set the previous value of the volumetric rate multiplier
             */

            floatType rateMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getPreviousTemperature( ) },
                                                                                  *getVolumetricTemperatureParameters( ) ) );

            set_previousVolumetricRateMultiplier( rateMultiplier );

        }

        void residual::setIsochoricRateMultiplier( ){
            /*!
             * Set the value of the isochoric rate multiplier
             */

            floatType rateMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getTemperature( ) },
                                                                                  *getIsochoricTemperatureParameters( ) ) );

            set_isochoricRateMultiplier( rateMultiplier );

        }

        void residual::setPreviousIsochoricRateMultiplier( ){
            /*!
             * Set the previous value of the isochoric rate multiplier
             */

            floatType rateMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH( rateMultiplier = computeRateMultiplier( { *hydra->getPreviousTemperature( ) },
                                                                                  *getIsochoricTemperatureParameters( ) ) );

            set_previousIsochoricRateMultiplier( rateMultiplier );

        }

        void residual::setdVolumetricRateMultiplierdT( ){
            /*!
             * Set the value of the derivative of the volumetric rate multiplier
             * with respect to the temperature
             */

            floatType dRateMultiplierdT;

            TARDIGRADE_ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getTemperature( ) },
                                                                                                *getVolumetricTemperatureParameters( ) )[ 0 ] );

            set_dVolumetricRateMultiplierdT( dRateMultiplierdT );

        }

        void residual::setdPreviousVolumetricRateMultiplierdPreviousT( ){
            /*!
             * Set the previous value of the derivative of the volumetric rate multiplier
             * with respect to the temperature
             */

            floatType dRateMultiplierdT;

            TARDIGRADE_ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getPreviousTemperature( ) },
                                                                                       *getVolumetricTemperatureParameters( ) )[ 0 ] );

            set_dPreviousVolumetricRateMultiplierdPreviousT( dRateMultiplierdT );

        }

        void residual::setdIsochoricRateMultiplierdT( ){
            /*!
             * Set the value of the derivative of the isochoric rate multiplier
             * with respect to the temperature
             */

            floatType dRateMultiplierdT;

            TARDIGRADE_ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getTemperature( ) },
                                                                                       *getIsochoricTemperatureParameters( ) )[ 0 ] );

            set_dIsochoricRateMultiplierdT( dRateMultiplierdT );

        }

        void residual::setdPreviousIsochoricRateMultiplierdPreviousT( ){
            /*!
             * Set the previous value of the isochoric rate multiplier
             */

            floatType dRateMultiplierdT;

            TARDIGRADE_ERROR_TOOLS_CATCH( dRateMultiplierdT = computedRateMultiplierdVariables( { *hydra->getPreviousTemperature( ) },
                                                                                       *getIsochoricTemperatureParameters( ) )[ 0 ] );

            set_dPreviousIsochoricRateMultiplierdPreviousT( dRateMultiplierdT );

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
             * Get the volumetric viscoelastic parameters prepared for tardigradeStressTools::linearViscoelasticity
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
             * Get the isochoric viscoelastic parameters prepared for tardigradeStressTools::linearViscoelasticity
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

            setPK2MeanStress( false );

        }

        void residual::setPreviousPK2MeanStress( ){
            /*!
             * Set the previous mean stress for the Second Piola-Kirchhoff stress
             */

            setPK2MeanStress( true );

        }

        void residual::setPK2MeanStress( const bool isPrevious ){
            /*!
             * Set the mean stress for the Second Piola-Kirchhoff stress
             * 
             * \param &isPrevious: Flag for if the previous (true) or current (false) stress should be calculated
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const floatType *Je;

            const secondOrderTensor *dJedFe;

            const floatType *previousJe = get_previousJe( );

            floatType time;

            floatType previousTime = *hydra->getTime( ) - *hydra->getDeltaTime( );

            const floatType *volumetricRateMultiplier;

            const floatType *dVolumetricRateMultiplierdT;

            const floatType *previousVolumetricRateMultiplier = get_previousVolumetricRateMultiplier( );

            if ( isPrevious ){

                time                        = previousTime;

                Je                          = get_previousJe( );

                dJedFe                      = get_previousdJedFe( );

                volumetricRateMultiplier    = get_previousVolumetricRateMultiplier( );

                dVolumetricRateMultiplierdT = get_dPreviousVolumetricRateMultiplierdPreviousT( );

            }
            else{

                time                     = *hydra->getTime( );

                Je                       = get_Je( );

                dJedFe                   = get_dJedFe( );

                volumetricRateMultiplier = get_volumetricRateMultiplier( );

                dVolumetricRateMultiplierdT = get_dVolumetricRateMultiplierdT( );
            }

            // Compute the strain measures

            floatVector volumetricStrain = { ( *Je - 1 ) };

            floatVector previousVolumetricStrain = { ( *previousJe - 1 ) };

            // Get the previous state variable values

            floatVector previousVolumetricStateVariables;

            floatVector previousIsochoricStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( decomposeStateVariableVector( previousVolumetricStateVariables,
                                                                        previousIsochoricStateVariables ) );

            floatVector PK2MeanStress;

            floatVector deltaPK2MeanStress;

            floatMatrix _dPK2MeanStressdJe;

            floatVector dPK2MeanStressdJe;

            floatVector dPK2MeanStressdRateModifier;

            floatMatrix _dPK2MeanStressdPreviousJe;

            floatVector dPK2MeanStressdPreviousJe;

            floatVector dPK2MeanStressdPreviousRateModifier;

            floatMatrix _dPK2MeanStressdPreviousVolumetricISVs;

            floatVector dPK2MeanStressdPreviousVolumetricISVs;

            floatMatrix _dISVsdJe;

            floatVector dISVsdJe;

            floatVector dISVsdRateModifier;

            floatMatrix _dISVsdPreviousJe;

            floatVector dISVsdPreviousJe;

            floatVector dISVsdPreviousRateModifier;

            floatMatrix _dISVsdPreviousVolumetricISVs;

            floatVector dISVsdPreviousVolumetricISVs;

            floatVector currentVolumetricStateVariables;

            // Compute the viscous mean stress

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::linearViscoelasticity( time,         volumetricStrain,
                                                                                                     previousTime, previousVolumetricStrain,
                                                                                                     *volumetricRateMultiplier,
                                                                                                     *previousVolumetricRateMultiplier,
                                                                                                     previousVolumetricStateVariables,
                                                                                                     getVolumetricViscoelasticParameters( ),
                                                                                                     *getIntegrationAlpha( ), deltaPK2MeanStress,
                                                                                                     PK2MeanStress, currentVolumetricStateVariables,
                                                                                                     _dPK2MeanStressdJe, dPK2MeanStressdRateModifier,
                                                                                                     _dPK2MeanStressdPreviousJe, dPK2MeanStressdPreviousRateModifier,
                                                                                                     _dPK2MeanStressdPreviousVolumetricISVs,
                                                                                                     _dISVsdJe, dISVsdRateModifier,
                                                                                                     _dISVsdPreviousJe, dISVsdPreviousRateModifier,
                                                                                                     _dISVsdPreviousVolumetricISVs ) );

            dPK2MeanStressdJe                     = tardigradeVectorTools::appendVectors( _dPK2MeanStressdJe );

            dPK2MeanStressdPreviousJe             = tardigradeVectorTools::appendVectors( _dPK2MeanStressdPreviousJe );

            dPK2MeanStressdPreviousVolumetricISVs = tardigradeVectorTools::appendVectors( _dPK2MeanStressdPreviousVolumetricISVs );

            dISVsdJe                              = tardigradeVectorTools::appendVectors( _dISVsdJe );

            dISVsdPreviousJe                      = tardigradeVectorTools::appendVectors( _dISVsdPreviousJe );

            dISVsdPreviousVolumetricISVs          = tardigradeVectorTools::appendVectors( _dISVsdPreviousVolumetricISVs );

            if( isPrevious ){

                set_previousPK2MeanStress( PK2MeanStress[ 0 ] );
    
                set_previousdPK2MeanStressdFe( dPK2MeanStressdJe[ 0 ] * ( *dJedFe ) );
    
                set_previousdPK2MeanStressdT( dPK2MeanStressdRateModifier[ 0 ] * ( *dVolumetricRateMultiplierdT ) );

            }
            else{

                set_PK2MeanStress( PK2MeanStress[ 0 ] );
    
                set_volumetricViscoelasticStateVariables( currentVolumetricStateVariables );
    
                set_dPK2MeanStressdFe( dPK2MeanStressdJe[ 0 ] * ( *dJedFe ) );

                set_dPK2MeanStressdT( dPK2MeanStressdRateModifier[ 0 ] * ( *dVolumetricRateMultiplierdT ) );

                set_dPK2MeanStressdPreviousFe( dPK2MeanStressdPreviousJe[ 0 ] * ( *get_previousdJedFe( ) ) );

                set_dPK2MeanStressdPreviousT( dPK2MeanStressdPreviousRateModifier[ 0 ] * ( *get_dPreviousVolumetricRateMultiplierdPreviousT( ) ) );


                floatVector dPK2MeanStressdPreviousISVs( previousVolumetricStateVariables.size( ) + previousIsochoricStateVariables.size( ), 0 );

                const unsigned int vol_isvs_size = previousVolumetricStateVariables.size( );

                const unsigned int iso_isvs_size = previousIsochoricStateVariables.size( );

                floatVector dVolumetricISVsdPreviousISVs( vol_isvs_size * ( vol_isvs_size + iso_isvs_size ), 0 );

                for ( unsigned int i = 0; i < vol_isvs_size; i++ ){

                    dPK2MeanStressdPreviousISVs[ i ] = dPK2MeanStressdPreviousVolumetricISVs[ i ];

                    for ( unsigned int j = 0; j < vol_isvs_size; j++ ){

                        dVolumetricISVsdPreviousISVs[ ( vol_isvs_size + iso_isvs_size ) * i + j ] = dISVsdPreviousVolumetricISVs[ ( vol_isvs_size ) * i + j ];

                    }

                }

                set_dPK2MeanStressdPreviousISVs( dPK2MeanStressdPreviousISVs );

                set_dVolumetricISVsdFe( tardigradeVectorTools::matrixMultiply( dISVsdJe,  *dJedFe, dISVsdJe.size( ), 1, 1, sot_dim ) );

                set_dVolumetricISVsdT( dISVsdRateModifier * ( *dVolumetricRateMultiplierdT ) );

                set_dVolumetricISVsdPreviousFe( tardigradeVectorTools::matrixMultiply( dISVsdPreviousJe, *get_previousdJedFe( ), dISVsdPreviousJe.size( ), 1, 1, sot_dim ) );

                set_dVolumetricISVsdPreviousT( dISVsdPreviousRateModifier * ( *get_dPreviousVolumetricRateMultiplierdPreviousT( ) ) );

                set_dVolumetricISVsdPreviousISVs( dVolumetricISVsdPreviousISVs );

            }

        }

        void residual::setdPK2MeanStressdT( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the temperature
             */

            get_PK2MeanStress( );

        }

        void residual::setdPK2MeanStressdFe( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the elastic deformation gradient
             */

            get_PK2MeanStress( );

        }

        void residual::setdPK2MeanStressdPreviousT( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the previous temperature
             */

            get_PK2MeanStress( );

        }

        void residual::setdPK2MeanStressdPreviousFe( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the previous elastic deformation gradient
             */

            get_PK2MeanStress( );

        }

        void residual::setdPK2MeanStressdPreviousISVs( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the previous state variables
             */

            get_PK2MeanStress( );

        }

        void residual::setdVolumetricISVsdT( ){
            /*!
             * Set the derivative of the volumetric isvs w.r.t. the temperature
             */

            get_PK2MeanStress( );

        }

        void residual::setdVolumetricISVsdFe( ){
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the elastic deformation gradient
             */

            get_PK2MeanStress( );

        }

        void residual::setdVolumetricISVsdPreviousT( ){
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the previous temperature
             */

            get_PK2MeanStress( );

        }

        void residual::setdVolumetricISVsdPreviousFe( ){
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the previous elastic deformation gradient
             */

            get_PK2MeanStress( );

        }

        void residual::setdVolumetricISVsdPreviousISVs( ){
            /*!
             * Set the derivative of the volumetric ISVs w.r.t. the previous state variables
             */

            get_PK2MeanStress( );

        }

        void residual::setPreviousdPK2MeanStressdT( ){
            /*!
             * Set the previous derivative of the PK2 mean stress w.r.t. the temperature
             */

            get_previousPK2MeanStress( );

        }

        void residual::setPreviousdPK2MeanStressdFe( ){
            /*!
             * Set the derivative of the PK2 mean stress w.r.t. the elastic deformation gradient
             */

            get_previousPK2MeanStress( );

        }

        void residual::setPK2IsochoricStress( ){
            /*!
             * Set the isochoric second Piola-Kirchhoff stress
             */

            setPK2IsochoricStress( false );

        }

        void residual::setPreviousPK2IsochoricStress( ){
            /*!
             * Set the previous isochoric second Piola-Kirchhoff stress
             */

            setPK2IsochoricStress( true );

        }

        void residual::setPK2IsochoricStress( const bool isPrevious ){
            /*!
             * Set the isochoric second Piola-Kirchhoff stress
             * 
             * \param &isPrevious: Flag for if the previous (true) or current (false) stress should be calculated
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const secondOrderTensor *Fehat;

            const secondOrderTensor *previousFehat = get_previousFehat( );

            const fourthOrderTensor *dFehatdFe;

            floatType time;

            floatType previousTime = *hydra->getTime( ) - *hydra->getDeltaTime( );

            const floatType *isochoricRateMultiplier;

            const floatType *dIsochoricRateMultiplierdT;

            const floatType *previousIsochoricRateMultiplier = get_previousIsochoricRateMultiplier( );

            setDataStorageBase< secondOrderTensor > PK2IsochoricStress;

            setDataStorageBase< fourthOrderTensor > dPK2IsochoricStressdFe;

            setDataStorageBase< secondOrderTensor > dPK2IsochoricStressdT;

            if ( isPrevious ){

                time                       = previousTime;

                Fehat                      = get_previousFehat( );

                dFehatdFe                  = get_previousdFehatdFe( );

                isochoricRateMultiplier    = get_previousIsochoricRateMultiplier( );

                dIsochoricRateMultiplierdT = get_dPreviousIsochoricRateMultiplierdPreviousT( );

                PK2IsochoricStress         = get_setDataStorage_previousPK2IsochoricStress( );

                dPK2IsochoricStressdFe     = get_setDataStorage_previousdPK2IsochoricStressdFe( );

                dPK2IsochoricStressdT      = get_setDataStorage_previousdPK2IsochoricStressdT( );

            }
            else{

                time                       = *hydra->getTime( );

                Fehat                      = get_Fehat( );

                dFehatdFe                  = get_dFehatdFe( );

                isochoricRateMultiplier    = get_isochoricRateMultiplier( );

                dIsochoricRateMultiplierdT = get_dIsochoricRateMultiplierdT( );

                PK2IsochoricStress         = get_setDataStorage_PK2IsochoricStress( );

                dPK2IsochoricStressdFe     = get_setDataStorage_dPK2IsochoricStressdFe( );

                dPK2IsochoricStressdT      = get_setDataStorage_dPK2IsochoricStressdT( );

            }

            // Compute the strain measures

            secondOrderTensor isochoricStrain, previousIsochoricStrain;
            fourthOrderTensor dEehatdFehat, previousdEehatdFehat;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( *Fehat, isochoricStrain, dEehatdFehat ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( *previousFehat, previousIsochoricStrain, previousdEehatdFehat ) );

            fourthOrderTensor dEehatdFe = tardigradeVectorTools::matrixMultiply( dEehatdFehat, *dFehatdFe, sot_dim, sot_dim, sot_dim, sot_dim );

            // Get the previous state variable values

            floatVector previousVolumetricStateVariables;

            floatVector previousIsochoricStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( decomposeStateVariableVector( previousVolumetricStateVariables,
                                                             previousIsochoricStateVariables ) );

            secondOrderTensor deltaPK2IsochoricStress;

            floatVector currentIsochoricStateVariables;

            floatMatrix _dPK2IsochoricStressdEe;

            fourthOrderTensor dPK2IsochoricStressdEe;

            secondOrderTensor dPK2IsochoricStressdRateMultiplier;

            floatMatrix _dPK2IsochoricStressdPreviousEe;

            fourthOrderTensor dPK2IsochoricStressdPreviousEe;

            secondOrderTensor dPK2IsochoricStressdPreviousRateMultiplier;

            floatMatrix _dPK2IsochoricStressdPreviousIsochoricISVs;

            floatVector dPK2IsochoricStressdPreviousIsochoricISVs;

            floatMatrix _dISVsdEe;

            floatVector dISVsdEe;

            floatVector dISVsdRateMultiplier;

            floatMatrix _dISVsdPreviousEe;

            floatVector dISVsdPreviousEe;

            floatVector dISVsdPreviousRateMultiplier;

            floatMatrix _dISVsdPreviousIsochoricISVs;

            floatVector dISVsdPreviousIsochoricISVs;

            // Compute the viscous isochoric stress
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::linearViscoelasticity( time, isochoricStrain,
                                                                                                     previousTime, previousIsochoricStrain,
                                                                                                     *isochoricRateMultiplier,
                                                                                                     *previousIsochoricRateMultiplier,
                                                                                                     previousIsochoricStateVariables,
                                                                                                     getIsochoricViscoelasticParameters( ),
                                                                                                     *getIntegrationAlpha( ),
                                                                                                     deltaPK2IsochoricStress, *PK2IsochoricStress.value, currentIsochoricStateVariables,
                                                                                                     _dPK2IsochoricStressdEe,         dPK2IsochoricStressdRateMultiplier,
                                                                                                     _dPK2IsochoricStressdPreviousEe, dPK2IsochoricStressdPreviousRateMultiplier,
                                                                                                     _dPK2IsochoricStressdPreviousIsochoricISVs,
                                                                                                     _dISVsdEe,                       dISVsdRateMultiplier,
                                                                                                     _dISVsdPreviousEe,               dISVsdPreviousRateMultiplier,
                                                                                                     _dISVsdPreviousIsochoricISVs ) );

            dPK2IsochoricStressdEe                    = tardigradeVectorTools::appendVectors( _dPK2IsochoricStressdEe );

            dPK2IsochoricStressdPreviousEe            = tardigradeVectorTools::appendVectors( _dPK2IsochoricStressdPreviousEe );

            dPK2IsochoricStressdPreviousIsochoricISVs = tardigradeVectorTools::appendVectors( _dPK2IsochoricStressdPreviousIsochoricISVs );

            dISVsdEe                                  = tardigradeVectorTools::appendVectors( _dISVsdEe );

            dISVsdPreviousEe                          = tardigradeVectorTools::appendVectors( _dISVsdPreviousEe );

            dISVsdPreviousIsochoricISVs               = tardigradeVectorTools::appendVectors( _dISVsdPreviousIsochoricISVs );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2IsochoricStressdEe( dPK2IsochoricStressdEe.data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dEehatdFe( dEehatdFe.data( ), sot_dim, sot_dim );

            dPK2IsochoricStressdFe.zero( sot_dim * sot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2IsochoricStressdFe( dPK2IsochoricStressdFe.value->data( ), sot_dim, sot_dim );

            map_dPK2IsochoricStressdFe = ( map_dPK2IsochoricStressdEe * map_dEehatdFe ).eval( );

            *dPK2IsochoricStressdT.value = dPK2IsochoricStressdRateMultiplier * ( *dIsochoricRateMultiplierdT );

            if ( !isPrevious ){

                set_isochoricViscoelasticStateVariables( currentIsochoricStateVariables );

                fourthOrderTensor previousdEehatdFe = tardigradeVectorTools::matrixMultiply( previousdEehatdFehat, *get_previousdFehatdFe( ), sot_dim, sot_dim, sot_dim, sot_dim );

                set_dPK2IsochoricStressdPreviousFe( tardigradeVectorTools::matrixMultiply( dPK2IsochoricStressdPreviousEe, previousdEehatdFe, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dPK2IsochoricStressdPreviousT( dPK2IsochoricStressdPreviousRateMultiplier * ( *get_dPreviousIsochoricRateMultiplierdPreviousT( ) ) );

                const unsigned int vol_isvs_size = previousVolumetricStateVariables.size( );

                const unsigned int iso_isvs_size = previousIsochoricStateVariables.size( );

                floatVector dPK2IsochoricStressdPreviousISVs( sot_dim * ( vol_isvs_size + iso_isvs_size ), 0 );

                floatVector dIsochoricISVsdPreviousISVs( iso_isvs_size * ( vol_isvs_size + iso_isvs_size ), 0 );

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < iso_isvs_size; j++ ){

                        dPK2IsochoricStressdPreviousISVs[ ( vol_isvs_size + iso_isvs_size ) * i + j + vol_isvs_size ] = dPK2IsochoricStressdPreviousIsochoricISVs[ iso_isvs_size * i + j ];

                    }

                }

                for ( unsigned int i = 0; i < iso_isvs_size; i++ ){

                    for ( unsigned int j = 0; j < iso_isvs_size; j++ ){

                        dIsochoricISVsdPreviousISVs[ ( vol_isvs_size + iso_isvs_size ) * i + j + vol_isvs_size ] = dISVsdPreviousIsochoricISVs[ iso_isvs_size * i + j ];

                    }

                }

                set_dPK2IsochoricStressdPreviousISVs( dPK2IsochoricStressdPreviousISVs );

                set_dIsochoricISVsdFe( tardigradeVectorTools::matrixMultiply( dISVsdEe, dEehatdFe, iso_isvs_size, sot_dim, sot_dim, sot_dim ) );

                set_dIsochoricISVsdT( dISVsdRateMultiplier * ( *dIsochoricRateMultiplierdT ) );

                set_dIsochoricISVsdPreviousFe( tardigradeVectorTools::matrixMultiply( dISVsdPreviousEe, previousdEehatdFe, iso_isvs_size, sot_dim, sot_dim, sot_dim ) );

                set_dIsochoricISVsdPreviousT( dISVsdPreviousRateMultiplier * ( *get_dPreviousIsochoricRateMultiplierdPreviousT( ) ) );

                set_dIsochoricISVsdPreviousISVs( dIsochoricISVsdPreviousISVs );

            }

        }

        void residual::setdPK2IsochoricStressdFe( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the elastic
             * deformation gradient
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdPK2IsochoricStressdT( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the temperature
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdPK2IsochoricStressdPreviousFe( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the previous elastic
             * deformation gradient
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdPK2IsochoricStressdPreviousT( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the previous temperature
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdPK2IsochoricStressdPreviousISVs( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the previous state variables
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdIsochoricISVsdFe( ){
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the elastic
             * deformation gradient
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdIsochoricISVsdT( ){
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the temperature
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdIsochoricISVsdPreviousFe( ){
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the previous elastic
             * deformation gradient
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdIsochoricISVsdPreviousT( ){
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the previous temperature
             */

            get_PK2IsochoricStress( );

        }

        void residual::setdIsochoricISVsdPreviousISVs( ){
            /*!
             * Set the derivative of the isochoric isvs w.r.t. the previous state variables
             */

            get_PK2IsochoricStress( );

        }

        void residual::setPreviousdPK2IsochoricStressdFe( ){
            /*!
             * Set the previous derivative of the isochoric PK2 stress w.r.t. the elastic
             * deformation gradient
             */

            get_previousPK2IsochoricStress( );

        }

        void residual::setPreviousdPK2IsochoricStressdT( ){
            /*!
             * Set the derivative of the isochoric PK2 stress w.r.t. the temperature
             */

            get_previousPK2IsochoricStress( );

        }

        void residual::setUpdatedVolumetricViscoelasticStateVariables( ){
            /*!
             * Set the updated values of the volumetric viscoelastic state variables
             */

           get_PK2MeanStress( );

        }

        void residual::setUpdatedIsochoricViscoelasticStateVariables( ){
            /*!
             * Set the updated values of the isochoric viscoelastic state variables
             */

           get_PK2IsochoricStress( );

        }

        void residual::setCurrentAdditionalStateVariables( ){
            /*!
             * Set the updated current additional state variables
             */

            auto currentAdditionalStateVariables = get_setDataStorage_currentAdditionalStateVariables( );

            TARDIGRADE_ERROR_TOOLS_CATCH( *currentAdditionalStateVariables.value = tardigradeVectorTools::appendVectors( { *get_volumetricViscoelasticStateVariables( ),
                                                                                                                           *get_isochoricViscoelasticStateVariables( ) } ) );

        }

        void residual::setPK2Stress( const bool isPrevious ){
            /*!
             * Set the PK2 stress
             * 
             * \param isPrevious: Flag for if to compute the current (false) or previous (true) PK2 stress
             */

            constexpr unsigned int dim = 3;

            const secondOrderTensor *isochoric;

            const floatType *mean;

            setDataStorageBase< secondOrderTensor > PK2Stress;

            if ( isPrevious ){

                isochoric = get_previousPK2IsochoricStress( );

                mean      = get_previousPK2MeanStress( );

                PK2Stress = get_setDataStorage_previousPK2Stress( );

            }
            else{

                isochoric = get_PK2IsochoricStress( );

                mean      = get_PK2MeanStress( );

                PK2Stress = get_setDataStorage_PK2Stress( );

            }

            *PK2Stress.value = *isochoric;

            for ( unsigned int i = 0; i < dim; i++ ){

                ( *PK2Stress.value )[ dim * i + i ] += *mean;

            }

        }

        void residual::setdPK2StressdFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             */

            setdPK2StressdFe( false );

        }

        void residual::setPreviousdPK2StressdFe( ){
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             */

            setdPK2StressdFe( true );

        }

        void residual::setdPK2StressdFe( const bool isPrevious ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             * 
             * \param isPrevious: Flag for whether to compute the derivative of the current (false) or previous (true) stress
             */

            const unsigned int dim = hydra->getDimension( );
            const unsigned int sot_dim = dim * dim;

            const fourthOrderTensor *dIsodFe;

            const secondOrderTensor *dMeandFe;

            setDataStorageBase< fourthOrderTensor > dPK2StressdFe;

            if ( isPrevious ){

                dIsodFe  = get_previousdPK2IsochoricStressdFe( );

                dMeandFe = get_previousdPK2MeanStressdFe( );

                dPK2StressdFe = get_setDataStorage_previousdPK2StressdFe( );

            }
            else{

                dIsodFe  = get_dPK2IsochoricStressdFe( );

                dMeandFe = get_dPK2MeanStressdFe( );

                dPK2StressdFe = get_setDataStorage_dPK2StressdFe( );

            }

            *dPK2StressdFe.value = *dIsodFe;

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int jk = 0; jk < sot_dim; jk++ ){

                    ( *dPK2StressdFe.value )[ dim * dim * dim * i + dim * dim * i + jk ] += ( *dMeandFe )[ jk ];

                }

            }

            if ( !isPrevious ){

                auto dPK2StressdPreviousFe = get_setDataStorage_dPK2StressdPreviousFe( );

                *dPK2StressdPreviousFe.value = *get_dPK2IsochoricStressdPreviousFe( );

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int jk = 0; jk < sot_dim; jk++ ){

                        ( *dPK2StressdPreviousFe.value )[ dim * dim * dim * i + dim * dim * i + jk ] += ( *get_dPK2MeanStressdPreviousFe( ) )[ jk ];

                    }

                }

            }

        }

        void residual::setdPK2StressdPreviousFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the previous elastic deformation gradient
             */

            setdPK2StressdFe( false );

        }

        void residual::setdPK2StressdT( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the temperature
             */

            constexpr unsigned int dim = 3;

            auto dPK2StressdT = get_setDataStorage_dPK2StressdT( );

            *dPK2StressdT.value = *get_dPK2IsochoricStressdT( );

            for ( unsigned int i = 0; i < dim; i++ ){

                ( *dPK2StressdT.value )[ dim * i + i ]
                    += ( *get_dPK2MeanStressdT( ) );

            }

        }

        void residual::setdPK2StressdPreviousT( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the previous temperature
             */

            constexpr unsigned int dim = 3;

            auto dPK2StressdPreviousT = get_setDataStorage_dPK2StressdPreviousT( );

            *dPK2StressdPreviousT.value = *get_dPK2IsochoricStressdPreviousT( );

            for ( unsigned int i = 0; i < dim; i++ ){

                ( *dPK2StressdPreviousT.value )[ dim * i + i ]
                    += ( *get_dPK2MeanStressdPreviousT( ) );

            }

        }

        void residual::setPreviousdPK2StressdT( ){
            /*!
             * Set the prevoius derivative of the second Piola-Kirchhoff stress w.r.t. the temperature
             */

            constexpr unsigned int dim = 3;

            auto previousdPK2StressdT = get_setDataStorage_previousdPK2StressdT( );

            *previousdPK2StressdT.value = *get_previousdPK2IsochoricStressdT( );

            for ( unsigned int i = 0; i < dim; i++ ){

                ( *previousdPK2StressdT.value )[ dim * i + i ]
                    += ( *get_previousdPK2MeanStressdT( ) );

            }

        }

        void residual::setdPK2StressdPreviousISVs( ){
            /*!
             * Set the prevoius derivative of the second Piola-Kirchhoff stress w.r.t. the previous ISVs
             */

            constexpr unsigned int dim = 3;

            const unsigned int num_isvs = get_dPK2MeanStressdPreviousISVs( )->size( ); 

            auto dPK2StressdPreviousISVs = get_setDataStorage_dPK2StressdPreviousISVs( );

            *dPK2StressdPreviousISVs.value = *get_dPK2IsochoricStressdPreviousISVs( );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < num_isvs; j++ ){

                    ( *dPK2StressdPreviousISVs.value )[ dim * num_isvs * i + num_isvs* i + j ]
                        += ( *get_dPK2MeanStressdPreviousISVs( ) )[ j ];

                }

            }

        }

        void residual::setdCauchyStressdT( ){
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the temperature
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            auto dCauchyStressdT = get_setDataStorage_dCauchyStressdT( );

            dCauchyStressdT.zero( sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdPK2Stress( get_dCauchyStressdPK2Stress( )->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Vector< floatType, sot_dim > > map_dPK2StressdT( get_dPK2StressdT( )->data( ), sot_dim );

            Eigen::Map< Eigen::Vector< floatType, sot_dim > > map_dCauchyStressdT( dCauchyStressdT.value->data( ), sot_dim );

            map_dCauchyStressdT = ( map_dCauchyStressdPK2Stress * map_dPK2StressdT ).eval( );

        }

        void residual::setdCauchyStressdPreviousT( ){
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the previous temperature
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            auto dCauchyStressdPreviousT = get_setDataStorage_dCauchyStressdPreviousT( );

            dCauchyStressdPreviousT.zero( sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdPK2Stress( get_dCauchyStressdPK2Stress( )->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Vector< floatType, sot_dim > > map_dPK2StressdPreviousT( get_dPK2StressdPreviousT( )->data( ), sot_dim );

            Eigen::Map< Eigen::Vector< floatType, sot_dim > > map_dCauchyStressdPreviousT( dCauchyStressdPreviousT.value->data( ), sot_dim );

            map_dCauchyStressdPreviousT = ( map_dCauchyStressdPK2Stress * map_dPK2StressdPreviousT ).eval( );

        }

        void residual::setdCauchyStressdPreviousISVs( ){
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the previous internal state variables
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_isvs = get_dPK2StressdPreviousISVs( )->size( ) / sot_dim;

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdPK2Stress( get_dCauchyStressdPK2Stress( )->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dPK2StressdPreviousISVs( get_dPK2StressdPreviousISVs( )->data( ), sot_dim, num_isvs );

            auto dCauchyStressdPreviousISVs = get_setDataStorage_dCauchyStressdPreviousISVs( );

            dCauchyStressdPreviousISVs.zero( sot_dim * num_isvs );

            Eigen::Map< Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dCauchyStressdPreviousISVs( dCauchyStressdPreviousISVs.value->data( ), sot_dim, num_isvs );

            map_dCauchyStressdPreviousISVs = ( map_dCauchyStressdPK2Stress * map_dPK2StressdPreviousISVs ).eval( );

        }

        void residual::setPreviousdCauchyStressdT( ){
            /*!
             * Set previous the derivative of the Cauchy stress w.r.t. the temperature
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            auto previousdCauchyStressdT = get_setDataStorage_previousdCauchyStressdT( );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdPK2Stress( get_dCauchyStressdPK2Stress( )->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Vector< floatType, sot_dim > > map_previousdPK2StressdT( get_previousdPK2StressdT( )->data( ), sot_dim );

            previousdCauchyStressdT.zero( sot_dim );

            Eigen::Map< Eigen::Vector< floatType, sot_dim > > map_previousdCauchyStressdT( previousdCauchyStressdT.value->data( ), sot_dim );

            map_previousdCauchyStressdT = ( map_dCauchyStressdPK2Stress * map_previousdPK2StressdT ).eval( );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature.
             */

            auto dRdT = get_setDataStorage_dRdT( );

            *dRdT.value = *get_dCauchyStressdT( );

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
