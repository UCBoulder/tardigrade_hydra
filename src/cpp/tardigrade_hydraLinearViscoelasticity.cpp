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
   
            const unsigned int *dim = hydra->getDimension( );
 
            if ( parameters.size( ) < 10 ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Parameter vector is expected to have a length of at least 10 but has a length of " + std::to_string( parameters.size( ) ) ) );
    
            }

            setNumVolumetricViscousTerms( ( unsigned int )( parameters[ 0 ] + 0.5 ) );

            setNumIsochoricViscousTerms( ( unsigned int )( parameters[ 1 ] + 0.5 ) );

            setNumStateVariables( *getNumVolumetricViscousTerms( ) + ( *dim ) * ( *dim ) * ( *getNumIsochoricViscousTerms( ) ) );

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

            floatVector Fehat;

            TARDIGRADE_ERROR_TOOLS_CATCH( decomposeDeformation( *get_Fe( ), Je, Fehat ) );

            set_Je( Je );

            set_Fehat( Fehat );

        }

        void residual::decomposePreviousElasticDeformation( ){
            /*!
             * Decompose the previous elastic deformation into volumetric and isochoric parts
             */

            floatType previousJe;

            floatVector previousFehat;

            TARDIGRADE_ERROR_TOOLS_CATCH( decomposeDeformation( *get_previousFe( ), previousJe, previousFehat ) );

            set_previousJe( previousJe );

            set_previousFehat( previousFehat );

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

            TARDIGRADE_ERROR_TOOLS_CATCH( J = tardigradeVectorTools::determinant( F, ( *dim ), ( *dim ) ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( Fhat = F / std::pow( J, 1./3 ) );

        }

        void residual::setdJedFe( const bool isPrevious ){
            /*!
             * Set the derivative of the elastic Jacobian of deformation w.r.t.
             * the elastic deformation gradient.
             * 
             * \param isPrevious: Flag for if the derivative is of the current (false) or previous (true) value
             */

            const unsigned int* dim = hydra->getDimension( );

            const floatVector *Fe;

            if ( isPrevious ){

                Fe = get_previousFe( );

            }
            else{

                Fe = get_Fe( );

            }

            floatVector dJedFe = tardigradeVectorTools::computeDDetADA( *Fe, ( *dim ), ( *dim ) );

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

            const unsigned int* dim = hydra->getDimension( );

            const floatType   *Je;

            const floatVector *Fe;

            const floatVector *dJedFe;

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

            floatMatrix dFehatdFe = tardigradeVectorTools::eye< floatType >( ( *dim ) * ( *dim ) ) * std::pow( ( *Je ), -1. / 3 );

            dFehatdFe -= tardigradeVectorTools::dyadic( *Fe, *dJedFe ) * std::pow( ( *Je ), -4. / 3 ) / 3.;

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

            const unsigned int *dim = hydra->getDimension( );

            unsigned int lb = *getViscoelasticISVLowerIndex( );
            unsigned int ub = lb + *getNumVolumetricViscousTerms( );

            volumetricStateVariables = floatVector( hydra->get_additionalStateVariables( )->begin( ) + lb,
                                                    hydra->get_additionalStateVariables( )->begin( ) + ub );

            lb = ub;

            ub = lb + ( *dim ) * ( *dim ) * *getNumIsochoricViscousTerms( );

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

            const floatType *Je;

            const floatVector *dJedFe;

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

            floatMatrix dPK2MeanStressdJe;

            floatVector dPK2MeanStressdRateModifier;

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
                                                                                                     dPK2MeanStressdJe, dPK2MeanStressdRateModifier ) );

            if( isPrevious ){

                set_previousPK2MeanStress( PK2MeanStress[ 0 ] );
    
                set_previousdPK2MeanStressdFe( dPK2MeanStressdJe[ 0 ][ 0 ] * ( *dJedFe ) );
    
                set_previousdPK2MeanStressdT( dPK2MeanStressdRateModifier[ 0 ] * ( *dVolumetricRateMultiplierdT ) );

            }
            else{

                set_PK2MeanStress( PK2MeanStress[ 0 ] );
    
                set_volumetricViscoelasticStateVariables( currentVolumetricStateVariables );
    
                set_dPK2MeanStressdFe( dPK2MeanStressdJe[ 0 ][ 0 ] * ( *dJedFe ) );
    
                set_dPK2MeanStressdT( dPK2MeanStressdRateModifier[ 0 ] * ( *dVolumetricRateMultiplierdT ) );

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

            const floatVector *Fehat;

            const floatVector *previousFehat = get_previousFehat( );

            const floatMatrix *dFehatdFe;

            floatType time;

            floatType previousTime = *hydra->getTime( ) - *hydra->getDeltaTime( );

            const floatType *isochoricRateMultiplier;

            const floatType *dIsochoricRateMultiplierdT;

            const floatType *previousIsochoricRateMultiplier = get_previousIsochoricRateMultiplier( );

            if ( isPrevious ){

                time                       = previousTime;

                Fehat                      = get_previousFehat( );

                dFehatdFe                  = get_previousdFehatdFe( );

                isochoricRateMultiplier    = get_previousIsochoricRateMultiplier( );

                dIsochoricRateMultiplierdT = get_dPreviousIsochoricRateMultiplierdPreviousT( );

            }
            else{

                time                       = *hydra->getTime( );

                Fehat                      = get_Fehat( );

                dFehatdFe                  = get_dFehatdFe( );

                isochoricRateMultiplier    = get_isochoricRateMultiplier( );

                dIsochoricRateMultiplierdT = get_dIsochoricRateMultiplierdT( );

            }

            // Compute the strain measures

            floatVector isochoricStrain, previousIsochoricStrain;
            floatMatrix dEehatdFehat;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( *Fehat, isochoricStrain, dEehatdFehat ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( *previousFehat, previousIsochoricStrain ) );

            floatMatrix dEehatdFe = tardigradeVectorTools::dot( dEehatdFehat, *dFehatdFe );

            // Get the previous state variable values

            floatVector previousVolumetricStateVariables;

            floatVector previousIsochoricStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( decomposeStateVariableVector( previousVolumetricStateVariables,
                                                             previousIsochoricStateVariables ) );

            floatVector PK2IsochoricStress;

            floatVector deltaPK2IsochoricStress;

            floatVector currentIsochoricStateVariables;

            floatMatrix dPK2IsochoricStressdEe;

            floatVector dPK2IsochoricStressdRateMultiplier;

            // Compute the viscous isochoric stress
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::linearViscoelasticity( time, isochoricStrain,
                                                                                                     previousTime, previousIsochoricStrain,
                                                                                                     *isochoricRateMultiplier,
                                                                                                     *previousIsochoricRateMultiplier,
                                                                                                     previousIsochoricStateVariables,
                                                                                                     getIsochoricViscoelasticParameters( ),
                                                                                                     *getIntegrationAlpha( ), deltaPK2IsochoricStress,
                                                                                                     PK2IsochoricStress, currentIsochoricStateVariables,
                                                                                                     dPK2IsochoricStressdEe, dPK2IsochoricStressdRateMultiplier ) );

            if ( isPrevious ){

                set_previousPK2IsochoricStress( PK2IsochoricStress );
    
                set_previousdPK2IsochoricStressdFe( tardigradeVectorTools::dot( dPK2IsochoricStressdEe, dEehatdFe ) );
    
                set_previousdPK2IsochoricStressdT( dPK2IsochoricStressdRateMultiplier * ( *dIsochoricRateMultiplierdT ) );

            }
            else{

                set_PK2IsochoricStress( PK2IsochoricStress );
    
                set_isochoricViscoelasticStateVariables( currentIsochoricStateVariables );
    
                set_dPK2IsochoricStressdFe( tardigradeVectorTools::dot( dPK2IsochoricStressdEe, dEehatdFe ) );
    
                set_dPK2IsochoricStressdT( dPK2IsochoricStressdRateMultiplier * ( *dIsochoricRateMultiplierdT ) );

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

            floatVector viscoelasticStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( viscoelasticStateVariables = tardigradeVectorTools::appendVectors( { *get_volumetricViscoelasticStateVariables( ),
                                                                                          *get_isochoricViscoelasticStateVariables( ) } ) );

            tardigradeHydra::residualBase::setCurrentAdditionalStateVariables( viscoelasticStateVariables );

        }

        void residual::setPK2Stress( ){
            /*!
             * Set the PK2 stress
             */

            const unsigned int* dim = hydra->getDimension( );

            floatVector eye( ( *dim ) * ( *dim ), 0 );
            tardigradeVectorTools::eye( eye );

            floatVector PK2Stress = ( *get_PK2IsochoricStress( ) ) + ( *get_PK2MeanStress( ) ) * eye;

            set_PK2Stress( PK2Stress );

        }

        void residual::setdPK2StressdFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic deformation gradient
             */

            const unsigned int* dim = hydra->getDimension( );

            floatVector eye( ( *dim ) * ( *dim ), 0 );
            tardigradeVectorTools::eye( eye );

            floatMatrix dPK2StressdFe = *get_dPK2IsochoricStressdFe( ) + tardigradeVectorTools::dyadic( eye, *get_dPK2MeanStressdFe( ) );

            set_dPK2StressdFe( dPK2StressdFe );

        }

        void residual::setdPK2StressdT( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the temperature
             */

            const unsigned int* dim = hydra->getDimension( );

            floatVector eye( ( *dim ) * ( *dim ), 0 );
            tardigradeVectorTools::eye( eye );

            floatVector dPK2StressdT = *get_dPK2IsochoricStressdT( ) + *get_dPK2MeanStressdT( ) * eye;

            set_dPK2StressdT( dPK2StressdT );

        }

        void residual::setdCauchyStressdT( ){
            /*!
             * Set the derivative of the Cauchy stress w.r.t. the temperature
             */

            floatVector dCauchyStressdT = tardigradeVectorTools::dot( *get_dCauchyStressdPK2Stress( ), *get_dPK2StressdT( ) );

            set_dCauchyStressdT( dCauchyStressdT );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature.
             */

            setdRdT( *get_dCauchyStressdT( ) );

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
