/**
  ******************************************************************************
  * \file tardigrade_hydraPerzynaViscodamage.cpp
  ******************************************************************************
  * An implementation of perzynaViscodamage using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraPerzynaViscodamage.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace perzynaViscodamage{


        void residual::setDamage( ){
            /*!
             * Get the value of the damage from the evolved state variables
             */

            auto damage = get_setDataStorage_damage( );
            *damage.value = ( *get_plasticStateVariables( ) )[ damageISVIndex ];

        }

        void residual::setDamageJacobians( const bool withPrevious){
            /*!
             * Set the damage along with the jacobians of the damage.
             * 
             * \param withPrevious: Flag for whether to include the derivatives w.r.t. the previous values.
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            if ( withPrevious ){

                auto dDamagedPreviousCauchyStress   = get_setDataStorage_dDamagedPreviousCauchyStress( );

                auto dDamagedPreviousF              = get_setDataStorage_dDamagedPreviousF( );

                auto dDamagedPreviousSubFs          = get_setDataStorage_dDamagedPreviousSubFs( );

                auto dDamagedPreviousT              = get_setDataStorage_dDamagedPreviousT( );

                auto dDamagedPreviousStateVariables = get_setDataStorage_dDamagedPreviousStateVariables( );

                *dDamagedPreviousCauchyStress.value   = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousCauchyStress( ), num_isvs, sot_dim, damageISVIndex );
                
                *dDamagedPreviousF.value              = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousF( ), num_isvs, sot_dim, damageISVIndex );

                *dDamagedPreviousSubFs.value          = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousSubFs( ), num_isvs, ( num_configs - 1 ) * sot_dim, damageISVIndex );

                *dDamagedPreviousT.value              = ( *get_dPlasticStateVariablesdPreviousT( ) )[ damageISVIndex ];

                *dDamagedPreviousStateVariables.value = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousStateVariables( ), num_isvs, num_isvs, damageISVIndex );

            }

            auto damage = get_setDataStorage_damage( );
            *damage.value = ( *get_plasticStateVariables( ) )[ damageISVIndex ];

            auto dDamagedCauchyStress   = get_setDataStorage_dDamagedCauchyStress( );

            auto dDamagedF              = get_setDataStorage_dDamagedF( );

            auto dDamagedSubFs          = get_setDataStorage_dDamagedSubFs( );

            auto dDamagedT              = get_setDataStorage_dDamagedT( );

            auto dDamagedStateVariables = get_setDataStorage_dDamagedStateVariables( );

            *dDamagedCauchyStress.value   = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdCauchyStress( ), num_isvs, sot_dim, damageISVIndex );
            
            *dDamagedF.value              = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdF( ), num_isvs, sot_dim, damageISVIndex );

            *dDamagedSubFs.value          = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdSubFs( ), num_isvs, ( num_configs - 1 ) * sot_dim, damageISVIndex );

            *dDamagedT.value              = ( *get_dPlasticStateVariablesdT( ) )[ damageISVIndex ];

            *dDamagedStateVariables.value = tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdStateVariables( ), num_isvs, num_isvs, damageISVIndex );

        }

        void residual::setDamageJacobians( ){
            /*!
             * Set the damage along with the jacobians of the damage.
             * This routine will not set the jacobians w.r.t. the previous values.
             */

            setDamageJacobians( false );

        }

        void residual::setAllDamageJacobians( ){
            /*!
             * Set the damage along with the jacobians of the damage.
             * This routine will set the jacobians w.r.t. the previous values.
             */

            setDamageJacobians( true );

        }

        void residual::setDamageDeformationGradient( ){
            /*!
             * Set the deformation gradient of the damage
             * 
             * We define the damage through the strain as
             * 
             * \f$ \bf{E}^d = D \bf{E}^r \f$
             * 
             * where \f$\bf{E}^d\f$ is the Green-Lagrange damage strain,
             * \f$D\f$ is the damage and \f$\bf{E}^r\f$ is the recoverable
             * strain and is defined as
             * 
             * \f$ \bf{E}^r = \bf{E}^e + \bf{E}^d \f$ where \f$\bf{E}^e\f$
             * is the elastic strain. We can solve for the damage deformation
             * gradient using the definition of the Green-Lagrange strain and
             * knowledge of the elastic configuration.
             */

            const unsigned int dim = hydra->getDimension( );
            const unsigned int sot_dim = hydra->getSOTDimension( );
            const unsigned int elastic_config_index = *getElasticConfigurationIndex( );

            auto Fd = get_setDataStorage_damageDeformationGradient( );

            // Get the elastic deformation gradient
            floatVector Fe = floatVector( hydra->get_configurations( )->begin( ) + sot_dim * elastic_config_index,
                                          hydra->get_configurations( )->begin( ) + sot_dim * ( elastic_config_index + 1 ) );
    
            // Compute the elastic Green-Lagrange strain
            floatVector Ee;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee ) );
    
            // Compute the damage strain
            floatVector Ed = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * Ee;
    
            // Compute the square root to solve for the damage deformation gradient
            floatVector eye( sot_dim );
            tardigradeVectorTools::eye( eye );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( *Fd.value = tardigradeVectorTools::matrixSqrt( 2.0 * Ed + eye, dim ) );

        }

        void residual::setDamageDeformationGradientJacobians( ){
            /*!
             * Set the jacobians for the damage deformation gradient
             */

            setDamageDeformationGradientJacobians( false );

        }

        void residual::setAllDamageDeformationGradientJacobians( ){
            /*!
             * Set all of the jacobians for the damage deformation gradient
             */

            setDamageDeformationGradientJacobians( true );

        }

        void residual::setDamageDeformationGradientJacobians( const bool withPrevious ){
            /*!
             * Set the deformation gradient of the damage
             * 
             * We define the damage through the strain as
             * 
             * \f$ \bf{E}^d = D \bf{E}^r \f$
             * 
             * where \f$\bf{E}^d\f$ is the Green-Lagrange damage strain,
             * \f$D\f$ is the damage and \f$\bf{E}^r\f$ is the recoverable
             * strain and is defined as
             * 
             * \f$ \bf{E}^r = \bf{E}^e + \bf{E}^d \f$ where \f$\bf{E}^e\f$
             * is the elastic strain. We can solve for the damage deformation
             * gradient using the definition of the Green-Lagrange strain and
             * knowledge of the elastic configuration.
             * 
             * \param withPrevious: Flag for whether to set the Jacobians w.r.t. the previous unknowns
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int elastic_config_index = *getElasticConfigurationIndex( );
    
            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            auto Fd = get_setDataStorage_damageDeformationGradient( );

            // Get the elastic deformation gradient
            floatVector Fe = floatVector( hydra->get_configurations( )->begin( ) + sot_dim * elastic_config_index,
                                          hydra->get_configurations( )->begin( ) + sot_dim * ( elastic_config_index + 1 ) );
    
            // Compute the elastic Green-Lagrange strain
            floatVector Ee;

            floatVector dEedFe;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee, dEedFe ) );

            floatVector dFedF( fot_dim, 0 );

            floatVector dFedSubFs( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            if ( elastic_config_index == 0 ){

                dFedF     = *hydra->get_dF1dF( );

                dFedSubFs = *hydra->get_dF1dFn( );

            }
            else{

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    dFedSubFs[ ( num_configs - 1 ) * sot_dim * i + i + elastic_config_index - 1 ] = 1;

                }

            }

            auto map_dEedFe    = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dEedFe.data( ) );

            auto map_dFedF     = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dFedF.data( ) );
 
            auto map_dFedSubFs = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dFedSubFs.data( ), ( num_configs - 1 ) * sot_dim );
 
            // Compute the damage strain
            floatVector Ed = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * Ee;

            floatVector dEddD = 1 / ( 1 - ( *get_damage( ) ) ) * ( 1 + ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) ) * Ee;

            fourthOrderTensor dEddF( fot_dim, 0 );
            auto map_dEddF = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dEddF.data( ) );

            floatVector dEddSubFs( ( num_configs - 1 ) * fot_dim, 0 );
            auto map_dEddSubFs = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dEddSubFs.data( ), ( num_configs - 1 ) * sot_dim );

            map_dEddF = ( ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * map_dEedFe * map_dFedF ).eval( );

            map_dEddSubFs = ( ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * map_dEedFe * map_dFedSubFs ).eval( );

            // Compute the square root to solve for the damage deformation gradient
            floatVector eye( sot_dim );
            tardigradeVectorTools::eye( eye );
    
            floatMatrix _dAdFe; //A = 2.0 * Ed + eye
                                //
            floatVector dAdFe; //A = 2.0 * Ed + eye

            TARDIGRADE_ERROR_TOOLS_CATCH( *Fd.value = tardigradeVectorTools::matrixSqrt( 2.0 * Ed + eye, dim, _dAdFe ) );

            dAdFe = tardigradeVectorTools::appendVectors( _dAdFe );
            auto map_dAdFe = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dAdFe.data( ) );

            fourthOrderTensor dFddEd( fot_dim, 0 );
            auto map_dFddEd = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dFddEd.data( ) );
            map_dFddEd = ( 2 * map_dAdFe.inverse( ) ).eval( );

            auto map_dEddD = getFixedSizeVectorMap< floatType, sot_dim >( dEddD.data( ) );

            secondOrderTensor dFddD( sot_dim, 0 );
            auto map_dFddD = getFixedSizeVectorMap< floatType, sot_dim >( dFddD.data( ) );

            fourthOrderTensor dFddF( fot_dim, 0 );
            auto map_dFddF = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dFddF.data( ) );

            fourthOrderTensor dFddSubFs( ( num_configs - 1 ) * fot_dim, 0 );
            auto map_dFddSubFs = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dFddSubFs.data( ), ( num_configs - 1 ) * sot_dim );

            map_dFddD = ( map_dFddEd * map_dEddD ).eval( );

            map_dFddF = ( map_dFddEd * map_dEddF ).eval( );

            map_dFddSubFs = ( map_dFddEd * map_dEddSubFs ).eval( );

            if ( withPrevious ){

                auto map_dDamagedPreviousCauchyStress   = getFixedSizeMatrixMap< floatType, 1, sot_dim >( get_dDamagedPreviousCauchyStress( )->data( ) );

                auto map_dDamagedPreviousF              = getFixedSizeMatrixMap< floatType, 1, sot_dim >( get_dDamagedPreviousF( )->data( ) );

                auto map_dDamagedPreviousSubFs          = getDynamicColumnSizeMatrixMap< floatType, 1 >( get_dDamagedPreviousSubFs( )->data( ), ( num_configs - 1 ) * sot_dim );

                auto map_dDamagedPreviousStateVariables = getDynamicColumnSizeMatrixMap< floatType, 1 >( get_dDamagedPreviousStateVariables( )->data( ), num_isvs );

                auto dDamageDeformationGradientdPreviousCauchyStress = get_setDataStorage_dDamageDeformationGradientdPreviousCauchyStress( );
                auto map_dDamageDeformationGradientdPreviousCauchyStress = dDamageDeformationGradientdPreviousCauchyStress.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dDamageDeformationGradientdPreviousF = get_setDataStorage_dDamageDeformationGradientdPreviousF( );
                auto map_dDamageDeformationGradientdPreviousF = dDamageDeformationGradientdPreviousF.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dDamageDeformationGradientdPreviousSubFs = get_setDataStorage_dDamageDeformationGradientdPreviousSubFs( );
                auto map_dDamageDeformationGradientdPreviousSubFs = dDamageDeformationGradientdPreviousSubFs.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

                auto dDamageDeformationGradientdPreviousT = get_setDataStorage_dDamageDeformationGradientdPreviousT( );

                auto dDamageDeformationGradientdPreviousStateVariables = get_setDataStorage_dDamageDeformationGradientdPreviousStateVariables( );
                auto map_dDamageDeformationGradientdPreviousStateVariables = dDamageDeformationGradientdPreviousStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

                map_dDamageDeformationGradientdPreviousCauchyStress = ( map_dFddD * map_dDamagedPreviousCauchyStress ).eval( );

                map_dDamageDeformationGradientdPreviousF            = ( map_dFddD * map_dDamagedPreviousF ).eval( );

                map_dDamageDeformationGradientdPreviousSubFs        = ( map_dFddD * map_dDamagedPreviousSubFs ).eval( );

                *dDamageDeformationGradientdPreviousT.value = dFddD * ( *get_dDamagedPreviousT( ) );

                map_dDamageDeformationGradientdPreviousStateVariables = ( map_dFddD * map_dDamagedPreviousStateVariables ).eval( );

            }

            auto map_dDamagedCauchyStress   = getFixedSizeMatrixMap< floatType, 1, sot_dim >( get_dDamagedCauchyStress( )->data( ) );

            auto map_dDamagedF              = getFixedSizeMatrixMap< floatType, 1, sot_dim >( get_dDamagedF( )->data( ) );

            auto map_dDamagedSubFs          = getDynamicColumnSizeMatrixMap< floatType, 1 >(  get_dDamagedSubFs( )->data( ), ( num_configs - 1 ) * sot_dim );

            auto map_dDamagedStateVariables = getDynamicColumnSizeMatrixMap< floatType, 1 >(  get_dDamagedStateVariables( )->data( ), num_isvs );

            auto dDamageDeformationGradientdCauchyStress = get_setDataStorage_dDamageDeformationGradientdCauchyStress( );
            auto map_dDamageDeformationGradientdCauchyStress = dDamageDeformationGradientdCauchyStress.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dDamageDeformationGradientdF = get_setDataStorage_dDamageDeformationGradientdF( );
            auto map_dDamageDeformationGradientdF = dDamageDeformationGradientdF.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dDamageDeformationGradientdSubFs = get_setDataStorage_dDamageDeformationGradientdSubFs( );
            auto map_dDamageDeformationGradientdSubFs = dDamageDeformationGradientdSubFs.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

            auto dDamageDeformationGradientdT = get_setDataStorage_dDamageDeformationGradientdT( );

            auto dDamageDeformationGradientdStateVariables = get_setDataStorage_dDamageDeformationGradientdStateVariables( );
            auto map_dDamageDeformationGradientdStateVariables = dDamageDeformationGradientdStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

            map_dDamageDeformationGradientdCauchyStress = ( map_dFddD * map_dDamagedCauchyStress ).eval( );

            map_dDamageDeformationGradientdF            = ( map_dFddD * map_dDamagedF ).eval( );
            *dDamageDeformationGradientdF.value        += dFddF;

            map_dDamageDeformationGradientdSubFs        = ( map_dFddD * map_dDamagedSubFs ).eval( );
            *dDamageDeformationGradientdSubFs.value     += dFddSubFs;

            *dDamageDeformationGradientdT.value         = dFddD * ( *get_dDamagedT( ) );

            map_dDamageDeformationGradientdStateVariables = ( map_dFddD * map_dDamagedStateVariables ).eval( );

        }

        void residual::setStateVariableEvolutionRates( const bool isPrevious ){
            /*!
             * Set the state variable evolution rates including the damage
             * 
             * \param &isPrevious: Whether to compute the previous values or not
             */

            tardigradeHydra::perzynaViscoplasticity::residual::setStateVariableEvolutionRates( isPrevious );

            const floatVector *stateVariableEvolutionRates;

            const floatType   *plasticMultiplier;

            setDataStorageBase< floatVector > evolutionRates;

            if ( isPrevious ){

                stateVariableEvolutionRates = get_previousStateVariableEvolutionRates( );

                plasticMultiplier           = get_previousPlasticMultiplier( );

                evolutionRates              = get_setDataStorage_previousStateVariableEvolutionRates( );

            }
            else{

                stateVariableEvolutionRates = get_stateVariableEvolutionRates( );

                plasticMultiplier           = get_plasticMultiplier( );

                evolutionRates              = get_setDataStorage_stateVariableEvolutionRates( );

            }

            *evolutionRates.value = { ( *stateVariableEvolutionRates )[ 0 ], *plasticMultiplier };

        }

        void residual::setStateVariableEvolutionRateDerivatives( const bool isPrevious ){
            /*!
             * Set the jacobians of the state variable evolution rates including the damage
             * 
             * \param &isPrevious: Whether to compute the previous values or not
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            tardigradeHydra::perzynaViscoplasticity::residual::setStateVariableEvolutionRateDerivatives( isPrevious );

            const floatVector *stateVariableEvolutionRates;

            const floatType   *plasticMultiplier;

            const floatVector *dPlasticMultiplierdCauchyStress;

            const floatVector *dPlasticMultiplierdF;

            const floatVector *dPlasticMultiplierdSubFs;

            const floatType   *dPlasticMultiplierdT;

            const floatVector *dPlasticMultiplierdStateVariables;

            setDataStorageBase< floatVector > evolutionRates;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdCauchyStress;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdF;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdSubFs;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdT;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdStateVariables;

            if ( isPrevious ){

                stateVariableEvolutionRates = get_previousStateVariableEvolutionRates( );

                plasticMultiplier           = get_previousPlasticMultiplier( );

                evolutionRates              = get_setDataStorage_previousStateVariableEvolutionRates( );

                dPlasticMultiplierdCauchyStress   = get_dPreviousPlasticMultiplierdPreviousCauchyStress( );

                dPlasticMultiplierdF              = get_dPreviousPlasticMultiplierdPreviousF( );

                dPlasticMultiplierdSubFs          = get_dPreviousPlasticMultiplierdPreviousSubFs( );

                dPlasticMultiplierdT              = get_dPreviousPlasticMultiplierdPreviousT( );

                dPlasticMultiplierdStateVariables = get_dPreviousPlasticMultiplierdPreviousStateVariables( );

                dStateVariableEvolutionRatesdCauchyStress   = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

                dStateVariableEvolutionRatesdF              = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousF( );

                dStateVariableEvolutionRatesdSubFs          = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousSubFs( );

                dStateVariableEvolutionRatesdT              = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousT( );

                dStateVariableEvolutionRatesdStateVariables = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( );

            }
            else{

                stateVariableEvolutionRates = get_stateVariableEvolutionRates( );

                plasticMultiplier           = get_plasticMultiplier( );

                evolutionRates              = get_setDataStorage_stateVariableEvolutionRates( );

                dPlasticMultiplierdCauchyStress   = get_dPlasticMultiplierdCauchyStress( );

                dPlasticMultiplierdF              = get_dPlasticMultiplierdF( );

                dPlasticMultiplierdSubFs          = get_dPlasticMultiplierdSubFs( );

                dPlasticMultiplierdT              = get_dPlasticMultiplierdT( );

                dPlasticMultiplierdStateVariables = get_dPlasticMultiplierdStateVariables( );

                dStateVariableEvolutionRatesdCauchyStress   = get_setDataStorage_dStateVariableEvolutionRatesdCauchyStress( );

                dStateVariableEvolutionRatesdF              = get_setDataStorage_dStateVariableEvolutionRatesdF( );

                dStateVariableEvolutionRatesdSubFs          = get_setDataStorage_dStateVariableEvolutionRatesdSubFs( );

                dStateVariableEvolutionRatesdT              = get_setDataStorage_dStateVariableEvolutionRatesdT( );

                dStateVariableEvolutionRatesdStateVariables = get_setDataStorage_dStateVariableEvolutionRatesdStateVariables( );

            }

            const unsigned int num_isvs = stateVariableEvolutionRates->size( );

            *evolutionRates.value = { ( *stateVariableEvolutionRates )[ 0 ], *plasticMultiplier };

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *dStateVariableEvolutionRatesdCauchyStress.value )[ sot_dim + i ] = ( *dPlasticMultiplierdCauchyStress )[ i ];

            }

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *dStateVariableEvolutionRatesdF.value )[ sot_dim + i ] = ( *dPlasticMultiplierdF )[ i ];

            }

            for ( unsigned int i = 0; i < ( num_configs - 1 ) * sot_dim; i++ ){

                ( *dStateVariableEvolutionRatesdSubFs.value )[ ( num_configs - 1 ) * sot_dim + i ] = ( *dPlasticMultiplierdSubFs )[ i ];

            }

            dStateVariableEvolutionRatesdT.value->resize( dStateVariableEvolutionRatesdT.value->size( ) );
            for ( unsigned int i = 0; i < 1; i++ ){

                ( *dStateVariableEvolutionRatesdT.value )[ i + num_isvs ] = *dPlasticMultiplierdT;

            }

            dStateVariableEvolutionRatesdStateVariables.value->resize( 2 * num_isvs + 2 );
            ( *dStateVariableEvolutionRatesdStateVariables.value )[ 1 ] = 0;
            ( *dStateVariableEvolutionRatesdStateVariables.value )[ 2 ] = ( *dPlasticMultiplierdStateVariables )[ 0 ];
            ( *dStateVariableEvolutionRatesdStateVariables.value )[ 3 ] = 0;

        }

        void residual::setResidual( ){
            /*!
             * Set the residual vector
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            auto residual = get_setDataStorage_residual( );
            residual.zero( *getNumEquations( ) );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *residual.value )[ i ] = hydra->getConfiguration( *getDamageConfigurationIndex( ) )[ i ] - ( *get_damageDeformationGradient( ) )[ i ];

            }

            for ( unsigned int i = 0; i < num_isvs; i++ ){

                ( *residual.value )[ sot_dim + i ] = ( *get_stateVariables( ) )[ i ] - ( *get_plasticStateVariables( ) )[ i ];

            }

        }

        void residual::setJacobian( ){
            /*!
             * Set the Jacobian matrix
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int damage_configuration_index = *getDamageConfigurationIndex( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            auto jacobian = get_setDataStorage_jacobian( );
            jacobian.zero( *getNumEquations( ) * num_unknowns );

            // Jacobians of the damage deformation gradient
            for ( unsigned int i = 0; i < sot_dim; i++ ){
                unsigned int row = i;

                // Jacobians w.r.t. the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dDamageDeformationGradientdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Jacobians w.r.t. the sub configurations
                ( *jacobian.value )[ num_unknowns * row + sot_dim * damage_configuration_index + i ] += 1;

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){
                    unsigned int col = sot_dim + j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dDamageDeformationGradientdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Jacobians w.r.t. the state variables
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dDamageDeformationGradientdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            // Jacobians of the damage hardening state variables
            for ( unsigned int i = 0; i < num_isvs; i++ ){
                unsigned int row = sot_dim + i;

                // Jacobians w.r.t. the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dPlasticStateVariablesdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Jacobians w.r.t. the sub configurations
                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){
                    unsigned int col = sot_dim + j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dPlasticStateVariablesdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Jacobians w.r.t. the state variables
                ( *jacobian.value )[ num_unknowns * row + sot_dim + ( num_configs - 1 ) * sot_dim + ( *getStateVariableIndices( ) )[ i ] ] += 1;
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dPlasticStateVariablesdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_isvs = get_dPlasticStateVariablesdT( )->size( );

            auto dRdT = get_setDataStorage_dRdT( );

            dRdT.zero( *getNumEquations( ) );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *dRdT.value )[ i ] -= ( *get_dDamageDeformationGradientdT( ) )[ i ];

            }

            for ( unsigned int i = 0; i < num_isvs; i++ ){

                ( *dRdT.value )[ sot_dim + i ] = -( *get_dPlasticStateVariablesdT( ) )[ i ];

            }

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            auto dRdF = get_setDataStorage_dRdF( );
            dRdF.zero( ( *getNumEquations( ) ) * sot_dim );

            *dRdF.value = tardigradeVectorTools::appendVectors( { *get_dDamageDeformationGradientdF( ),
                                                                  *get_dPlasticStateVariablesdF( ) } );

            std::transform( dRdF.value->cbegin( ), dRdF.value->cend( ), dRdF.value->begin( ), std::negate<floatType>( ) );

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the incoming parameter vector
             * 
             * \param &parameters: The incoming parameter vector
             * 
             * The form has the same interpretation as the base viscoplastic case
             * except we have a different method of hardening for the damage state
             * variable which we account for here.
             */

            tardigradeHydra::perzynaViscoplasticity::residual::decomposeParameters( parameters );

            //Setting the contribution of damage to the calculation of the drag stress to zero
            set_dragStressParameters( { parameters[ 1 ], parameters[  2 ], 0. } );

            //Setting the contribution of damage to the hardening of the damage state variable to zero
            set_hardeningParameters(  { parameters[ 9 ], parameters[ 10 ], 0. } );

        }

        void residual::addParameterizationInfo( std::string &parameterization_info ){
            /*!
             * Add the parameterization info to the incoming string
             * 
             * \param &parameterization_info: The incoming string
             */

            std::stringstream ss;

            parameterization_info += "class: tardigradeHydra::perzynaViscodamage::residual\n\n";
            parameterization_info += "Getting information from parent class\n\n";
            tardigradeHydra::perzynaViscoplasticity::residual::addParameterizationInfo( parameterization_info );
            ss.precision(9);
            ss << std::scientific;
            ss << "\n\n";
            ss << "Modifying the drag-stress and hardening parameters\n\n";
            ss << "name,                              description,  units, current value\n";
            ss << "  q0,                  the initial drag stress, stress, " << ( *get_dragStressParameters( ) )[ 0 ] << "\n";
            ss << "  Ei,    the drag stress isv hardening modulus, stress, " << ( *get_dragStressParameters( ) )[ 1 ] << "\n";
            ss << "  ED, the drag stress damage hardening modulus, stress, " << ( *get_dragStressParameters( ) )[ 2 ] << "\n";
            ss << " hi0,               isv initial evolution rate,   none, " <<  ( *get_hardeningParameters( ) )[ 0 ] << "\n";
            ss << " hi1,            isv linear isv evolution rate,   none, " <<  ( *get_hardeningParameters( ) )[ 1 ] << "\n"; 
            ss << " hd1,         isv linear damage evolution rate,   none, " <<  ( *get_hardeningParameters( ) )[ 2 ] << "\n";

            ss.unsetf(std::ios_base::floatfield);
            parameterization_info.append(ss.str());

        }

    }

}
