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

namespace tardigradeHydra{

    namespace linearViscoelasticity{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length 2 which defines lambda and mu.
             */
    
            if ( parameters.size( ) < 4 ){
    
                ERROR_TOOLS_CATCH( throw std::runtime_error( "Parameter vector is expected to have a length of at least 4 but has a length of " + std::to_string( parameters.size( ) ) ) );
    
            }

            setNumVolumetricViscousTerms( ( unsigned int )( parameters[ 0 ] + 0.5 ) );

            setNumIsochoricViscousTerms( ( unsigned int )( parameters[ 1 ] + 0.5 ) );

            setKinf( parameters[ 2 ] );

            setGinf( parameters[ 3 ] );

            unsigned int parameterCount = 4 + 2 * ( *getNumVolumetricViscousTerms( ) )
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

            unsigned int lb = 4;
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

            const unsigned int *dim = hydra->getDimension( );

            floatType Je;
            floatVector Fehat;

            ERROR_TOOLS_CATCH( Je = vectorTools::determinant( ( *hydra->getConfigurations( ) )[ 0 ], ( *dim ), ( *dim ) ) );

            ERROR_TOOLS_CATCH( Fehat = ( *hydra->getConfigurations( ) )[ 0 ] / std::pow( Je, 1./3 ) );
            setJe( Je );

            setFehat( Fehat );

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

        void residual::setPK2Stress( ){
            /*!
             * Set the second Piola-Kirchhoff stress
             */

        }

    }

}
