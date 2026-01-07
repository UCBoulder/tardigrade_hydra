/**
******************************************************************************
* \file tardigrade_CustomErrors.h
******************************************************************************
* Custom errors for the hydra library
******************************************************************************
*/

#ifndef TARDIGRADE_CUSTOMERRORS_H
#define TARDIGRADE_CUSTOMERRORS_H

#include<exception>
#include<string>

namespace tardigradeHydra{

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

}

#include"tardigrade_CustomErrors.cpp"

#endif
