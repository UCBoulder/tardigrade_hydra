/**
 ******************************************************************************
 * \file tardigrade_HyperelasticBase.tpp
 ******************************************************************************
 * Template definitions for HyperelasticBase
 ******************************************************************************
 */

namespace tardigradeHydra {

    /*!
     * Compute the first invariant of the deformation gradient
     *
     * \param isPrevious: Whether to compute this for the current (false) or previous (true) configuration
     */
    template<typename T>
    T HyperelasticBase::compute_I1(const bool isPrevious){

        constexpr unsigned int dim = 3; //TODO: Replace with constant value from container

        const secondOrderTensor *Fe;

        if (isPrevious){
            Fe = get_previousFe();
        }else{
            Fe = get_Fe();
        }

        T I1 = T();

        for ( unsigned int i = 0; i < dim; ++i ){

            for ( unsigned int j = 0; j < dim; ++j ){

                I1 += (*Fe)[dim * i + j] * (*Fe)[dim * i + j];

            }

        }

        return I1;

    }

    /*!
     * Compute the derivative of the first invariant of the deformation gradient
     *
     * \param isPrevious: Whether to compute this for the current (false) or previous (true) configuration
     * \param dI1dFe_begin: The starting iterator of the container of the derivative of I1 with respect to the deformation gradient
     * \param dI1dFe_end: The stopping iterator of the container of the derivative of I1 with respect to the deformation gradient
     */
    template<class dI1dFe_iter>
    void HyperelasticBase::compute_dI1dFe(const bool isPrevious, dI1dFe_iter dI1dFe_begin, dI1dFe_iter dI1dFe_end){

        using dI1dFe_type = typename std::iterator_traits<dI1dFe_iter>::value_type;

        constexpr unsigned int dim = 3; //TODO: Replace with constant value from container

        const secondOrderTensor *Fe;

        if (isPrevious){
            Fe = get_previousFe();
        }else{
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dI1dFe_end - dI1dFe_begin) == dim * dim, "dI1dFe has a size of " + std::to_string((unsigned int)(dI1dFe_end - dI1dFe_begin)) + " but it should have a size of " + std::to_string(dim*dim))

        std::fill(dI1dFe_begin, dI1dFe_end, dI1dFe_type());

        for ( unsigned int i = 0; i < dim; ++i ){

            for ( unsigned int j = 0; j < dim; ++j ){

                *(dI1dFe_begin + dim * i + j) = 2 * (*Fe)[dim * i + j];

            }

        }

    }

    /*!
     * Compute the second derivative of the first invariant of the deformation gradient
     *
     * \param d2I1dFe2_begin: The starting iterator of the container of the second derivative of I1 with respect to the deformation gradient
     * \param d2I1dFe2_end: The stopping iterator of the container of the second derivative of I1 with respect to the deformation gradient
     */
    template<class d2I1dFe2_iter>
    void HyperelasticBase::compute_d2I1dFe2(d2I1dFe2_iter d2I1dFe2_begin, d2I1dFe2_iter d2I1dFe2_end){

        using d2I1dFe2_type = typename std::iterator_traits<d2I1dFe2_iter>::value_type;

        constexpr unsigned int dim = 3; //TODO: Replace with constant value from container

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(d2I1dFe2_end - d2I1dFe2_begin) == dim * dim * dim * dim, "d2I1dFe2 has a size of " + std::to_string((unsigned int)(d2I1dFe2_end - d2I1dFe2_begin)) + " but it should have a size of " + std::to_string(dim*dim*dim*dim))

        std::fill(d2I1dFe2_begin, d2I1dFe2_end, d2I1dFe2_type());

        for ( unsigned int i = 0; i < dim; ++i ){

            for ( unsigned int j = 0; j < dim; ++j ){

                *(d2I1dFe2_begin + dim * dim * dim * i + dim * dim * j + dim * i + j) += 2;

            }

        }

    }

}
