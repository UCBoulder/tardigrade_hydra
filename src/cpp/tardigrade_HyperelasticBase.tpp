/**
 ******************************************************************************
 * \file tardigrade_HyperelasticBase.tpp
 ******************************************************************************
 * Template definitions for HyperelasticBase
 ******************************************************************************
 */

#include<numeric>
#include<functional>

namespace tardigradeHydra {

    /*!
     * Compute the first invariant of the deformation gradient
     *
     * \param isPrevious: Whether to compute this for the current (false) or previous (true) configuration
     */
    template<typename T>
    T HyperelasticBase::compute_I1(const bool isPrevious){

        const secondOrderTensor *Fe;

        if (isPrevious){
            Fe = get_previousFe();
        }else{
            Fe = get_Fe();
        }

        T I1 = std::inner_product(std::begin(*Fe),std::end(*Fe),std::begin(*Fe),T());

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

        TARDIGRADE_ERROR_TOOLS_EVAL(
        constexpr unsigned int dim = 3; //TODO: Replace with constant value from container
        )

        const secondOrderTensor *Fe;

        if (isPrevious){
            Fe = get_previousFe();
        }else{
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dI1dFe_end - dI1dFe_begin) == dim * dim, "dI1dFe has a size of " + std::to_string((unsigned int)(dI1dFe_end - dI1dFe_begin)) + " but it should have a size of " + std::to_string(dim*dim))


        std::transform(std::begin(*Fe), std::end(*Fe), dI1dFe_begin,
               std::bind(std::multiplies<>(), std::placeholders::_1, 2));

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

        for ( unsigned int ij = 0; ij < dim * dim; ++ij ){

                *(d2I1dFe2_begin + dim * dim * ij + ij) += 2;

        }

    }

}
