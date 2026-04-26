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
     * Compute the first invariant of the deformation
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
     * Compute the derivative of the first invariant of the deformation
     *
     * \param isPrevious: Whether to compute this for the current (false) or previous (true) configuration
     * \param dI1dFe_begin: The starting iterator of the container of the derivative of I1 with respect to the deformation
     * \param dI1dFe_end: The stopping iterator of the container of the derivative of I1 with respect to the deformation
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
     * \param d2I1dFe2_begin: The starting iterator of the container of the second derivative of I1 with respect to the deformation
     * \param d2I1dFe2_end: The stopping iterator of the container of the second derivative of I1 with respect to the deformation
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

    /*!
     * Compute the Right Cauchy-Green deformation tensor
     *
     * \f$ C_{IJ} = F_{iI} F_{iJ} \f$
     *
     * \param isPrevious: Whether to compute the current or previous value
     * \param C_begin: The starting iterator of the right Cauchy-Green deformation tensor
     * \param C_end: The stopping iterator of the right Cauchy-Green deformation tensor
     */
    template<typename C_iter>
    void HyperelasticBase::compute_right_cauchy_green_deformation_tensor(const bool isPrevious, C_iter C_begin, C_iter C_end){

        constexpr unsigned int dim = 3; //TODO: Get this from ResidualBase

        const secondOrderTensor *Fe;

        if (isPrevious){
            Fe = get_previousFe();
        }else{
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(std::end(*Fe) - std::begin(*Fe)) == dim * dim, "The storage array for F is of size " + std::to_string((unsigned int)(std::end(*Fe) - std::begin(*Fe))) + " but must be of size " + std::to_string(dim*dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(C_end - C_begin) == dim * dim, "The storage array for C is of size " + std::to_string((unsigned int)(C_end - C_begin)) + " but must be of size " + std::to_string(dim*dim))

        for ( unsigned int i = 0; i < dim; ++i ){
            for ( unsigned int I = 0; I < dim; ++I ){
                for ( unsigned int J = 0; J < dim; ++J ){
                    *(C_begin + dim * I + J) += (*Fe)[dim * i + I] * (*Fe)[dim * i + J];
                }
            }
        }
    }

    /*!
     * Compute the second invariant of the deformation
     *
     * \param isPrevious: Whether to compute this for the current (false) or previous (true) configuration
     */
    template<typename T>
    T HyperelasticBase::compute_I2(const bool isPrevious){

        constexpr unsigned int dim = 3; //TODO: Get this from ResidualBase

        T I1 = compute_I1<T>(isPrevious);

        T I2 = 0.5 * I1 * I1;

        std::array<T, dim*dim> C{};
        compute_right_cauchy_green_deformation_tensor(isPrevious, std::begin(C), std::end(C));

        I2 -= 0.5 * std::inner_product(std::begin(C), std::end(C), std::begin(C), T());

        return I2;

    }

    /*!
     * Compute the derivative of the second invariant of the deformation with respect to the elastic deformation gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param dI2dFe_begin: The starting iterator of the derivative of the second inveriant of the deformation with respect ot the elastic deformation gradient
     * \param dI2dFe_end: The stopping iterator of the derivative of the second inveriant of the deformation with respect ot the elastic deformation gradient
     */
    template<class dI2dFe_iter>
    void HyperelasticBase::compute_dI2dFe(const bool isPrevious, dI2dFe_iter dI2dFe_begin, dI2dFe_iter dI2dFe_end){

        using dI2dFe_type = typename std::iterator_traits<dI2dFe_iter>::value_type;

        constexpr unsigned int dim = 3; //TODO: Get this value from ResidualBase

        const secondOrderTensor *Fe;

        if (isPrevious){
            Fe = get_previousFe();
        }else{
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dI2dFe_end - dI2dFe_begin) == dim * dim, "The size of dI2dFe is " + std::to_string((unsigned int)(dI2dFe_end - dI2dFe_begin)) + " but it must be " + std::to_string(dim*dim));

        dI2dFe_type I1 = compute_I1<dI2dFe_type>(isPrevious);

        std::transform(std::begin(*Fe), std::end(*Fe), dI2dFe_begin,
                std::bind(std::multiplies<>(), std::placeholders::_1, I1));

        std::array<dI2dFe_type,dim*dim> C{};
        compute_right_cauchy_green_deformation_tensor(isPrevious, std::begin(C), std::end(C));

        for ( unsigned int a = 0; a < dim; ++a){
            for ( unsigned int I = 0; I < dim; ++I){
                for ( unsigned int A = 0; A < dim; ++A){
                    *(dI2dFe_begin + dim * a + A) -= (*Fe)[dim * a + I] * C[dim * I + A];
                }
            }
        }
        
        std::transform(dI2dFe_begin, dI2dFe_end, dI2dFe_begin,
                std::bind(std::multiplies<>(), std::placeholders::_1, 2));

    }

    /*!
     * Compute the second derivative of the second invariant of the deformation with respect to the elastic deformation gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param d2I2dFe2_begin: The starting iterator of the second derivative of the second inveriant of the deformation with respect ot the elastic deformation gradient
     * \param d2I2dFe2_end: The stopping iterator of the second derivative of the second inveriant of the deformation with respect ot the elastic deformation gradient
     */
    template<class d2I2dFe2_iter>
    void HyperelasticBase::compute_d2I2dFe2(const bool isPrevious, d2I2dFe2_iter d2I2dFe2_begin, d2I2dFe2_iter d2I2dFe2_end){

        using d2I2dFe2_type = typename std::iterator_traits<d2I2dFe2_iter>::value_type;

        constexpr unsigned int dim = 3; //TODO: Get this value from ResidualBase

        const secondOrderTensor *Fe;

        if (isPrevious){
            Fe = get_previousFe();
        }else{
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(d2I2dFe2_end - d2I2dFe2_begin) == dim * dim * dim * dim, "The size of d2I2dFe2 is " + std::to_string((unsigned int)(d2I2dFe2_end - d2I2dFe2_begin)) + " but it must be " + std::to_string(dim*dim*dim*dim));

        d2I2dFe2_type I1 = compute_I1<d2I2dFe2_type>(isPrevious);

        std::array<d2I2dFe2_type,dim*dim> C{};
        compute_right_cauchy_green_deformation_tensor(isPrevious, std::begin(C), std::end(C));
        std::fill(d2I2dFe2_begin, d2I2dFe2_end, d2I2dFe2_type());

        for (unsigned int aA = 0; aA < dim * dim; ++aA){
            *(d2I2dFe2_begin + dim * dim * aA + aA) += I1;
        }

        for (unsigned int a = 0; a < dim; ++a){
            for ( unsigned int A = 0; A < dim; ++A){
                for ( unsigned int b = 0; b < dim; ++b){
                    *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * a + b) -= C[dim*A+b];
                    for ( unsigned int B = 0; B < dim; ++B){
                        *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * b + B) += 2 * (*Fe)[dim*a+A] * (*Fe)[dim*b+B];
                        *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * b + B) -= (*Fe)[dim*a+B] * (*Fe)[dim*b+A];
                        *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * b + A) -= (*Fe)[dim*a+B] * (*Fe)[dim*b+B];
                    }
                }
            }
        }

        std::transform(d2I2dFe2_begin, d2I2dFe2_end, d2I2dFe2_begin,
                std::bind(std::multiplies<>(), std::placeholders::_1, 2));
    }

}
