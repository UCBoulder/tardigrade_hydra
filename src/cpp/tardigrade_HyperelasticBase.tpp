/**
 ******************************************************************************
 * \file tardigrade_HyperelasticBase.tpp
 ******************************************************************************
 * Template definitions for HyperelasticBase
 ******************************************************************************
 */

#include <functional>
#include <numeric>

namespace tardigradeHydra {

    /*!
     * Compute the first invariant of the deformation
     *
     * \param isPrevious: Whether to compute this for the current (false) or previous (true) configuration
     */
    template <typename T>
    T HyperelasticBase::compute_I1(const bool isPrevious) {
        const secondOrderTensor *Fe;

        if (isPrevious) {
            Fe = get_previousFe();
        } else {
            Fe = get_Fe();
        }

        T I1 = std::inner_product(std::begin(*Fe), std::end(*Fe), std::begin(*Fe), T());

        return I1;
    }

    /*!
     * Compute the derivative of the first invariant of the deformation
     *
     * \param isPrevious: Whether to compute this for the current (false) or previous (true) configuration
     * \param dI1dFe_begin: The starting iterator of the container of the derivative of I1 with respect to the
     * deformation
     * \param dI1dFe_end: The stopping iterator of the container of the derivative of I1 with respect to the deformation
     */
    template <class dI1dFe_iter>
    void HyperelasticBase::compute_dI1dFe(const bool isPrevious, dI1dFe_iter dI1dFe_begin, dI1dFe_iter dI1dFe_end) {
        TARDIGRADE_ERROR_TOOLS_EVAL(constexpr unsigned int dim = 3;  // TODO: Replace with constant value from container
        )

        const secondOrderTensor *Fe;

        if (isPrevious) {
            Fe = get_previousFe();
        } else {
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dI1dFe_end - dI1dFe_begin) == dim * dim,
                                     "dI1dFe has a size of " +
                                         std::to_string((unsigned int)(dI1dFe_end - dI1dFe_begin)) +
                                         " but it should have a size of " + std::to_string(dim * dim))

        std::transform(std::begin(*Fe), std::end(*Fe), dI1dFe_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, 2));
    }

    /*!
     * Compute the second derivative of the first invariant of the deformation gradient
     *
     * \param d2I1dFe2_begin: The starting iterator of the container of the second derivative of I1 with respect to the
     * deformation
     * \param d2I1dFe2_end: The stopping iterator of the container of the second derivative of I1 with respect to the
     * deformation
     */
    template <class d2I1dFe2_iter>
    void HyperelasticBase::compute_d2I1dFe2(d2I1dFe2_iter d2I1dFe2_begin, d2I1dFe2_iter d2I1dFe2_end) {
        using d2I1dFe2_type = typename std::iterator_traits<d2I1dFe2_iter>::value_type;

        constexpr unsigned int dim = 3;  // TODO: Replace with constant value from container

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(d2I1dFe2_end - d2I1dFe2_begin) == dim * dim * dim * dim,
                                     "d2I1dFe2 has a size of " +
                                         std::to_string((unsigned int)(d2I1dFe2_end - d2I1dFe2_begin)) +
                                         " but it should have a size of " + std::to_string(dim * dim * dim * dim))

        std::fill(d2I1dFe2_begin, d2I1dFe2_end, d2I1dFe2_type());

        for (unsigned int ij = 0; ij < dim * dim; ++ij) {
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
    template <typename C_iter>
    void HyperelasticBase::compute_right_cauchy_green_deformation_tensor(const bool isPrevious, C_iter C_begin,
                                                                         C_iter C_end) {
        constexpr unsigned int dim = 3;  // TODO: Get this from ResidualBase

        const secondOrderTensor *Fe;

        if (isPrevious) {
            Fe = get_previousFe();
        } else {
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(std::end(*Fe) - std::begin(*Fe)) == dim * dim,
                                     "The storage array for F is of size " +
                                         std::to_string((unsigned int)(std::end(*Fe) - std::begin(*Fe))) +
                                         " but must be of size " + std::to_string(dim * dim))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(C_end - C_begin) == dim * dim,
                                     "The storage array for C is of size " +
                                         std::to_string((unsigned int)(C_end - C_begin)) + " but must be of size " +
                                         std::to_string(dim * dim))

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                for (unsigned int J = 0; J < dim; ++J) {
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
    template <typename T>
    T HyperelasticBase::compute_I2(const bool isPrevious) {
        constexpr unsigned int dim = 3;  // TODO: Get this from ResidualBase

        T I1 = compute_I1<T>(isPrevious);

        T I2 = 0.5 * I1 * I1;

        std::array<T, dim * dim> C{};
        compute_right_cauchy_green_deformation_tensor(isPrevious, std::begin(C), std::end(C));

        I2 -= 0.5 * std::inner_product(std::begin(C), std::end(C), std::begin(C), T());

        return I2;
    }

    /*!
     * Compute the derivative of the second invariant of the deformation with respect to the elastic deformation
     * gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param dI2dFe_begin: The starting iterator of the derivative of the second invariant of the deformation with
     * respect to the elastic deformation gradient
     * \param dI2dFe_end: The stopping iterator of the derivative of the second invariant of the deformation with
     * respect to the elastic deformation gradient
     */
    template <class dI2dFe_iter>
    void HyperelasticBase::compute_dI2dFe(const bool isPrevious, dI2dFe_iter dI2dFe_begin, dI2dFe_iter dI2dFe_end) {
        using dI2dFe_type = typename std::iterator_traits<dI2dFe_iter>::value_type;

        constexpr unsigned int dim = 3;  // TODO: Get this value from ResidualBase

        const secondOrderTensor *Fe;

        if (isPrevious) {
            Fe = get_previousFe();
        } else {
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dI2dFe_end - dI2dFe_begin) == dim * dim,
                                     "The size of dI2dFe is " +
                                         std::to_string((unsigned int)(dI2dFe_end - dI2dFe_begin)) +
                                         " but it must be " + std::to_string(dim * dim));

        dI2dFe_type I1 = compute_I1<dI2dFe_type>(isPrevious);

        std::transform(std::begin(*Fe), std::end(*Fe), dI2dFe_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, I1));

        std::array<dI2dFe_type, dim * dim> C{};
        compute_right_cauchy_green_deformation_tensor(isPrevious, std::begin(C), std::end(C));

        for (unsigned int a = 0; a < dim; ++a) {
            for (unsigned int I = 0; I < dim; ++I) {
                for (unsigned int A = 0; A < dim; ++A) {
                    *(dI2dFe_begin + dim * a + A) -= (*Fe)[dim * a + I] * C[dim * I + A];
                }
            }
        }

        std::transform(dI2dFe_begin, dI2dFe_end, dI2dFe_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, 2));
    }

    /*!
     * Compute the second derivative of the second invariant of the deformation with respect to the elastic deformation
     * gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param d2I2dFe2_begin: The starting iterator of the second derivative of the second invariant of the deformation
     * with respect to the elastic deformation gradient
     * \param d2I2dFe2_end: The stopping iterator of the second derivative of the second invariant of the deformation
     * with respect to the elastic deformation gradient
     */
    template <class d2I2dFe2_iter>
    void HyperelasticBase::compute_d2I2dFe2(const bool isPrevious, d2I2dFe2_iter d2I2dFe2_begin,
                                            d2I2dFe2_iter d2I2dFe2_end) {
        using d2I2dFe2_type = typename std::iterator_traits<d2I2dFe2_iter>::value_type;

        constexpr unsigned int dim = 3;  // TODO: Get this value from ResidualBase

        const secondOrderTensor *Fe;

        if (isPrevious) {
            Fe = get_previousFe();
        } else {
            Fe = get_Fe();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(d2I2dFe2_end - d2I2dFe2_begin) == dim * dim * dim * dim,
                                     "The size of d2I2dFe2 is " +
                                         std::to_string((unsigned int)(d2I2dFe2_end - d2I2dFe2_begin)) +
                                         " but it must be " + std::to_string(dim * dim * dim * dim));

        d2I2dFe2_type I1 = compute_I1<d2I2dFe2_type>(isPrevious);

        std::array<d2I2dFe2_type, dim * dim> C{};
        compute_right_cauchy_green_deformation_tensor(isPrevious, std::begin(C), std::end(C));
        std::fill(d2I2dFe2_begin, d2I2dFe2_end, d2I2dFe2_type());

        for (unsigned int aA = 0; aA < dim * dim; ++aA) {
            *(d2I2dFe2_begin + dim * dim * aA + aA) += I1;
        }

        for (unsigned int a = 0; a < dim; ++a) {
            for (unsigned int A = 0; A < dim; ++A) {
                for (unsigned int b = 0; b < dim; ++b) {
                    *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * a + b) -= C[dim * A + b];
                    for (unsigned int B = 0; B < dim; ++B) {
                        *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * b + B) +=
                            2 * (*Fe)[dim * a + A] * (*Fe)[dim * b + B];
                        *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * b + B) -=
                            (*Fe)[dim * a + B] * (*Fe)[dim * b + A];
                        *(d2I2dFe2_begin + dim * dim * dim * a + dim * dim * A + dim * b + A) -=
                            (*Fe)[dim * a + B] * (*Fe)[dim * b + B];
                    }
                }
            }
        }

        std::transform(d2I2dFe2_begin, d2I2dFe2_end, d2I2dFe2_begin,
                       std::bind(std::multiplies<>(), std::placeholders::_1, 2));
    }

    /*!
     * Compute the isochoric first invariant of the deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     */
    template <typename T>
    inline T HyperelasticBase::compute_Ibar1(const bool isPrevious) {
        const floatType *Je;

        if (isPrevious) {
            Je = get_previousJe();
        } else {
            Je = get_Je();
        }

        T I1 = compute_I1<T>(isPrevious);

        T Ibar1 = I1 / std::pow(*Je, 2. / 3.);

        return Ibar1;
    }

    /*!
     * Compute the derivative of the isochoric first invariant of the deformation with respect to the elastic
     * deformation gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param dIbar1dFe_begin: The starting iterator of the derivative of the isochoric first invariant of the
     * deformation with respect to the elastic deformation gradient
     * \param dIbar1dFe_end: The stopping iterator of the derivative of the isochoric first invariant of the deformation
     * with respect to the elastic deformation gradient
     */
    template <class dIbar1dFe_iter>
    void HyperelasticBase::compute_dIbar1dFe(const bool isPrevious, dIbar1dFe_iter dIbar1dFe_begin,
                                             dIbar1dFe_iter dIbar1dFe_end) {
        using dIbar1dFe_type = typename std::iterator_traits<dIbar1dFe_iter>::value_type;

        constexpr unsigned int dim = 3;  // TODO: Get this value from ResidualBase

        constexpr double f = 2. / 3;

        const floatType *Je;

        const secondOrderTensor *dJedFe;

        if (isPrevious) {
            Je     = get_previousJe();
            dJedFe = get_dPreviousJedPreviousFe();
        } else {
            Je     = get_Je();
            dJedFe = get_dJedFe();
        }

        dIbar1dFe_type I1 = compute_I1<dIbar1dFe_type>(isPrevious);

        compute_dI1dFe(isPrevious, dIbar1dFe_begin, dIbar1dFe_end);

        for (unsigned int iI = 0; iI < dim * dim; ++iI) {
            *(dIbar1dFe_begin + iI) -= f * I1 * (*dJedFe)[iI] / (*Je);
        }

        std::transform(dIbar1dFe_begin, dIbar1dFe_end, dIbar1dFe_begin,
                       std::bind(std::divides<>(), std::placeholders::_1, std::pow(*Je, f)));
    }

    /*!
     * Compute the second derivative of the isochoric first invariant of the deformation with respect to the elastic
     * deformation gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param d2Ibar1dFe2_begin: The starting iterator of the second derivative of the isochoric second invariant of the
     * deformation with respect to the elastic deformation gradient
     * \param d2Ibar1dFe2_end: The stopping iterator of the second derivative of the isochoric second invariant of the
     * deformation with respect to the elastic deformation gradient
     */
    template <class d2Ibar1dFe2_iter>
    void HyperelasticBase::compute_d2Ibar1dFe2(const bool isPrevious, d2Ibar1dFe2_iter d2Ibar1dFe2_begin,
                                               d2Ibar1dFe2_iter d2Ibar1dFe2_end) {
        using d2Ibar1dFe2_type = typename std::iterator_traits<d2Ibar1dFe2_iter>::value_type;

        constexpr unsigned int dim = 3;  // TODO: Get this value from ResidualBase

        constexpr double f = 2. / 3.;

        const floatType *Je;

        const secondOrderTensor *dJedFe;

        const fourthOrderTensor *d2JedFe2;

        if (isPrevious) {
            Je       = get_previousJe();
            dJedFe   = get_dPreviousJedPreviousFe();
            d2JedFe2 = get_d2PreviousJedPreviousFe2();
        } else {
            Je       = get_Je();
            dJedFe   = get_dJedFe();
            d2JedFe2 = get_d2JedFe2();
        }

        d2Ibar1dFe2_type I1 = compute_I1<d2Ibar1dFe2_type>(isPrevious);

        std::array<d2Ibar1dFe2_type, dim * dim>             dI1dFe{};
        std::array<d2Ibar1dFe2_type, dim * dim>             dIbar1dFe{};
        std::array<d2Ibar1dFe2_type, dim * dim * dim * dim> d2I1dFe2{};

        compute_dI1dFe(isPrevious, std::begin(dI1dFe), std::end(dI1dFe));
        compute_dIbar1dFe(isPrevious, std::begin(dIbar1dFe), std::end(dIbar1dFe));
        compute_d2I1dFe2(std::begin(d2I1dFe2), std::end(d2I1dFe2));

        std::fill(d2Ibar1dFe2_begin, d2Ibar1dFe2_end, d2Ibar1dFe2_type());

        for (unsigned int aA = 0; aA < dim * dim; ++aA) {
            for (unsigned int bB = 0; bB < dim * dim; ++bB) {
                *(d2Ibar1dFe2_begin + dim * dim * aA + bB) +=
                    (d2I1dFe2[dim * dim * aA + bB] - f / (*Je) * (*dJedFe)[aA] * dI1dFe[bB] +
                     f * I1 / ((*Je) * (*Je)) * (*dJedFe)[aA] * (*dJedFe)[bB] -
                     f * I1 / (*Je) * (*d2JedFe2)[dim * dim * aA + bB]) /
                        std::pow((*Je), f) -
                    f / (*Je) * dIbar1dFe[aA] * (*dJedFe)[bB];
            }
        }
    }

    /*!
     * Compute the isochoric second invariant of the deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     */
    template <typename T>
    inline T HyperelasticBase::compute_Ibar2(const bool isPrevious) {
        constexpr double f = 4. / 3;

        const floatType *Je;

        if (isPrevious) {
            Je = get_previousJe();
        } else {
            Je = get_Je();
        }

        T I2 = compute_I2<T>(isPrevious);

        T Ibar2 = I2 / std::pow(*Je, f);

        return Ibar2;
    }
    /*!
     * Compute the derivative of the isochoric second invariant of the deformation with respect to the elastic
     * deformation gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param dIbar2dFe_begin: The starting iterator of the derivative of the isochoric second invariant of the
     * deformation with respect to the elastic deformation gradient
     * \param dIbar2dFe_end: The stopping iterator of the derivative of the isochoric second invariant of the
     * deformation with respect to the elastic deformation gradient
     */
    template <class dIbar2dFe_iter>
    void HyperelasticBase::compute_dIbar2dFe(const bool isPrevious, dIbar2dFe_iter dIbar2dFe_begin,
                                             dIbar2dFe_iter dIbar2dFe_end) {
        using dIbar2dFe_type = typename std::iterator_traits<dIbar2dFe_iter>::value_type;

        constexpr unsigned int dim = 3;  // TODO: Get this value from ResidualBase

        constexpr double f = 4. / 3;

        const floatType *Je;

        const secondOrderTensor *dJedFe;

        if (isPrevious) {
            Je     = get_previousJe();
            dJedFe = get_dPreviousJedPreviousFe();
        } else {
            Je     = get_Je();
            dJedFe = get_dJedFe();
        }

        dIbar2dFe_type I2 = compute_I2<dIbar2dFe_type>(isPrevious);

        compute_dI2dFe(isPrevious, dIbar2dFe_begin, dIbar2dFe_end);

        for (unsigned int iI = 0; iI < dim * dim; ++iI) {
            *(dIbar2dFe_begin + iI) -= f * I2 * (*dJedFe)[iI] / (*Je);
        }

        std::transform(dIbar2dFe_begin, dIbar2dFe_end, dIbar2dFe_begin,
                       std::bind(std::divides<>(), std::placeholders::_1, std::pow(*Je, f)));
    }

    /*!
     * Compute the second derivative of the isochoric second invariant of the deformation with respect to the elastic
     * deformation gradient
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true) value
     * \param d2Ibar2dFe2_begin: The starting iterator of the second derivative of the isochoric second invariant of the
     * deformation with respect to the elastic deformation gradient
     * \param d2Ibar2dFe2_end: The stopping iterator of the second derivative of the isochoric second invariant of the
     * deformation with respect to the elastic deformation gradient
     */
    template <class d2Ibar2dFe2_iter>
    void HyperelasticBase::compute_d2Ibar2dFe2(const bool isPrevious, d2Ibar2dFe2_iter d2Ibar2dFe2_begin,
                                               d2Ibar2dFe2_iter d2Ibar2dFe2_end) {
        using d2Ibar2dFe2_type = typename std::iterator_traits<d2Ibar2dFe2_iter>::value_type;

        constexpr unsigned int dim = 3;  // TODO: Get this value from ResidualBase

        constexpr double f = 4. / 3.;

        const floatType *Je;

        const secondOrderTensor *dJedFe;

        const fourthOrderTensor *d2JedFe2;

        if (isPrevious) {
            Je       = get_previousJe();
            dJedFe   = get_dPreviousJedPreviousFe();
            d2JedFe2 = get_d2PreviousJedPreviousFe2();
        } else {
            Je       = get_Je();
            dJedFe   = get_dJedFe();
            d2JedFe2 = get_d2JedFe2();
        }

        d2Ibar2dFe2_type I2 = compute_I2<d2Ibar2dFe2_type>(isPrevious);

        std::array<d2Ibar2dFe2_type, dim * dim>             dI2dFe{};
        std::array<d2Ibar2dFe2_type, dim * dim>             dIbar2dFe{};
        std::array<d2Ibar2dFe2_type, dim * dim * dim * dim> d2I2dFe2{};

        compute_dI2dFe(isPrevious, std::begin(dI2dFe), std::end(dI2dFe));
        compute_dIbar2dFe(isPrevious, std::begin(dIbar2dFe), std::end(dIbar2dFe));
        compute_d2I2dFe2(isPrevious, std::begin(d2I2dFe2), std::end(d2I2dFe2));

        std::fill(d2Ibar2dFe2_begin, d2Ibar2dFe2_end, d2Ibar2dFe2_type());

        for (unsigned int aA = 0; aA < dim * dim; ++aA) {
            for (unsigned int bB = 0; bB < dim * dim; ++bB) {
                *(d2Ibar2dFe2_begin + dim * dim * aA + bB) +=
                    (d2I2dFe2[dim * dim * aA + bB] - f / (*Je) * (*dJedFe)[aA] * dI2dFe[bB] +
                     f * I2 / ((*Je) * (*Je)) * (*dJedFe)[aA] * (*dJedFe)[bB] -
                     f * I2 / (*Je) * (*d2JedFe2)[dim * dim * aA + bB]) /
                        std::pow((*Je), f) -
                    f / (*Je) * dIbar2dFe[aA] * (*dJedFe)[bB];
            }
        }
    }

}  // namespace tardigradeHydra
