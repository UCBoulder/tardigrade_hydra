/**
 ******************************************************************************
 * \file tardigrade_DeformationEvolutionBase.tpp
 ******************************************************************************
 * The base class for defining an evolution equation for deformation
 ******************************************************************************
 */

namespace tardigradeHydra {

    /*!
     * Compute the deformation
     *
     * \param &dt: The change in time
     * \param &Lt_begin: The starting iterator of the previous velocity gradient
     * \param &Lt_end: The stopping iterator of the previous velocity gradient
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Ft_begin: The starting iterator of the previous deformation
     * \param &Ft_end: The stopping iterator of the previous deformation
     * \param &Ftp1_begin: The starting iterator of the deformation
     * \param &Ftp1_end: The stopping iterator of the deformation
     */
    template<class container, int size>
    template<
    typename dt, class Lt_iterator, class Ltp1_iterator, class Ft_iterator, class Ftp1_iterator
    >
    void DeformationEvolutionBase<container, size>::computeDeformation(
        const dt_type &dt,
        const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end,
        const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
        const Ft_iterator &Ft_begin, const Ft_iterator &Ft_end,
        Ftp1_iterator Ftp1_begin, Ftp1_iterator Ftp1_end){

        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;
        using Lt_type   = typename std::iterator_traits<Lt_iterator>::value_type;
        using Ft_type   = typename std::iterator_traits<Ft_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Lt_end - Lt_begin) == size*size,
                "The size of Lt is " + std::to_string((unsigned int)(Lt_end - Lt_begin)) + " but it should be " + std::to_string(size*size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ltp1_end - Ltp1_begin) == size*size,
                "The size of Ltp1 is " + std::to_string((unsigned int)(Ltp1_end - Ltp1_begin)) + " but it should be " + std::to_string(size*size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ft_end - Ft_begin) == size*size,
                "The size of Ft is " + std::to_string((unsigned int)(Ft_end - Ft_begin)) + " but it should be " + std::to_string(size*size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ftp1_end - Ftp1_begin) == size*size,
                "The size of Ftp1 is " + std::to_string((unsigned int)(Ftp1_end - Ftp1_begin)) + " but it should be " + std::to_string(size*size))

        //Form the left and right hand sides
        std::array<Ltp1_type, size*size> LHS;
        std::transform(Ltp1_begin, Ltp1_end, std::begin(LHS),
                std::bind(std::multiplies<>(), std::placeholders::_1, -dt * integration_parameter));
        std::array<Lt_type, size> previous_dF{};
        std::array<Ft_type, size> RHS{};
        std::transform(Lt_begin, Lt_end, std::begin(previous_dF),
                std::bind(std::multiplies<>(), std::placeholders::_1, dt * (1 - integration_parameter)));

        for ( unsigned int i = 0; i < size, ++i ){
            LHS[size*i+i] += 1;
            previous_dF[size*i+i] += 1;
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int I = 0; I < size; ++I ){
                    RHS[size*i+I] += previous_dF[dim * i + j] * (*(Ft_begin + dim* j + I));
                }
            }
        }
    }
}
