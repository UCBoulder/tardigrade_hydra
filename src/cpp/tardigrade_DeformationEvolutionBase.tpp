/**
 ******************************************************************************
 * \file tardigrade_DeformationEvolutionBase.tpp
 ******************************************************************************
 * The base class for defining an evolution equation for deformation
 ******************************************************************************
 */


    /*!
     * Compute the deformation
     *
     * \param &Lt_begin: The starting iterator of the previous velocity gradient
     * \param &Lt_end: The stopping iterator of the previous velocity gradient
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Ft_begin: The starting iterator of the previous deformation
     * \param &Ft_end: The stopping iterator of the previous deformation
     * \param &Ftp1_begin: The starting iterator of the deformation
     * \param &Ftp1_end: The stopping iterator of the deformation
     */
    template<container>
    template<
    class Lt_iterator, class Ltp1_iterator, class Ft_iterator, class Ftp1_iterator
    >
    void DeformationEvolutionBase<container>::computeDeformation(
        const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end,
        const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
        const Ft_iterator &Ft_begin, const Ft_iterator &Ft_end,
        Ftp1_iterator Ftp1_begin, Ftp1_iterator Ftp1_end){

    }
