/**
 ******************************************************************************
 * \file tardigrade_hydra.tpp
 ******************************************************************************
 * Template definitions for tardigrade_hydra
 ******************************************************************************
 */

namespace tardigradeHydra {

    /*!
     * Add a general non-iterable object to the output string
     *
     * \param &v: The value to add
     * \param add_endline: Whether to append a newline character after the value
     */
    template <class v_type>
    void hydraBase::addToFailureOutput(const v_type &v, bool add_endline) {
        _failure_output << v;
        if (add_endline) {
            _failure_output << "\n";
        }
    }

    /*!
     * Add a general iterable object to the output string
     *
     * \param &v_begin: The starting iterator of the value vector
     * \param &v_end: The stopping iterator of the value vector
     * \param add_endline: Whether to append a newline character after the value
     */
    template <class v_iterator>
    void hydraBase::addToFailureOutput(const v_iterator &v_begin, const v_iterator &v_end, bool add_endline) {
        for (auto v = v_begin; v != v_end; ++v) {
            _failure_output << *v << ", ";
        }
        if (add_endline) {
            _failure_output << "\n";
        }
    }
}
