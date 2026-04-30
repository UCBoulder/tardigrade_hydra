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
    template <class configuration>
    template <class v_type>
    void HydraBase<configuration>::addToFailureOutput(const v_type &v, bool add_endline) {
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
    template <class configuration>
    template <class v_iterator>
    void HydraBase<configuration>::addToFailureOutput(const v_iterator &v_begin, const v_iterator &v_end,
                                                      bool add_endline) {
        for (auto v = v_begin; v != v_end; ++v) {
            _failure_output << *v << ", ";
        }
        if (add_endline) {
            _failure_output << "\n";
        }
    }

    /*!
     * Add a string to the failure output string
     *
     * \param &value: The string to append to the output
     * \param add_endline: A boolean for if the endline character should be added after the value
     */
    template <class configuration>
    void HydraBase<configuration>::addToFailureOutput(const std::string &value, bool add_endline) {
        addToFailureOutput<std::string>(value, add_endline);
    }

    /*!
     * Add a floatVector to the output string
     *
     * \param &value: The vector to add to the output string
     * \param add_endline: A boolean for if the endline character should be added after the value
     */
    template <class configuration>
    void HydraBase<configuration>::addToFailureOutput(const floatVector &value, bool add_endline) {
        addToFailureOutput(std::begin(value), std::end(value), add_endline);
    }

    /*!
     * Add a vector of booleans to the output string
     *
     * \param &value: The vector to add to the output string
     * \param add_endline: A boolean for if the endline character should be added after the value
     */
    template <class configuration>
    void HydraBase<configuration>::addToFailureOutput(const std::vector<bool> &value, bool add_endline) {
        addToFailureOutput(std::begin(value), std::end(value), add_endline);
    }

    /*!
     * Add a floating point value to the output string
     *
     * \param &value: The value to add to the output string
     * \param add_endline: A boolean for if the endline character should be added after the value
     */
    template <class configuration>
    void HydraBase<configuration>::addToFailureOutput(const floatType &value, bool add_endline) {
        addToFailureOutput<floatType>(value, add_endline);
    }

}  // namespace tardigradeHydra
