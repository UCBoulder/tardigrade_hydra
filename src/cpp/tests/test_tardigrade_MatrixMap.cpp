/**
 * \file test_tardigrade_MatrixMap.cpp
 *
 * Tests for tardigrade_MatrixMap
 */

#include <tardigrade_MatrixMap.h>

#define BOOST_TEST_MODULE test_tardigrade_MatrixMap
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

bool tolerantCheck(const std::vector<double> &v1, const std::vector<double> &v2, double eps = 1e-6, double tol = 1e-9) {
    if (v1.size() != v2.size()) {
        return false;
    }

    BOOST_CHECK(v1.size() == v2.size());

    const unsigned int len = v1.size();

    for (unsigned int i = 0; i < len; i++) {
        if (std::fabs(v1[i]) < tol) {
            if (std::fabs(v1[i] - v2[i]) > eps) {
                return false;
            }

            BOOST_CHECK(std::fabs(v1[i] - v2[i]) <= eps);

        } else {
            if ((std::fabs(v1[i] - v2[i]) / std::fabs(v1[i]) > eps) ||
                (std::fabs(v1[i] - v2[i]) / std::fabs(v2[i]) > eps)) {
                return false;
            }

            BOOST_TEST(v1[i] == v2[i]);
        }
    }

    return true;
}

BOOST_AUTO_TEST_CASE(test_matrixMap, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    const std::vector<double> A = {+1.000000000e-01, +1.473684211e-01, +1.947368421e-01, +2.421052632e-01,
                                   +2.894736842e-01, +3.368421053e-01, +3.842105263e-01, +4.315789474e-01,
                                   +4.789473684e-01, +5.263157895e-01, +5.736842105e-01, +6.210526316e-01,
                                   +6.684210526e-01, +7.157894737e-01, +7.631578947e-01, +8.105263158e-01,
                                   +8.578947368e-01, +9.052631579e-01, +9.526315789e-01, +1.000000000e+00};

    const std::vector<double> B = {
        +1.200000000e+00, +1.227586207e+00, +1.255172414e+00, +1.282758621e+00, +1.310344828e+00, +1.337931034e+00,
        +1.365517241e+00, +1.393103448e+00, +1.420689655e+00, +1.448275862e+00, +1.475862069e+00, +1.503448276e+00,
        +1.531034483e+00, +1.558620690e+00, +1.586206897e+00, +1.613793103e+00, +1.641379310e+00, +1.668965517e+00,
        +1.696551724e+00, +1.724137931e+00, +1.751724138e+00, +1.779310345e+00, +1.806896552e+00, +1.834482759e+00,
        +1.862068966e+00, +1.889655172e+00, +1.917241379e+00, +1.944827586e+00, +1.972413793e+00, +2.000000000e+00};

    const std::vector<double> b = {+1.200000000e+00, +1.288888889e+00, +1.377777778e+00, +1.466666667e+00,
                                   +1.555555556e+00, +1.644444444e+00, +1.733333333e+00, +1.822222222e+00,
                                   +1.911111111e+00, +2.000000000e+00};

    std::vector<double> C_answer = {+5.247549909e+00, +5.333938294e+00, +5.420326679e+00,
                                    +1.269582577e+01, +1.291288566e+01, +1.312994555e+01};

    std::vector<double> c_answer = {+5.357894737e+00, +1.293684211e+01};

    // Fixed size map test
    {
        std::vector<double> C_result(6, 0);

        // Fixed size map test
        auto A_map = tardigradeHydra::getFixedSizeMatrixMap<double, 2, 10>(A.data());
        auto B_map = tardigradeHydra::getFixedSizeMatrixMap<double, 10, 3>(B.data());
        auto C_map = tardigradeHydra::getFixedSizeMatrixMap<double, 2, 3>(C_result.data());

        C_map = (A_map * B_map).eval();

        BOOST_TEST(C_result == C_answer, CHECK_PER_ELEMENT);
    }

    // Dynamic column size map test
    {
        std::vector<double> C_result(6, 0);

        // Fixed size map test
        auto A_map = tardigradeHydra::getDynamicColumnSizeMatrixMap<double, 2>(A.data(), 10);
        auto B_map = tardigradeHydra::getDynamicColumnSizeMatrixMap<double, 10>(B.data(), 3);
        auto C_map = tardigradeHydra::getDynamicColumnSizeMatrixMap<double, 2>(C_result.data(), 3);

        C_map = (A_map * B_map).eval();

        BOOST_TEST(C_result == C_answer, CHECK_PER_ELEMENT);
    }

    // Dynamic size map test
    {
        std::vector<double> C_result(6, 0);

        // Fixed size map test
        auto A_map = tardigradeHydra::getDynamicSizeMatrixMap<double>(A.data(), 2, 10);
        auto B_map = tardigradeHydra::getDynamicSizeMatrixMap<double>(B.data(), 10, 3);
        auto C_map = tardigradeHydra::getDynamicSizeMatrixMap<double>(C_result.data(), 2, 3);

        C_map = (A_map * B_map).eval();

        BOOST_TEST(C_result == C_answer, CHECK_PER_ELEMENT);
    }

    // Fixed size vector map test
    {
        std::vector<double> c_result(2, 0);

        // Fixed size map test
        auto A_map = tardigradeHydra::getFixedSizeMatrixMap<double, 2, 10>(A.data());
        auto b_map = tardigradeHydra::getFixedSizeVectorMap<double, 10>(b.data());
        auto c_map = tardigradeHydra::getFixedSizeVectorMap<double, 2>(c_result.data());

        c_map = (A_map * b_map).eval();

        BOOST_TEST(c_result == c_answer, CHECK_PER_ELEMENT);
    }

    // Dynamic size vector map test
    {
        std::vector<double> c_result(2, 0);

        // Fixed size map test
        auto A_map = tardigradeHydra::getFixedSizeMatrixMap<double, 2, 10>(A.data());
        auto b_map = tardigradeHydra::getFixedSizeVectorMap<double, 10>(b.data());
        auto c_map = tardigradeHydra::getDynamicSizeVectorMap<double>(c_result.data(), 2);

        c_map = (A_map * b_map).eval();

        BOOST_TEST(c_result == c_answer, CHECK_PER_ELEMENT);
    }
}
