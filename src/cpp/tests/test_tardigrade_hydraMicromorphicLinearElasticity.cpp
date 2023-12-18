/**
  * \file test_tardigrade_hydraMicromorphicLinearElasticity.cpp
  *
  * Tests for tardigrade_hydraMicromorphicLinearElasticity
  */

#include<tardigrade_hydraMicromorphicLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphicLinearElasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::micromorphicLinearElasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::micromorphicLinearElasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::micromorphicLinearElasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

typedef tardigradeHydra::micromorphicLinearElasticity::variableType variableType; //!< Redefinition of the variable type
typedef tardigradeHydra::micromorphicLinearElasticity::variableVector variableVector; //!< Redefinition of the vector of variable types
typedef tardigradeHydra::micromorphicLinearElasticity::variableMatrix variableMatrix; //!< Redefinition of the matrix of variable types

typedef tardigradeHydra::micromorphicLinearElasticity::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicLinearElasticity::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicLinearElasticity::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicLinearElasticity::constantVector constantVector; //!< Redefinition of the vector of constants
typedef tardigradeHydra::micromorphicLinearElasticity::constantMatrix constantMatrix; //!< Redefinition of the matrix of constants

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( testComputeReferenceHigherOrderStress ){
    /*!
     * Test the computation of the higher order stress in the reference configuration.
     *
     */

    variableVector Gamma = { -0.31120922, -0.3563267 , -0.36573233, -0.0771914 , -0.24252804,
                             -0.4738459 , -0.35937075, -0.01781817, -0.17465609, -0.82225557,
                             -0.36719542, -0.86494826, -0.92750732, -1.18214541, -1.00423785,
                             -0.43125133, -0.19543115, -0.49736256, -0.67098335, -0.98433811,
                             -1.0183107 , -0.12645195, -0.79818076, -1.23318541, -0.95577138,
                              0.1274431 , -0.47648617 };

    parameterVector C = { 0.73168423, 0.20787282, 0.44597068, 0.13971472, 0.67623962,
                          0.01024113, 0.26235898, 0.52638878, 0.55822266, 0.79169357,
                          0.79556295, 0.26274848, 0.22158735, 0.77447856, 0.09905053,
                          0.24537506, 0.82422833, 0.13912553, 0.13704714, 0.60418098,
                          0.97916951, 0.96975567, 0.15156735, 0.33820056, 0.40984014,
                          0.03316931, 0.07805217, 0.09086063, 0.0587449 , 0.93973718,
                          0.2088402 , 0.27030923, 0.90893679, 0.45495913, 0.38558114,
                          0.89599555, 0.99117823, 0.68663521, 0.47807759, 0.47462775,
                          0.63545614, 0.76359775, 0.58543208, 0.67889697, 0.23400169,
                          0.90355814, 0.81277492, 0.27582141, 0.34052255, 0.09211961,
                          0.38139711, 0.40588835, 0.72562255, 0.46843548, 0.53552493,
                          0.10983976, 0.70222301, 0.33678326, 0.5755877 , 0.25931049,
                          0.4508404 , 0.61309158, 0.93828205, 0.54132972, 0.74656135,
                          0.66254458, 0.79255496, 0.15105507, 0.73739719, 0.8699253 ,
                          0.30037259, 0.18970023, 0.7351555 , 0.19264707, 0.66209567,
                          0.96806941, 0.64975045, 0.54480589, 0.11455497, 0.48495469,
                          0.76967885, 0.23472007, 0.43837476, 0.95659971, 0.66906068,
                          0.06233808, 0.89998454, 0.61301419, 0.83415149, 0.08989232,
                          0.68242323, 0.5311455 , 0.04401578, 0.64268437, 0.17462098,
                          0.91385774, 0.3254543 , 0.05379172, 0.77842952, 0.53732987,
                          0.83781747, 0.7615337 , 0.62005118, 0.73293008, 0.82346972,
                          0.2179367 , 0.36317936, 0.55619014, 0.09312705, 0.4694652 ,
                          0.69716073, 0.02440034, 0.53107043, 0.97073677, 0.87477045,
                          0.68304994, 0.78714746, 0.23263201, 0.42970537, 0.72939955,
                          0.12454156, 0.63073136, 0.31116734, 0.39634253, 0.1281059 ,
                          0.23014438, 0.13811186, 0.23720232, 0.16020445, 0.10732273,
                          0.90941219, 0.89549355, 0.71971555, 0.01863396, 0.63512323,
                          0.57834974, 0.74377666, 0.67809585, 0.69763216, 0.75809662,
                          0.34523363, 0.42721884, 0.06343445, 0.4707797 , 0.91267953,
                          0.9522137 , 0.429135  , 0.94506974, 0.42002389, 0.48045774,
                          0.1877623 , 0.93541748, 0.10358528, 0.90585229, 0.94482345,
                          0.85486823, 0.71201071, 0.82729667, 0.67685002, 0.89074951,
                          0.13603059, 0.29549921, 0.5829021 , 0.85710379, 0.33597495,
                          0.2635317 , 0.54822056, 0.13518258, 0.07510343, 0.57277576,
                          0.66026008, 0.10590873, 0.40988651, 0.73181046, 0.21849923,
                          0.68193615, 0.4861005 , 0.90062638, 0.49503759, 0.53109181,
                          0.31197913, 0.8260051 , 0.56845431, 0.20510746, 0.48927707,
                          0.84796951, 0.57021869, 0.32869802, 0.00649644, 0.89085066,
                          0.58793337, 0.13725509, 0.49166181, 0.79467837, 0.55550476,
                          0.86168924, 0.26284446, 0.34931772, 0.69039842, 0.04226658,
                          0.91252659, 0.8532767 , 0.15745086, 0.11244899, 0.35188228,
                          0.66119509, 0.88971845, 0.90199259, 0.53564388, 0.08103036,
                          0.89537074, 0.43988547, 0.39234971, 0.90744335, 0.87819375,
                          0.25940274, 0.48165619, 0.08404158, 0.16900508, 0.20502448,
                          0.00336955, 0.94376888, 0.89722214, 0.06817336, 0.35272289,
                          0.34452052, 0.23363246, 0.79650105, 0.8107239 , 0.94490429,
                          0.26741852, 0.87105166, 0.25525768, 0.26586211, 0.6449152 ,
                          0.10839033, 0.6871309 , 0.59008043, 0.07558712, 0.99527881,
                          0.13052048, 0.81075174, 0.38967993, 0.25408067, 0.78035165,
                          0.48123955, 0.97775619, 0.50408867, 0.51411035, 0.17947261,
                          0.99740746, 0.84538866, 0.62373254, 0.38782162, 0.55585207,
                          0.24743969, 0.25980163, 0.50272755, 0.76170535, 0.338618  ,
                          0.33580793, 0.58798537, 0.13328799, 0.01026525, 0.13839967,
                          0.31984267, 0.72693472, 0.84737434, 0.97859975, 0.61637914,
                          0.23018791, 0.89805651, 0.69772024, 0.7491404 , 0.3818782 ,
                          0.50000777, 0.71398283, 0.29910862, 0.36270529, 0.93178041,
                          0.03156497, 0.57674924, 0.13573152, 0.59758916, 0.47467419,
                          0.21707829, 0.36305461, 0.58480959, 0.18659161, 0.20999611,
                          0.23732489, 0.11099326, 0.49055309, 0.52547794, 0.01722654,
                          0.19637688, 0.03560497, 0.89843994, 0.34941756, 0.10796044,
                          0.07166564, 0.67297414, 0.34139877, 0.56321003, 0.13224438,
                          0.58789568, 0.05265614, 0.93254668, 0.41326988, 0.67692951,
                          0.27922074, 0.38788297, 0.24478052, 0.29147   , 0.80741949,
                          0.67936156, 0.7442339 , 0.00343505, 0.97756934, 0.02268554,
                          0.56302353, 0.10718293, 0.11474464, 0.10465633, 0.04846449,
                          0.33695467, 0.43787266, 0.15092164, 0.80017919, 0.40017523,
                          0.40391072, 0.65025117, 0.3018835 , 0.15825793, 0.02963411,
                          0.85526189, 0.29678796, 0.23667277, 0.4013067 , 0.76988912,
                          0.06110263, 0.66297631, 0.79827956, 0.70776264, 0.07467447,
                          0.89814767, 0.00308201, 0.7823472 , 0.38646676, 0.75957091,
                          0.47684411, 0.7398732 , 0.09206989, 0.02529722, 0.1859329 ,
                          0.8380139 , 0.33920514, 0.0689636 , 0.07697459, 0.88068415,
                          0.1229827 , 0.89652486, 0.13968141, 0.38800301, 0.7728669 ,
                          0.16905682, 0.47354036, 0.90279152, 0.62008568, 0.82696116,
                          0.43869547, 0.94940345, 0.12938034, 0.20523408, 0.11727954,
                          0.54592836, 0.82115919, 0.96255349, 0.04999854, 0.15256932,
                          0.62537849, 0.15516518, 0.25683723, 0.85702076, 0.7925628 ,
                          0.46399241, 0.10106241, 0.73089281, 0.46200846, 0.24160109,
                          0.01349364, 0.94349643, 0.70886053, 0.06715038, 0.95042753,
                          0.38413263, 0.61285658, 0.97690412, 0.07900655, 0.21037925,
                          0.03351281, 0.36733596, 0.05601802, 0.50752553, 0.62088055,
                          0.94638543, 0.31649186, 0.19788369, 0.59813263, 0.31156879,
                          0.84129622, 0.18756002, 0.80252603, 0.44583102, 0.08424927,
                          0.8055779 , 0.89467745, 0.32244817, 0.5244238 , 0.38246742,
                          0.17552342, 0.09374914, 0.02755403, 0.86455687, 0.45570292,
                          0.58901182, 0.11888058, 0.65051228, 0.9634849 , 0.72370701,
                          0.6882061 , 0.30785926, 0.13060746, 0.29416438, 0.87322017,
                          0.26415365, 0.41275749, 0.44246432, 0.53266346, 0.0943344 ,
                          0.30480514, 0.37707017, 0.41691054, 0.94780656, 0.48190006,
                          0.55313378, 0.34750865, 0.4482111 , 0.62723585, 0.72810975,
                          0.05039657, 0.27579548, 0.03891394, 0.25236345, 0.53330415,
                          0.73508523, 0.07895291, 0.17533722, 0.35439847, 0.35161594,
                          0.56198773, 0.09715776, 0.65962064, 0.93017981, 0.69473252,
                          0.36962137, 0.01256946, 0.76358204, 0.71278475, 0.90981899,
                          0.12157694, 0.72489524, 0.35040793, 0.67236183, 0.70497179,
                          0.65668193, 0.42266533, 0.45514152, 0.13569985, 0.70177277,
                          0.47524764, 0.61467852, 0.7338517 , 0.54271988, 0.0514032 ,
                          0.75535566, 0.65658309, 0.01198156, 0.99577374, 0.22789359,
                          0.54174149, 0.25282717, 0.08002466, 0.26660684, 0.07484457,
                          0.87263561, 0.39011271, 0.93760461, 0.80665246, 0.38182704,
                          0.97697037, 0.89624951, 0.83816689, 0.39032672, 0.68852691,
                          0.85427299, 0.66773948, 0.9883756 , 0.06231242, 0.87223773,
                          0.12315628, 0.6709966 , 0.69840404, 0.19659599, 0.89613321,
                          0.79136648, 0.223493  , 0.8223016 , 0.81438973, 0.89321441,
                          0.82446078, 0.29317171, 0.23707863, 0.61798678, 0.85119219,
                          0.48374624, 0.18401164, 0.91021346, 0.18338856, 0.33581014,
                          0.51203974, 0.27373114, 0.26494583, 0.3909069 , 0.69500561,
                          0.44669268, 0.57792494, 0.92083106, 0.3760086 , 0.92086847,
                          0.98019746, 0.00930665, 0.91236066, 0.10589444, 0.98849831,
                          0.9937718 , 0.88560407, 0.20240646, 0.53980435, 0.21197806,
                          0.73511026, 0.91531094, 0.24419261, 0.79892127, 0.76623351,
                          0.32237578, 0.74371048, 0.89283081, 0.99471695, 0.59218713,
                          0.14229807, 0.26866107, 0.61418273, 0.53238885, 0.2847934 ,
                          0.33879263, 0.15419413, 0.81211755, 0.55982182, 0.33033445,
                          0.98925566, 0.21407401, 0.75933437, 0.50981508, 0.84659468,
                          0.27123332, 0.30602554, 0.78974943, 0.15961765, 0.75269879,
                          0.88404004, 0.25359787, 0.67575388, 0.10753205, 0.52492257,
                          0.2276367 , 0.57348205, 0.55631533, 0.48828726, 0.80950892,
                          0.68959411, 0.06038109, 0.3730253 , 0.44658293, 0.12323353,
                          0.90588169, 0.13484593, 0.58743073, 0.60592698, 0.67315081,
                          0.59887062, 0.3524358 , 0.47446065, 0.98078295, 0.31889862,
                          0.88225427, 0.81911728, 0.53942069, 0.6742203 , 0.73162166,
                          0.34597118, 0.70844054, 0.25029322, 0.29910746, 0.35906746,
                          0.53989701, 0.36776386, 0.04466711, 0.09399784, 0.53547235,
                          0.64757585, 0.03797524, 0.66378485, 0.34908186, 0.79601534,
                          0.10335962, 0.30468185, 0.57992791, 0.97139889, 0.40799129,
                          0.72985792, 0.65705408, 0.48913045, 0.46000752, 0.99260624,
                          0.52711571, 0.9383317 , 0.87126459, 0.16266698, 0.17428769,
                          0.11933665, 0.15703581, 0.17907467, 0.32411207, 0.56666047,
                          0.80868794, 0.49001672, 0.77590492, 0.63564239, 0.92169564,
                          0.5098178 , 0.40486284, 0.42978249, 0.4228141 , 0.33155423,
                          0.92695095, 0.40583767, 0.6635788 , 0.93854154, 0.48653097,
                          0.9305096 , 0.96009097, 0.03602708, 0.16874548, 0.4733966 ,
                          0.42363966, 0.18261131, 0.42653311, 0.48740795, 0.40008523,
                          0.35099519, 0.9641026 , 0.93447868, 0.16069199, 0.63925304,
                          0.17770585, 0.53600886, 0.72037108, 0.8454436 , 0.182311  ,
                          0.97860041, 0.41959913, 0.45109368, 0.24313804, 0.17884554,
                          0.04705525, 0.83247529, 0.89877392, 0.57362423, 0.27708354,
                          0.93649503, 0.43493419, 0.5422893 , 0.85565473, 0.86814896,
                          0.75788182, 0.02102082, 0.62643473, 0.24955471, 0.12775553,
                          0.96546452, 0.11001835, 0.82845919, 0.84811548, 0.21516077,
                          0.88871084, 0.55331041, 0.99744447, 0.70181741, 0.09492537,
                          0.18130881, 0.45487527, 0.82986703, 0.31207231, 0.08682494,
                          0.90971212, 0.40231716, 0.95428082, 0.10105085, 0.7243062 ,
                          0.87386255, 0.28549753, 0.90084605, 0.91034781, 0.44687279,
                          0.99318239, 0.58953929, 0.73074143, 0.05055378 };

    variableVector answer = { -6.34180841, -8.44235442, -7.66602685, -6.62791667, -6.30156652,
                              -7.81093903, -9.08319118, -7.62283755, -8.7120047 , -7.96533995,
                              -7.29110914, -7.63480242, -6.0360827 , -6.66816385, -6.38308499,
                              -8.06776472, -7.29777722, -7.77952498, -7.6470537 , -9.94159411,
                              -7.65257834, -5.90193479, -6.5591572 , -8.12839975, -8.56024681,
                              -7.40823637, -8.875604 };

    variableVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::computeReferenceHigherOrderStress( Gamma, C, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test the Jacobian

    variableVector resultJ;
    variableMatrix dMdGamma;
    error = tardigradeHydra::micromorphicLinearElasticity::computeReferenceHigherOrderStress( Gamma, C, resultJ, dMdGamma );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test dMdGamma
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < Gamma.size(); i++ ){
        constantVector delta( Gamma.size(), 0 );
        delta[i] = eps * fabs( delta[i] ) + eps;

        variableVector result_P, result_M;

        error = tardigradeHydra::micromorphicLinearElasticity::computeReferenceHigherOrderStress( Gamma + delta, C, result_P );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::computeReferenceHigherOrderStress( Gamma - delta, C, result_M );

        BOOST_CHECK( !error );

        constantVector gradCol = ( result_P - result_M ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dMdGamma[j][i] ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testComputeLinearElasticTerm3 ){
    /*!
     * Test the computation of the third term for micromorphic linear elasticity.
     *
     */

    variableVector invCGamma = { 3.50845589, -1.63297374,  1.45320064,  5.96175258,  3.87979583,
                                 1.87176717, -0.9248906 ,  2.62094337,  1.33218361,  0.32076702,
                                 1.4736958 ,  1.12865617, -0.64009695,  0.59763605,  1.30927157,
                                 1.3281096 , -0.5869566 ,  0.44018241,  0.57800689, -0.28019742,
                                 0.10974229,  0.84147417,  0.33054528,  0.22469373, -0.15308855,
                                 0.50800348,  0.14458868 };

    variableVector referenceHigherOrderStress = { -6.34180841, -8.44235442, -7.66602685, -6.62791667, -6.30156652,
                                                  -7.81093903, -9.08319118, -7.62283755, -8.7120047 , -7.96533995,
                                                  -7.29110914, -7.63480242, -6.0360827 , -6.66816385, -6.38308499,
                                                  -8.06776472, -7.29777722, -7.77952498, -7.6470537 , -9.94159411,
                                                  -7.65257834, -5.90193479, -6.5591572 , -8.12839975, -8.56024681,
                                                  -7.40823637, -8.875604 };

    variableVector answer = { -121.37119446,  -44.30225806,  -15.29816588,
                              -122.96816543,  -40.25141813,  -15.71355416,
                              -120.88746809,  -47.45291642,  -15.17132276 };

    variableVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress,
                                                                                               result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test the Jacobians

    variableVector resultJ;
    variableMatrix dTerm3dInvCGamma, dTerm3dM;

    error = tardigradeHydra::micromorphicLinearElasticity::computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress,
                                                                                      resultJ, dTerm3dInvCGamma, dTerm3dM );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test dTerm3dInvCGamma
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < invCGamma.size(); i++ ){
        constantVector delta( invCGamma.size(), 0 );
        delta[i] = eps * fabs( invCGamma[i] ) + eps;

        variableVector result_P, result_M;

        error = tardigradeHydra::micromorphicLinearElasticity::computeLinearElasticTerm3( invCGamma + delta, referenceHigherOrderStress,
                                                                                          result_P );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::computeLinearElasticTerm3( invCGamma - delta, referenceHigherOrderStress,
                                                                                         result_M );

        BOOST_CHECK( !error );

        constantVector gradCol = ( result_P - result_M ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dTerm3dInvCGamma[j][i] ) );
        }
    }

    //Test dTerm3dM
    for ( unsigned int i = 0; i < referenceHigherOrderStress.size(); i++ ){
        constantVector delta( referenceHigherOrderStress.size(), 0 );
        delta[i] = eps * fabs( referenceHigherOrderStress[i] ) + eps;

        variableVector result_P, result_M;

        error = tardigradeHydra::micromorphicLinearElasticity::computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress + delta,
                                                                                          result_P );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress - delta,
                                                                                          result_M );

        BOOST_CHECK( !error );

        constantVector gradCol = ( result_P - result_M ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dTerm3dM[j][i] ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testComputeInvRCGPsi ){
    /*!
     * Test the computation of the invRCG Psi product
     *
     */

    variableVector RCG = { 0.34852835, 0.47540122, 1.11252634,
                           0.47540122, 1.49184663, 1.57435946,
                           1.11252634, 1.57435946, 3.68235756 };

    variableVector Psi = { -0.02303102, -0.41101265, -0.36040573,
                           -0.14715403, -0.18045474, -0.8814645 ,
                           -0.36637526, -1.08887072, -1.44707636 };

    variableVector invRCG = tardigradeVectorTools::inverse( RCG, 3, 3 );

    variableVector answer = { 7.06496448, -6.60478112,  6.18226067,
                              0.01374041,  0.34618158, -0.31907041,
                             -2.23986034,  1.55175262, -2.12436531 };

    variableVector result;
    errorOut error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGPsi( invRCG, Psi, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobians

    variableVector resultJ;
    variableMatrix dInvRCGPsidRCG, dInvRCGPsidPsi;

    error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGPsi( invRCG, Psi, resultJ, dInvRCGPsidRCG, dInvRCGPsidPsi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    // Test dInvRCGPsidRCG
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < RCG.size(); i++ ){
        constantVector delta( RCG.size(), 0 );
        delta[i] = eps * fabs( delta[i] ) + eps;

        variableVector result_P, result_M;

        variableVector invRCG_P = tardigradeVectorTools::inverse( RCG + delta, 3, 3 );

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGPsi( invRCG_P, Psi, result_P );

        BOOST_CHECK( !error );

        variableVector invRCG_M = tardigradeVectorTools::inverse( RCG - delta, 3, 3 );

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGPsi( invRCG_M, Psi, result_M );

        BOOST_CHECK( !error );

        constantVector gradCol = ( result_P - result_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){ 
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dInvRCGPsidRCG[j][i] ) );
        }
    }

    // Test dInvRCGPsidPsi
    for ( unsigned int i = 0; i < Psi.size(); i++ ){
        constantVector delta( Psi.size(), 0 );
        delta[i] = eps * fabs( delta[i] ) + eps;

        variableVector result_P, result_M;

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGPsi( invRCG, Psi + delta, result_P );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGPsi( invRCG, Psi - delta, result_M );

        BOOST_CHECK( !error );

        constantVector gradCol = ( result_P - result_M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dInvRCGPsidPsi[j][i] ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testComputeInvRCGGamma ){
    /*!
     * Test the computation of the invRCG Gamma product
     *
     */

    variableVector RCG = { 0.34852835, 0.47540122, 1.11252634,
                           0.47540122, 1.49184663, 1.57435946,
                           1.11252634, 1.57435946, 3.68235756 };

    variableVector Gamma = { -0.31120922, -0.3563267 , -0.36573233, -0.0771914 , -0.24252804,
                             -0.4738459 , -0.35937075, -0.01781817, -0.17465609, -0.82225557,
                             -0.36719542, -0.86494826, -0.92750732, -1.18214541, -1.00423785,
                             -0.43125133, -0.19543115, -0.49736256, -0.67098335, -0.98433811,
                             -1.0183107 , -0.12645195, -0.79818076, -1.23318541, -0.95577138,
                              0.1274431 , -0.47648617 };

    variableVector invRCG = tardigradeVectorTools::inverse( RCG, 3, 3 );


    variableVector answer = { -8.75662978, -4.74843358, -4.6911341 , -3.16355851, -0.13179485,
                              -8.17350013, -5.69120721, -4.54526807, -2.48129262, -0.65657733,
                               0.06407926, -0.52611055, -1.06777508, -1.02709034, -0.58509241,
                              -0.02936234, -0.30663195, -0.35940994,  2.74408079,  1.13990438,
                               1.36569754,  1.37796288,  0.26218363,  2.38466645,  1.47244621,
                               1.53893867,  0.77392204 };

    variableVector result;
    errorOut error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGGamma( invRCG, Gamma, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test Jacobians

    variableVector resultJ;
    variableMatrix dInvRCGGammadRCG, dInvRCGGammadGamma;

    error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGGamma( invRCG, Gamma, resultJ, dInvRCGGammadRCG, dInvRCGGammadGamma );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    // Test dInvRCGGammadRCG
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < RCG.size(); i++ ){
        constantVector delta( RCG.size(), 0 );
        delta[i] = eps * fabs( delta[i] ) + eps;

        variableVector result_P, result_M;

        variableVector invRCG_P = tardigradeVectorTools::inverse( RCG + delta, 3, 3 );

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGGamma( invRCG_P, Gamma, result_P );

        BOOST_CHECK( !error );

        variableVector invRCG_M = tardigradeVectorTools::inverse( RCG - delta, 3, 3 );

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGGamma( invRCG_M, Gamma, result_M );

        BOOST_CHECK( !error );

        constantVector gradCol = ( result_P - result_M ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dInvRCGGammadRCG[j][i] ) );
        }
    }

    // Test dInvRCGGammadPsi
    for ( unsigned int i = 0; i < Gamma.size(); i++ ){
        constantVector delta( Gamma.size(), 0 );
        delta[i] = eps * fabs( delta[i] ) + eps;

        variableVector result_P, result_M;

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGGamma( invRCG, Gamma + delta, result_P );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::computeInvRCGGamma( invRCG, Gamma - delta, result_M );

        BOOST_CHECK( !error );

        constantVector gradCol = ( result_P - result_M ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dInvRCGGammadGamma[j][i] ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testFormIsotropicA ){
    /*!
     * Test the formation of the isotropic A stiffness tensor.
     *
     */

    parameterType lambda = 4;
    parameterType mu = 7;

    parameterVector answer = { 18.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,  7.,  0.,  7.,
                                0.,  0.,  0.,  0.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  7.,  0.,
                                0.,  0.,  7.,  0.,  7.,  0.,  0.,  0.,  0.,  0.,  4.,  0.,  0.,
                                0., 18.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  0.,  0.,  7.,  0.,
                                7.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.,
                                0.,  0.,  0.,  7.,  0.,  7.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,
                                0.,  0., 18. };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicA( lambda, mu, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testFormIsotropicB ){
    /*!
     * Test the formation of the isotropic B stiffness tensor.
     *
     */

    parameterType eta = 3;
    parameterType tau = 5;
    parameterType kappa = 6;
    parameterType nu = 8;
    parameterType sigma = 4;

    parameterVector answer = {  4.,  0.,  0.,  0., -2.,  0.,  0.,  0., -2.,  0.,  2.,  0.,  4.,
                                0.,  0.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  4.,  0.,
                                0.,  0.,  4.,  0.,  2.,  0.,  0.,  0.,  0.,  0., -2.,  0.,  0.,
                                0.,  4.,  0.,  0.,  0., -2.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,
                                4.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  0.,
                                0.,  0.,  0.,  4.,  0.,  2.,  0., -2.,  0.,  0.,  0., -2.,  0.,
                                0.,  0.,  4. };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicB( eta, tau, kappa, nu, sigma, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testFormIsotropicC ){
    /*!
     * Test the formation of the isotropic C stiffness tensor.
     *
     */

    parameterVector taus = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    parameterVector answer = {
        97.,  0.,  0.,  0., 13.,  0.,  0.,  0., 13.,  0., 16.,  0.,  9.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  9.,  0.,
        0.,  0., 23.,  0., 22.,  0.,  0.,  0.,  0.,  0., 23.,  0.,  0.,
        0.,  9.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,
        4.,  0.,  0.,  0., 23.,  0.,  0.,  0., 22.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  4.,  0.,  3.,  0., 23.,  0.,  0.,  0.,  2.,  0.,
        0.,  0.,  9.,  0., 22.,  0., 27.,  0.,  0.,  0.,  0.,  0., 26.,
        0.,  0.,  0., 16.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,
        7.,  0.,  3.,  0., 13.,  0.,  0.,  0., 23.,  0.,  0.,  0.,  5.,
        0., 26.,  0., 23.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,
        0.,  0.,  2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  8.,  0., 10.,
        0.,  0.,  0., 11.,  0.,  0.,  0.,  9.,  0.,  0.,  0.,  9.,  0.,
       12.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 22.,  0.,  0.,  0., 27.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,  7.,  0., 26.,  0.,
        0.,  0.,  6.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  0.,  0., 10.,
        0.,  8.,  0.,  0.,  0.,  9.,  0.,  0.,  0., 12.,  0.,  0.,  0.,
       11.,  0.,  9.,  0.,  0.,  0.,  0.,  0., 13.,  0.,  0.,  0.,  5.,
        0.,  0.,  0., 23.,  0.,  6.,  0.,  2.,  0.,  0.,  0.,  0.,  0.,
        0.,  0., 26.,  0.,  0.,  0., 23.,  0.,  0.,  0., 23.,  0., 26.,
        0.,  0.,  0.,  0.,  0., 23.,  0.,  0.,  0., 13.,  0.,  0.,  0.,
        5.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,  2.,  0., 16.,  0.,  0.,
        0., 26.,  0.,  0.,  0.,  6.,  0., 27.,  0., 22.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,
        0.,  0.,  0., 11.,  0.,  9.,  0.,  0.,  0.,  8.,  0.,  0.,  0.,
       10.,  0.,  0.,  0., 12.,  0.,  9.,  0.,  0.,  0.,  0.,  0.,  9.,
        0.,  0.,  0., 23.,  0.,  0.,  0.,  2.,  0., 22.,  0., 23.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  4.,  0.,  0.,
        0.,  9.,  0., 16.,  0.,  0.,  0.,  0.,  0., 13.,  0.,  0.,  0.,
       97.,  0.,  0.,  0., 13.,  0.,  0.,  0.,  0.,  0., 16.,  0.,  9.,
        0.,  0.,  0.,  4.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,
        0.,  0., 23.,  0., 22.,  0.,  2.,  0.,  0.,  0., 23.,  0.,  0.,
        0.,  9.,  0.,  0.,  0.,  0.,  0.,  9.,  0., 12.,  0.,  0.,  0.,
       10.,  0.,  0.,  0.,  8.,  0.,  0.,  0.,  9.,  0., 11.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  7.,  0.,  0.,  0.,
        0.,  0.,  0.,  0., 22.,  0., 27.,  0.,  6.,  0.,  0.,  0., 26.,
        0.,  0.,  0., 16.,  0.,  2.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,
        5.,  0.,  0.,  0., 13.,  0.,  0.,  0., 23.,  0.,  0.,  0.,  0.,
        0., 26.,  0., 23.,  0.,  0.,  0., 23.,  0.,  0.,  0., 26.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,  6.,  0., 23.,  0.,  0.,
        0.,  5.,  0.,  0.,  0., 13.,  0.,  0.,  0.,  0.,  0.,  9.,  0.,
       11.,  0.,  0.,  0., 12.,  0.,  0.,  0.,  9.,  0.,  0.,  0.,  8.,
        0., 10.,  0.,  0.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  6.,  0.,
        0.,  0., 26.,  0.,  7.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  0.,
        0., 27.,  0.,  0.,  0., 22.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
       12.,  0.,  9.,  0.,  0.,  0.,  9.,  0.,  0.,  0., 11.,  0.,  0.,
        0., 10.,  0.,  8.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,
        0.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 23.,  0., 26.,
        0.,  5.,  0.,  0.,  0., 23.,  0.,  0.,  0., 13.,  0.,  3.,  0.,
        7.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,  0.,  0., 16.,  0.,  0.,
        0., 26.,  0.,  0.,  0.,  0.,  0., 27.,  0., 22.,  0.,  9.,  0.,
        0.,  0.,  2.,  0.,  0.,  0., 23.,  0.,  3.,  0.,  4.,  0.,  0.,
        0.,  0.,  0.,  0.,  0., 22.,  0.,  0.,  0., 23.,  0.,  0.,  0.,
        4.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  9.,
        0.,  0.,  0., 23.,  0.,  0.,  0.,  0.,  0., 22.,  0., 23.,  0.,
        0.,  0.,  9.,  0.,  0.,  0., 16.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  9.,  0., 16.,  0., 13.,  0.,  0.,  0., 13.,  0.,  0.,  0.,
       97.
    };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicC( taus, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testFormIsotropicD ){
    /*!
     * Test the formation of the isotropic D stiffness tensor.
     *
     */

    parameterType tau = 5;
    parameterType sigma = 4;

    parameterVector answer = { 13.,  0.,  0.,  0.,  5.,  0.,  0.,  0.,  5.,  0.,  4.,  0.,  4.,
                                0.,  0.,  0.,  0.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,
                                0.,  0.,  4.,  0.,  4.,  0.,  0.,  0.,  0.,  0.,  5.,  0.,  0.,
                                0., 13.,  0.,  0.,  0.,  5.,  0.,  0.,  0.,  0.,  0.,  4.,  0.,
                                4.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  4.,  0.,  0.,  0.,  0.,
                                0.,  0.,  0.,  4.,  0.,  4.,  0.,  5.,  0.,  0.,  0.,  5.,  0.,
                                0.,  0., 13. };

    parameterVector result;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicD( tau, sigma, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testAssembleFundamentalDeformationMeasures ){
    /*!
     * Assemble the fundamental deformation measures from the degrees of freedom.
     *
     */

    double grad_u[ 3 ][ 3 ] = { { 1, 2, 3 },
                                { 4, 5, 6 },
                                { 7, 8, 9 } };

    double phi[ 9 ] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    double grad_phi[ 9 ][ 3 ] = { {  1,  2,  3 },
                                  {  4,  5,  6 },
                                  {  7,  8,  9 },
                                  { 10, 11, 12 },
                                  { 13, 14, 15 },
                                  { 16, 17, 18 },
                                  { 19, 20, 21 },
                                  { 22, 23, 24 },
                                  { 25, 26, 27 } };

    variableVector answerDeformationGradient = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerMicroDeformation = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerGradientMicroDeformation = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                     10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                     19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector resultF, resultChi, resultGradChi;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                                                            resultF, resultChi, resultGradChi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultF, answerDeformationGradient ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultChi, answerMicroDeformation ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultGradChi, answerGradientMicroDeformation ) );

    //Test the Jacobians
    variableVector resultFJ, resultChiJ, resultGradChiJ;
    variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

    error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                                                   resultFJ, resultChiJ, resultGradChiJ,
                                                                                                   dFdGradU, dChidPhi, dGradChidGradPhi );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultFJ, answerDeformationGradient ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultChiJ, answerMicroDeformation ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultGradChiJ, answerGradientMicroDeformation ) );

    //Test the jacobians w.r.t. the gradient of the displacement
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta( 3, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_u[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 3 ][ 3 ] =
        {
                { grad_u[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] + delta[ 2 ][ 2 ] }
        };

        double negative_perturb[ 3 ][ 3 ] =
        {
                { grad_u[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] - delta[ 2 ][ 2 ] }
        };

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( positive_perturb, phi, grad_phi,
                                                                                                       FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( negative_perturb, phi, grad_phi,
                                                                                                       FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dFdGradU[ j ][ i ] ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 9; i++ ){
        constantVector delta( 9, 0 );

        delta[ i ] = eps * fabs( phi[ i ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ] = { phi[ 0 ] + delta[ 0 ], phi[ 1 ] + delta[ 1 ], phi[ 2 ] + delta[ 2 ],
                                         phi[ 3 ] + delta[ 3 ], phi[ 4 ] + delta[ 4 ], phi[ 5 ] + delta[ 5 ],
                                         phi[ 6 ] + delta[ 6 ], phi[ 7 ] + delta[ 7 ], phi[ 8 ] + delta[ 8 ] };

        double negative_perturb[ 9 ] = { phi[ 0 ] - delta[ 0 ], phi[ 1 ] - delta[ 1 ], phi[ 2 ] - delta[ 2 ],
                                         phi[ 3 ] - delta[ 3 ], phi[ 4 ] - delta[ 4 ], phi[ 5 ] - delta[ 5 ],
                                         phi[ 6 ] - delta[ 6 ], phi[ 7 ] - delta[ 7 ], phi[ 8 ] - delta[ 8 ] };

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, positive_perturb, grad_phi,
                                                                                                       FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, negative_perturb, grad_phi,
                                                                                                       FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dChidPhi[ j ][ i ] ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_u[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ][ 3 ] =
        {
                { grad_phi[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] + delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] + delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] + delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] + delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] + delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] + delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] + delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] + delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] + delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] + delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] + delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] + delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] + delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] + delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] + delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] + delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] + delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] + delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] + delta[ 8 ][ 2 ] }
        };

        double negative_perturb[ 9 ][ 3 ] =
        {
                { grad_phi[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] - delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] - delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] - delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] - delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] - delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] - delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] - delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] - delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] - delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] - delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] - delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] - delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] - delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] - delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] - delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] - delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] - delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] - delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] - delta[ 8 ][ 2 ] }
        };

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, positive_perturb,
                                                                                      FP, chiP, gradChiP );

        BOOST_CHECK( !error );

        error = tardigradeHydra::micromorphicLinearElasticity::assembleFundamentalDeformationMeasures( grad_u, phi, negative_perturb,
                                                                                                       FM, chiM, gradChiM );

        BOOST_CHECK( !error );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dGradChidGradPhi[ j ][ i ] ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testExtractMaterialParameters ){
    /*!
     * Test the extraction of the material parameters.
     *
     */

    std::vector< double > fparams = { 2, 1.7, 1.8,
                                      5, 2.8, .76, .15, 9.8, 5.4,
                                      11, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.,
                                      2, .76, 5.4 };

    parameterVector Amatrix;
    parameterVector Bmatrix;
    parameterVector Cmatrix;
    parameterVector Dmatrix;

    parameterVector answerAmatrix;
    parameterVector answerBmatrix;
    parameterVector answerCmatrix;
    parameterVector answerDmatrix;

    errorOut error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicA( 1.7, 1.8, answerAmatrix );
    BOOST_CHECK( !error );

    error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicB( 2.8, 0.76, 0.15, 9.8, 5.4, answerBmatrix );
    BOOST_CHECK( !error );

    error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicC( { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.}, answerCmatrix );
    BOOST_CHECK( !error );

    error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicD( 0.76, 5.4, answerDmatrix );
    BOOST_CHECK( !error );

    error = tardigradeHydra::micromorphicLinearElasticity::extractMaterialParameters( fparams, Amatrix, Bmatrix, Cmatrix, Dmatrix );

    BOOST_CHECK( !error );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( Amatrix, answerAmatrix ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( Bmatrix, answerBmatrix ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( Cmatrix, answerCmatrix ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( Dmatrix, answerDmatrix ) );
}
