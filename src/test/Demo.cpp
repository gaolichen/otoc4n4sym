//Link to Boost
#define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../StateGenerator.h"
#include "../PowerSeries.h"
#include "../HamOperator.h"
#include "../ScriptGenerator.h"
#include "../Spectrum.h"


// test suite
BOOST_FIXTURE_TEST_SUITE(Demo_suite, SimpleTestFixture, * utf::label("Demo"))

BOOST_AUTO_TEST_CASE(GenerateStatesDemo, * utf::disabled())
{
    StateGenerator generator;
	for (int i = 1; i <= MAX_BIT_TO_COUNT; i++)
	{
	    cout << i << "\t" << generator.SingleStateNumber(i) << "\t" << generator.StateNumber(i) << endl;
	}

	generator.GenerateAllStates();
	generator.InitStateCollection(StateCollection::Inst());
	for (int n = 4; n < 5; n++) {
	    cout << "n=" << n << endl;
	    for (int i = 0; i < generator.SingleStateNumber(n); i++) {
	        cout << i << ' ' << generator.SingleTraceState(n, i) << endl;
	    }
	    cout << "=======" << endl;
	    for (int i = 0; i < generator.StateNumber(n); i++) {
	        cout << i << ' ' << generator.State(n, i) << endl;
	    }
	}
}

BOOST_AUTO_TEST_CASE(HamOperatorDemo, * utf::disabled())
{
    DilatationOperator ham;
    int L = 10;
    int M = 4;
    const TraceState& state = *StateCollection::Inst()->Begin(L, M);
//    TraceState state;
//    state.AddTrace(SingleTrace("ZX"));
//    state.AddTrace(SingleTrace("ZX"));
    MixState ms;
    cout << "state=" << state << endl;
    ham.ApplyOn(state, ms);
    cout << "ms=" << ms << endl;
}

BOOST_AUTO_TEST_CASE(HamMatrixDemo, * utf::disabled())
{
    DilatationOperator ham;
    int input[1][2] = {{4,2}};
    
    for (int n = 0; n < 1; n++) {
        int L = input[n][0];
        int M = input[n][1];
        vector<vector<NExpansionSeries> > mat;    
        ham.Matrix(L, M, mat);
        for (int i = 0; i < mat.size(); i++) {
            for (int j = 0; j < mat.size(); j++) {
                if (j) cout << '\t';
                cout << mat[i][j];
            }
            cout << endl;
        }
    }
}

BOOST_AUTO_TEST_CASE(ScriptGeneratorDemo, * utf::disabled())
{
    string rootFolder = "/home/gaolichen/gitroot/mywork/chaos/build";
    ScriptGenerator gen(rootFolder);
//    int input[9][2] = {{4,2}, {6,2}, {8,2}, {8,4}, {10,2}, {10,4}, {12,2}, {12, 4}, {12, 6}};
    int input[3][2] = {{16, 2}, {16, 3}, {16, 4}};
    
    for (int n = 0; n < 3; n++) {
        int L = input[n][0];
        int M = input[n][1];
        gen.SaveHamMatrix(L, M);
        gen.SaveStates(L, M);
    }
}

BOOST_AUTO_TEST_CASE(SingleTraceDemo, * utf::disabled())
{
    SingleTrace single("XXZZXZ");
    SingleTrace a, b, c;
    single.Split(3, 1, a, b, c);
//    single.Split(2, a, b);
    cout << "single=" << single << endl;
    cout << "a=" << a <<", b=" << b << ", c=" << c << endl;
}

BOOST_AUTO_TEST_CASE(SpectrumDemo, * utf::disabled())
{
    Spectrum spectrum;
    int N = 10;
    f_type beta = 0.9;
    int L = 4;
    int M = 2;
    std::vector<var_t> eigens;
    spectrum.Eigenvalues(L, M, 10, 0.9, eigens);
    std::cout << "eigens=" << eigens << std::endl;
    
    int zeros = 0;
    for (int i = 0; i < eigens.size(); i++) {
        if (std::abs(eigens[i]) < EPS) zeros++;
    }
    std::cout << "number of eigenvalues = " << eigens.size() << ", # of zero eigenvalues = " << zeros << std::endl;

}

BOOST_AUTO_TEST_SUITE_END()

