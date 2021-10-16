#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include "BitUtility.h"
#include "StateGenerator.h"
#include "ScriptGenerator.h"
#include "Spectrum.h"
#include "NormCalculator.h"
#include "otoc.h"

using namespace std;

struct Options {
    int L = -1;
    int M = -1;
    f_type N = std::numeric_limits<f_type>::infinity();
    f_type beta = std::numeric_limits<f_type>::infinity();
    f_type bt = std::numeric_limits<f_type>::infinity();
    f_type tMin;
    f_type tMax;
    f_type tSteps;
    std::string W;
    std::string V;
    std::string Op;
    char ComputeNorm = 0;
    char ExcludeZeroStates = 0;
    char Parts = 1;
    f_type alpha = 1.0;
    f_type sigma = 0.0;
};

void TestBlockMatrixExp() {
    IOperator* op = IOperator::Create("XXXZ");
    int L = 11;
    int N = 13;
    BlockMatrix res(L);
    op->ToMatrixN(L, N, res);
    Stopwatch watch;
    BlockMatrix exp1 = res.exp();
    std::cout << "elapsed: " << watch.Elapsed() << std::endl;

    Matrix mat = res.matrix();
    watch.Start();
    Matrix exp2 = mat.exp();
    std::cout << "elapsed: " << watch.Elapsed() << std::endl;

    std::cout << "norm=" << exp2.norm() << std::endl;
    std::cout << "diff=" << (exp1.matrix() - exp2).norm() << std::endl;

    delete op;
}

void TestMatrixExp() {
    DilatationOperator ham;
    int L = 13;
    int M = 6;
    f_type N = 10;
    f_type beta = 0.9;
    int size = StateCollection::Inst()->StateNumber(L, M);
    Matrix res(size, size);
    ham.ToMatrixN(L, M, N, beta, res);
    
    Stopwatch watch;
/*    watch.Start();
    std::vector<var_t> v1;
    for (f_type t = 0.0; t <= 10; t += 0.5) {
        std::cout << "t=" << t << " ";
        var_t z(0, -t);
        Matrix exp1 = (res * z).exp();
        v1.push_back(exp1.trace());
    }

    std::cout << endl;
    std::cout << "v1=" << v1 << std::endl;

    std::cout << "time used:" << watch.Stop() << " seconds." << std::endl << std::endl;;*/

    watch.Start();
    Eigen::ComplexEigenSolver<Matrix> solver(res, true);
    Eigen::ArrayXcd diag = solver.eigenvalues().array();
//    Vector diag = solver.eigenvalues();
    Matrix V = solver.eigenvectors();
    Matrix Vinv = V.inverse();

//    std::cout << "diag=" << diag << std::endl;
    std::vector<var_t> v2;
    std::cout << "initial time= " << watch.Elapsed() << std::endl;

    for (f_type t = 0.0; t <= 10; t += 0.5) {
        std::cout << "t=" << t << " ";
        var_t z(0, -t);
        Matrix exp2 = V * (diag * z).exp().matrix().asDiagonal() * Vinv;
//        Matrix exp2 = V * (diag.array() * z).exp().matrix().asDiagonal() * Vinv;
//        std::cout << "(diag*z).exp()=" << (diag * z).exp() << std::endl;
        v2.push_back(exp2.trace());
    }
    std::cout << std::endl;
    std::cout << "v2=" << v2 << std::endl;

    std::cout << "time used:" << watch.Stop() << " seconds." << std::endl;
}

void TestOtoc() {
    int L = 4;
    int M = -1;
    f_type N = 17.0;
    f_type beta = 0.9;
    f_type b = 0.5;
    OtocParameters params = {"XZ", "ZX", L, M, N, beta, b, 1};
    std::cout << "test otoc" << std::endl;
    IOtoc* otoc = IOtoc::Create(params);
    otoc->Init();
    std::cout << "finished initialization..." << std::endl;
    for (int i = 0; i < 10; i++) {
        std::cout << i << ' ' << otoc->Compute((f_type)i) << std::endl;
    }
}

void TestForAnything() {
    Matrix m;
    std::cout << m.rows() << " " << m.cols() << std::endl;
}

void TestOp() {
    FourFieldOperator ffp("XZXZ");
    TraceState ts;
    ts.AddTrace(SingleTrace("ZXZX"));
    MixState ms;
    ffp.ApplyOn(ts, ms);
    std::cout << ms << std::endl;
}

void TestNormCalculator() {
    StateGenerator generator;
    generator.GenerateAllStates();
    generator.InitStateCollection(StateCollection::Inst());

    SingleTrace trace("ZZXX");
    TraceState ts;
    ts.AddTrace(trace);
    BruteForceCalculator calc;
    PowerSeries<i64> res = calc.Calculate(ts, ts);
    for (int i = 0; i <= res.HighestOrder(); i++) {
        std::cout << res.Coef(i) << " ";
    }
    std::cout << std::endl;
}

void TestNormCalculator2() {
    StateGenerator generator;
    generator.GenerateAllStates();
    generator.InitStateCollection(StateCollection::Inst());
    std::cout << "init done" << std::endl;
    /*    StateId id(0, 0);
    std::cout << StateCollection::Inst()->GetState(id) << std::endl;*/

    SingleTrace trace("ZZ");
    SingleTrace trace2("ZZXX");
    SingleTrace trace3("ZZXZZX");
    trace.Normalize();
    trace2.Normalize();
    trace3.Normalize();
    TraceState ts;
    ts.AddTrace(trace);
    ts.AddTrace(trace2);
    ts.QuickNormalize();
    TraceState ts2;
    ts2.AddTrace(trace3);
    RecursiveNormCalculator calc;
    PowerSeries<nType> res = calc.Calculate(ts, ts2);
    std::cout << "<" << trace << "^" << trace << ">=";
    std::cout << res << std::endl;

    BruteForceCalculator calc2;
    PowerSeries<nType> res2 = calc2.Calculate(ts, ts2);
    std::cout << "res=" << res2 << std::endl;
}

void ComputeOtoc(OtocParameters params, f_type tMin, f_type tMax, int steps, bool parts = false) {
    string rootFolder = "";
    ScriptGenerator gen(rootFolder);
    gen.SaveOtoc(params, tMin, tMax, steps, parts);
}

void GenerateHamMatrix(int L, int M)
{
    string rootFolder = "";
    ScriptGenerator gen(rootFolder);
    if (M >= 0 && M <= L) {
        gen.SaveHamMatrix(L, M);
        gen.SaveStates(L, M);
    }
    else {
        for (int i = 0; i <= L; i++) {
            gen.SaveHamMatrix(L, i);
            gen.SaveStates(L, i);
        }
    }
}

void GenerateOpMatrix(std::string opStr, int L, int M = -1)
{
    string rootFolder = "";
    ScriptGenerator gen(rootFolder);
    if (M >= 0) {
        gen.SaveOpMatrix(opStr, L, M);
    } else {
        IOperator* op = dynamic_cast<ISimpleOperator*>(IOperator::Create(opStr));
        if (!op->MagnonFixed()) {
            gen.SaveOpMatrix(opStr, L);
        }
        else {
            for (int m = 0; m <= L; m++) {
                gen.SaveOpMatrix(opStr, L, m);
            }
        }
        delete op;
    }
}

void GenerateNormMatrix(int L, int M, bool su)
{
    string rootFolder = "";
    ScriptGenerator gen(rootFolder);
    if (M < 0) {
        if (su) {
            gen.SaveNormToDataFile2(L);
        }
        else {
            gen.SaveNormToDataFile(L);
        }
        return;
    }
    if (su) {
        gen.SaveNormToDataFile2(L, M);
    }
    else {
        gen.SaveNormToDataFile(L, M);
    }
}

void ComputeSpectrum(int L, int M, int N, f_type beta)
{
    string rootFolder = "";
    ScriptGenerator gen(rootFolder);
    if (M >= 0) {
        gen.SaveSpectrum(L, M, N, beta);
    }
    else {
        for (int i = 0; i <= L; i++) {
            gen.SaveSpectrum(L, i, N, beta);
        }
    }
}

std::vector<std::string> optionDesc {
    "-L <integer number>", "Length of spin chain, 2 <= L <= " + ToString(MAX_BIT_TO_GENERATE),
    "-M <integer number>", "Number of magnons.",
    "-N <real number>", "Rank of SU(N) color group. In general it's integer but could be real.",
    "-b <real number>", "The beta parameter of beta-deformed N=4 SYM.",
    "-bt <real number>", "The reciprocal of temperature.",
    "-tmin <real number>", "The initial time of OTOC.",
    "-tmax <real number>", "The final time of OTOC.",
    "-tstep <integer number>", "Number of steps go from initial time to final time.",
    "-o <string>", "The operator to generate matrix for. The values can be NTr, NTrX, ", 
    " ", "or any two-letter or four-letter string built out of X and Z.",
    " ", "NTr indicates the operator counts number of traces;",
    " ", "NTrX indicates the operator counts number of traces containing X;",
    " ", "XZ indicates the trace operator Tr(X d/dZ); XZXZ indicates the operator Tr(XZ d/dX d/dZ)",
    "-W <string>", "The W operator of OTOC. Its scope of values are the same as -o option.",
    "-V <string>", "The V operator of OTOC. Its scope of values are the same as -o option.",
    "-n", "The option only used in 'otoc' command. If not spefied, the norm matrix will be computed;",
    " ", "otherwise, the norm matrix will be loaded form local files generated by 'norm' command",
    "-nozero", "If the option is specified, zero energy states is excluded in C(t) computation.",
    "-alpha <real number>", "The alpha parameter of the OTOC regularization",
    "-sigma <real number>", "The sigma parameter of the OTOC regularization"
};

void DisplayOptions() {
    std::cout << "General Options:" << std::endl;
    for (size_t i = 0; i < optionDesc.size(); i += 2) {
        std::cout << "  " << std::left << std::setfill(' ') << std::setw(25) << optionDesc[i];
        std::cout << "\t\t" << optionDesc[i + 1] << std::endl;
    }
}

std::vector<std::string> exampleDesc {
    "Generate Hamiltonian matrix and trace states for L=4, M=2",                         // 1
    "Generate all Hamiltonian matrix and trace states for L=4",                          // 2
    "Generate norm matrix for L=4, M=2",                                //  3
    "Generate all norm matrix for L=4",                                 // 4
    "Compute energy spectrum for L=4, M=2, N=17, beta=0.9",             // 5
    "Generate matrix for operator counting traces containing X for L=4",
    "Compute OTOC for W=Tr(Xd/dZ), V=Tr(Zd/dX), L=4, N=5, beta=0.9, beta of temperature=0.5,\n \
\t\tand time run from 0.0 to 5.0 with step size 0.1. Note that to run the example one needs\n\
\t\tto run \"otoc norm - L 4\" first to generate all norm matrices for L = 4 and then copy the\n \
\t\tgenerated files to otocdata folder",
    "Compute OTOC for W=Tr(XZ d/dX d/dZ), V=Tr(ZX d/dX d/dZ), L=4, M=2, beta=0.9, beta of temperature=0.5,\n \
\t\tand N=Infinity. As -n option is specified, it will compute norm matrix rather than load it from local files.",
    "Compute C(t) for W=Tr(XX d/dX d/dZ), V=Tr(ZZ d/dZ d/dX), L=4, N=5. Zero energy states are excluded.",
    "Compute C(t) for W=Tr(XX d/dX d/dZ), V=Tr(ZZ d/dZ d/dX), L=4, N=5. With regularization parameters alpha=0.5 and sigma=0.25.",
    "Example of computing OTOC for composite operators. Currently, addition (+), substraction(-) , \n\
\t\tmultiplication(*), and exponentiate (exp), and normal ordered (~) are supported. No division (/) \n\
\t\tand space is allowd in the expression. Complex number should write in the form (real,imag).",
};

std::vector<std::string> examples {
    "otoc ham -L 4 -M 2",                       // 1
    "otoc ham -L 4",                            // 2
    "otoc norm -L 4 -M 2",                      // 3
    "otoc norm -L 4",                           // 4
    "otoc spec -L 4 -M 2 -N 17 -b 0.9",         // 5
    "otoc operator -o NTrX -L 4",                          // 6
    "otoc otoc -W XZ -V ZX -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50",
    "otoc otoc -W XZXZ -V XZXZ -L 4 -M 2 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50 -n",
    "otoc otoc -W XXXZ -V ZZZX -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50 -nozero",
    "otoc otoc -W XXXZ -V ZZZX -alpha 0.5 -sigma 0.25 -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50",
    "otoc otoc -W -0.2*exp(0.2*XXXZ+0.3*XZ)+(0.3,-0.5) -V ~(-XZXZ+0.5*ZZZX-0.2*ZX) -L 4 -N 5 -b 0.9 -bt 0.5 -tmin 0.0 -tmax 5.0 -tstep 50",
};

void DisplayExample() {
    std::cout << "Examples:" << std::endl;
    for (int i = 0; i < examples.size(); i++) {
        std::cout << "  Example " << i + 1 << ":\t" << exampleDesc[i] << std::endl;
        std::cout << "\t\t\t" << examples[i] << std::endl << std::endl;
    }
}

std::vector<std::string> commands{ "ham", "norm", "spec", "operator", "otoc" };

void DisplayHelp() {
    std::cout << "Usage:\n  otoc <command> [options]" << std::endl << std::endl;
    std::cout << "Commands:" << std::endl;

    std::vector<std::string> commandDesc {
        "Generate and save Hamiltonian matrix and trace states.",
        "Compute and save spectrum of Hamiltonian.",
        "Generate and save norm matrix of trace states.",
        "Generate and save matrix of operators.",
        "Compute and save out-of-time-order correlator."
    };

    
    for (int i = 0; i < commands.size(); i++) {
        std::cout << "  " << std::left << std::setfill(' ') << std::setw(25) << commands[i];
        std::cout << "\t\t" << commandDesc[i] << std::endl;
    }

    std::cout << std::endl;
    DisplayOptions();

    std::cout << std::endl;
    DisplayExample();
}

std::vector<string> Tokenize(std::string str, char seperator) {
    std::vector<std::string> ret;
    stringstream ss(str);
    std::string token;
    while (getline(ss, token, seperator)) {
        ret.push_back(token);
    }

    return ret;
}

void InitStates(int maxL) {
    StateGenerator generator(maxL);
    generator.GenerateAllStates();
    generator.InitStateCollection(StateCollection::Inst());
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        DisplayHelp();
        return -1;
    }

    std::string command = argv[1];
    if (command == "test") {
        InitStates(16);
//        TestForAnything();
//        TestOtoc();
        //TestOp();
//        TestMatrixExp();
        TestBlockMatrixExp();
        return 0;
    }

    if (std::find(commands.begin(), commands.end(), command) == commands.end()) {
        DisplayHelp();
        return -1;
    }

    Options options;
    bool displayHelp = false;

    for (int i = 2; i < argc; i++) {
        string cmd = argv[i];
        if (cmd == "-L") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.L = ToExpression(val, -1.0L);
        }
        else if (cmd == "-M") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.M = ToExpression(val, -1);
        }
        else if (cmd == "-N") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.N = ToExpression(val, -1.0);
        }
        else if (cmd == "-b") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.beta = ToExpression(val, 1000000.0);
            beta_nsym = options.beta;
        }
        else if (cmd == "-bt") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.bt = ToExpression(val, options.bt);
            inverseTempeture = options.bt;
        }
        else if (cmd == "-tmax") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.tMax = ToExpression(val, options.tMax);
        }
        else if (cmd == "-tmin") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.tMin = ToExpression(val, options.tMin);
        }
        else if (cmd == "-tstep") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            options.tSteps = ToExpression(val, options.tSteps);
        } else if (cmd == "-o") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            options.Op = argv[i];
        }
        else if (cmd == "-W") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            options.W = argv[i];
        }
        else if (cmd == "-V") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            options.V = argv[i];
        }
        else if (cmd == "-alpha") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            options.alpha = ToExpression(argv[i], options.alpha);
        }
        else if (cmd == "-sigma") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            options.sigma = ToExpression(argv[i], options.sigma);;
        }
        else if (cmd == "-n") {
            options.ComputeNorm = 1;
        }
        else if (cmd == "-nozero") {
            options.ExcludeZeroStates = 1;
        }
//        else if (cmd == "-parts") {
//            options.Parts = 1;
 //       }
        else {
            std::cerr << "invalid option " << cmd << std::endl;
            return -1;
        }
    }

    if (displayHelp) {
        DisplayHelp();
        return -1;
    }

    if (options.L > MAX_BIT_TO_GENERATE || options.L <= 1 || options.M > options.L || options.M < -1) {
        displayHelp = true;
    }

    if (options.N < 1) {
        displayHelp = true;
    }

    if (displayHelp) {
        DisplayHelp();
        return -1;
    }

    Stopwatch watch;
    watch.Start();
    LOG("Initializing trace states...", Verbos);
    InitStates(options.L);

    try {
        if (command == commands[0]) {
            GenerateHamMatrix(options.L, options.M);
        }
        else if (command == commands[1]) {
            GenerateNormMatrix(options.L, options.M, true);
        }
        else if (command == commands[2]) {
            ComputeSpectrum(options.L, options.M, options.N, options.beta);
        }
        else if (command == commands[3]) {
            GenerateOpMatrix(options.Op, options.L, options.M);
        }
        else if (command == commands[4]) {
            OtocParameters params = { options.W, options.V, options.L, options.M, options.N,
                options.beta, options.bt, options.ComputeNorm, options.ExcludeZeroStates,
                options.alpha, options.sigma };
            ComputeOtoc(params, options.tMin, options.tMax, options.tSteps, options.Parts);
        }

        LOG("The task is completed in " << watch.Stop() << " seconds.", Info);
    }
    catch (OtocException& ex){
        LOG(ex.what(), Error);
        LOG("Elapsed time: " << watch.Stop() << " seconds.", Info);
    }
}
