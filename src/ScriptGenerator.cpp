#pragma warning(disable:4018)
#include <fstream>
#include <iomanip>
#include "ScriptGenerator.h"
#include "Spectrum.h"
#include "NormCalculator.h"
using namespace std;

ScriptGenerator::ScriptGenerator(string rootFolder)
{
    this->rootFolder = rootFolder;
}

void ScriptGenerator::SaveMatrix(vector<vector<NExpansionSeries> >& mat, ofstream& ofs) {
	for (int i = 0; i < mat.size(); i++) {
		bool hasOutput = false;
		for (int j = 0; j < mat.size(); j++) {
			if (mat[i][j].IsZero()) continue;
			if (!hasOutput) {
				ofs << i + 1 << endl;
				hasOutput = true;
			}
			ofs << j + 1 << '\t' << mat[i][j] << endl;
		}
	}
}

void ScriptGenerator::SaveHamMatrix(int chainLength, int magnonNumber)
{
    string file = CombinePath(rootFolder, "ham_L" + ToString(chainLength) + "M" + ToString(magnonNumber) + ".txt");
	cout << "Generating hamiltonian matreix for L=" << chainLength << ", M=" << magnonNumber << endl;
	
	ofstream ofs(file.c_str());
	DilatationOperator ham;
    vector<vector<NExpansionSeries> > mat;
    ham.ToMatrix(chainLength, magnonNumber, mat);
	SaveMatrix(mat, ofs);
	ofs.close();
	cout << "Hamiltonian matrix is saved to file " << file << endl;
}

void ScriptGenerator::SaveOpMatrix(std::string op, int chainLength) {
	StateCollection* inst = StateCollection::Inst();
	int n = inst->StateNumber(chainLength);
	string file = CombinePath(rootFolder, "op" + op + "_L" + ToString(chainLength) + ".txt");
	ISimpleOperator* opr = ISimpleOperator::Create(op);

	ofstream ofs(file.c_str());
	ofs << n << std::endl;
	cout << "Generating matreix for operator " << op << " L=" << chainLength << endl;

	for (int i = 0; i < n; i++)
	{
		bool hasOutput = false;
		MixState ms;
		opr->ApplyOn(inst->GetState(chainLength, i), ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			if (!hasOutput) {
				ofs << i + 1 << std::endl;
				hasOutput = true;
			}
			ofs << it2->first.Index + 1 << '\t' << it2->second << std::endl;
		}
	}

	ofs.close();
	cout << "Op matrix is saved to file " << file << endl;
	delete opr;
}

void ScriptGenerator::SaveOpMatrix(std::string op, int L, int M) {
	ISimpleOperator* opr = ISimpleOperator::Create(op);

	StateCollection* inst = StateCollection::Inst();
	auto beginIt = inst->Begin(L, M);
	int beginIndex = inst->GetId(*beginIt).Index;
	int n = inst->StateNumber(L, M);

	string file = CombinePath(rootFolder, "op" + op + "_L" + ToString(L) + "M" + ToString(M) + ".txt");
	ofstream ofs(file.c_str());
	ofs << n << std::endl;
	cout << "Generating matreix for operator " << op << " L=" << L << ", M=" << M << endl;

	for (auto it = beginIt; it != inst->End(L, M); ++it)
	{
		bool hasOutput = false;
		MixState ms;
		opr->ApplyOn(*it, ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			if (!hasOutput) {
				ofs << it - beginIt + 1 << std::endl;
				hasOutput = true;
			}
			ofs << it2->first.Index - beginIndex + 1 << '\t' << it2->second << std::endl;
		}
	}

	ofs.close();
	cout << "Op matrix is saved to file " << file << endl;
	delete opr;
}

/*
void ScriptGenerator::SaveHamMatrix(int chainLength, int magnonNumber)
{
    vector<ofstream> ofs;
    for (int i = 0; i < 4; i++) {
        string file = CombinePath(rootFolder, "ham" + ToString(i) + "_L" + ToString(chainLength) + "M" + ToString(magnonNumber) + ".txt");
        ofs.push_back(ofstream(file.c_str()));
    }
	cout << "saving hamiltonian matreix for L=" << chainLength << ", M=" << magnonNumber << endl;
	
	DilatationOperator ham;
    vector<vector<NExpansionSeries> > mat;
    ham.Matrix(chainLength, magnonNumber, mat);
    for (int i = 0; i < mat.size(); i++) {
        vector<char> hasOutput{0, 0, 0, 0};
        
        for (int j = 0; j < mat.size(); j++) {
            if (mat[i][j].IsZero()) continue;
            for (int k = 0; k < ofs.size(); k++) {
                const ExpBetaSeries& ebs = mat[i][j].Coef(-k);
                if (ebs.IsZero()) continue;
                if (!hasOutput[k]) {
                    hasOutput[k] = 1;
                    ofs[k] << i + 1 << endl;
                }
                ofs[k] << j + 1 << '\t' << ebs << endl;
            }
        }
    }

    for (int k = 0; k < ofs.size(); k++) {
	    ofs[k].close();
	}
	cout << "Finished saving Hamiltonian matrix." << endl;
}*/

void ScriptGenerator::SaveStates(int L, int M)
{
    string file = CombinePath(rootFolder, "sta_L" + ToString(L) + "M" + ToString(M) + ".txt");
	cout << "Saving trace states for L=" << L << ", M=" << M << endl;
	ofstream ofs(file.c_str());

    StateCollection* inst = StateCollection::Inst();
    int i = 0;
    for (auto it = inst->Begin(L, M); it != inst->End(L, M); it++)
	{
	    ofs << ++i << '\t' << *it << endl;
	}
	cout << "Trace states are saved to file " << file << endl;
}

void ScriptGenerator::SaveSpectrum(int L, int M, int N, f_type beta)
{
    string file = CombinePath(rootFolder, "spec_L" + ToString(L) + "M" + ToString(M) + "N" + ToString(N) + "b" + ToString(beta) + ".txt");
	cout << "computing energy spectrum for L=" << L << ", M=" << M << endl;
	ofstream ofs(file.c_str());
	
	Spectrum spectrum;
    std::vector<var_t> eigens;
    spectrum.Eigenvalues(L, M, N, beta, eigens);
    
    int zeros = 0;
    for (int i = 0; i < eigens.size(); i++) {
        ofs << i + 1 << '\t' << std::setprecision(15) << Chop(eigens[i].real()) << endl;
        if (std::abs(eigens[i]) < EPS) zeros++;
    }
    std::cout << "number of eigenvalues = " << eigens.size() << ", # of zero eigenvalues = " << zeros << std::endl;
	cout << "Energy eigenvalues are saved to file " << file << endl;
}

void ScriptGenerator::SaveNormToDataFile(int bits, int fNumber) {
	BruteForceCalculator calc;
	SaveNormToDataFile(bits, fNumber,calc, 0);
}

void ScriptGenerator::SaveNormToDataFile(int bits) {
	BruteForceCalculator calc;
	for (int i = 0; i + i <= bits; i++) {
		SaveNormToDataFile(bits, i, calc, 1);
		if (i + i != bits) {
			SaveNormToDataFile(bits, bits - i, calc, 0);
		}
	}
}

void ScriptGenerator::SaveNormToDataFile2(int bits, int fNumber) {
//	RecursiveNormCalculator calc;
	DfsNormCalculator calc;
	SaveNormToDataFile(bits, fNumber, calc, 1);
}

void ScriptGenerator::SaveNormToDataFile2(int bits) {
//	RecursiveNormCalculator calc;
	DfsNormCalculator calc;
	for (int i = 0; i + i <= bits; i++) {
		SaveNormToDataFile(bits, i, calc, 1);
		if (i + i != bits) {
			SaveNormToDataFile(bits, bits - i, calc, 1);
		}
	}
}

void ScriptGenerator::SaveNormToDataFile(int bits, int fNumber, NormCalculator& calc, int format)
{
	StateCollection* sc = StateCollection::Inst();
	vector<int > stateIndex;
	int tot = sc->StateNumber(bits);

	for (int i = 0; i < tot; i++)
	{
		TraceState state = sc->GetState(bits, i);
		if (state.MagnonNumber() == fNumber)
		{
			stateIndex.push_back(i);
		}
	}

	if (stateIndex.size() == 0)
	{
		cout << "no states found for L= " << bits << " M=" << fNumber << endl;
		return;
	}

	string name = "norm_L" + ToString(bits) + "M" + ToString(fNumber);
	if (format == 1) {
		name += "_su";
	}
	string file = CombinePath(rootFolder, name + ".txt");
	cout << "computing norm matrix for L=" << bits << ", M=" << fNumber << endl;
	//cout << "file =" << file << endl;
	ofstream ofs(file.c_str());

	ofs << stateIndex.size() << endl << endl;

	for (int i = 0; i < stateIndex.size(); i++)
	{
		TraceState state = sc->GetState(bits, stateIndex[i]);
		for (int j = i; j < stateIndex.size(); j++)
		{
			PowerSeries<nType> poly = calc.Calculate(state, sc->GetState(bits, stateIndex[j]));
//			if (format == 1)
//			{
//				ofs << poly << std::endl;
//			}
//			else
//			{
				if (poly.IsZero())
				{
					ofs << 0 << endl;
				}
				else
				{
					int maxPow = poly.HighestOrder();
					int offset = (bits - maxPow) % 2;
					if (offset == 0) { ofs << 0; }
					else { ofs << 1; }

					int lowPow = 0;
					if (format == 1) {
						lowPow = -(bits - 2);
					}

					for (int k = bits - offset; k >= lowPow; k -= 2)
					{
						ofs << ' ' << poly.Coef(k);
					}
					ofs << endl;
				}
//			}
		}

		ofs << endl;
	}

	ofs.close();
	cout << "Norm data are saved to file " << file << endl;
}

std::string pow2String(std::string var, float_t pow) {
	if (pow == 0) {
		return "";
	}
	else if (pow == 1) {
		return var;
	}
	else if (pow < 0) {
		return var + "^(" + ToString(pow) + ")";
	}
	else {
		return var + "^" + ToString(pow);
	}
}

void ScriptGenerator::SaveOtoc(OtocParameters params, f_type tMin, f_type tMax, int steps, bool parts) {
	IOtoc* otoc = IOtoc::Create(params);

	if (params.ExcludeZeroStates) {
		LOG("computing OTOC for W=" << params.W << ", V=" << params.V 
			<< ", L=" << params.L << " excluding zero energy states...", Verbos);
	}
	else {
		LOG("computing OTOC for W=" << params.W << ", V=" << params.V
			<< ", L=" << params.L, Verbos);
	}

	otoc->Init();

	var_t w = Chop(otoc->ExpectationValueW());
	var_t v = Chop(otoc->ExpectationValueV());
	var_t ww = Chop(otoc->AsymptoticWW());
	var_t vv = Chop(otoc->AsymptoticVV());

	var_t ww1 = Chop(otoc->WW1());
	var_t ww2 = Chop(otoc->WW2());
	var_t vv1 = Chop(otoc->VV1());
	var_t vv2 = Chop(otoc->VV2());

	int precision = std::numeric_limits<f_type>::digits10;
	LOG("<W>= Tr(rho*W) \t=\t" << std::setprecision(precision) << w, Info);
	LOG("<V>= Tr(rho*V) \t=\t" << std::setprecision(precision) << v, Info);
	LOG("<WW> = Tr(rho*W*W^{\\dag}) \t=\t" << std::setprecision(precision) << ww, Info);
	LOG("<VV> = Tr(rho*V*V^{\\dag}) \t=\t" << std::setprecision(precision) << vv, Info);
	LOG("<WW1> = Tr(" << pow2String("rho", params.alpha - params.sigma) << "*W*"
		<< pow2String("rho", 1 - params.alpha + params.sigma) << "*W^{\\dag})\t=\t"
		<< std::setprecision(precision) << ww1, Info);
	LOG("<VV1> = Tr(" << pow2String("rho", params.alpha + params.sigma) << "*V*"
		<< pow2String("rho", 1 - params.alpha - params.sigma) << "*V^{\\dag})\t=\t"
		<< std::setprecision(precision) << vv1, Info);
	LOG("<WW2> = Tr(" << pow2String("rho", params.alpha + params.sigma) << "*W*"
		<< pow2String("rho", 1 - params.alpha - params.sigma) << "*W^{\\dag})\t=\t"
		<< std::setprecision(precision) << ww2, Info);
	LOG("<VV2> = Tr(" << pow2String("rho", params.alpha - params.sigma) << "*V*"
		<< pow2String("rho", 1 - params.alpha + params.sigma) << " V^{\\dag})\t=\t"
		<< std::setprecision(precision) << vv2, Info);

	LOG("Asymptotic value of C(t) = <WW1> * <VV1> + <WW2> * <VV2> = "
		<< std::setprecision(precision) << Chop(ww1 * vv1 + ww2 * vv2) << std::setprecision(6), Info);

	string name = "otoc";
	if (params.alpha != 1.0 || params.sigma != 0.0) {
		name += "_alpha" + ToString(params.alpha) + "sigma" + ToString(params.sigma);
	}

	name += "_" + params.W + "_" + params.V + "_L" + ToString(params.L);

	if (params.M >= 0) {
		name += "M" + ToString(params.M);
	}
	if (params.N != std::numeric_limits<f_type>::infinity()) {
		name += "N" + ToString(params.N);
	}
	name += "beta" + ToString(params.Beta) + "b" + ToString(params.B);

	if (params.ExcludeZeroStates) {
		name += "_nozero";
	}

	LOG("Computing otoc...", Verbos);
	string file = CombinePath(rootFolder, name + ".txt");
	ofstream ofs;
	ofs.open(file.c_str(), std::ofstream::out | std::ofstream::app);

	if (std::abs(tMin) < 1e-10) {
		ofs << "<WW1>" << "\t" << std::setprecision(precision) << ww1.real() << std::endl;
		ofs << "<VV1>" << "\t" << std::setprecision(precision) << vv1.real() << std::endl;
		ofs << "<WW2>" << "\t" << std::setprecision(precision) << ww2.real() << std::endl;
		ofs << "<VV2>" << "\t" << std::setprecision(precision) << vv2.real() << std::endl;
	}

	f_type delta = (tMax - tMin) / steps;
	for (int i = 0; i <= steps; i++) {
		f_type t = tMin + i * delta;
		var_t res;
		if (!parts) {
			var_t res = Chop(otoc->Compute(t));
			LOG("t=" << t << "\tC(t)=" << res.real(), Info);
			ofs << t << "\t" << std::setprecision(precision) << res.real() << std::endl;
		}
		else {
			var_t toc1, otoc1;
			otoc->ComputeParts(t, toc1, otoc1);
			res = toc1 * 2.0 - 2 * otoc1.real();
			toc1 = Chop(toc1);
//			toc2 = Chop(toc2);
			otoc1 = Chop(otoc1);
//			otoc2 = Chop(otoc2);
			res = Chop(res);
			LOG("t=" << std::setw(4) << std::left << t
				<< " C(t)=" << std::setw(8) << std::left << res.real()
				<< " TOC1=" << std::setw(8) << std::left << toc1.real()
//				<< " TOC2=" << std::setw(8) << std::left << toc2.real()
				<< " OTOC1=" << std::setw(6) << std::left << otoc1
//				<< " OTOC2=" << std::setw(6) << std::left  << otoc2
				, Info);
			ofs << std::setw(8) << std::left << t
				<< ' ' << std::setw(20) << std::left << std::setprecision(precision) << res.real()
				<< '\t' << std::setw(20) << std::left << 2 * toc1.real()
				<< '\t' << std::setw(20) << std::left << 2 * otoc1.real() << std::endl;
			if (toc1.imag() != .0) {
				LOG("Imaginary part of TOC(t) is not 0. TOC1=" << toc1, Warning);
			}
		}
		if (res.imag() != .0L) {
			LOG("Imaginary part of C(t) is not 0. otoc=" << res, Warning);
		}
		ofs.flush();
	}

	ofs.close();
	LOG("otoc data are saved to file " << file, Info);
	delete otoc;
}