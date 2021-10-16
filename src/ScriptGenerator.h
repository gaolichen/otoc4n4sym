#pragma once
#include <string>
#include <vector>
#include <fstream>
#include "HamOperator.h"
#include "otoc.h"
#include "NormCalculator.h"
using namespace std;

class ScriptGenerator
{
private:
	string rootFolder;
	void SaveNormToDataFile(int chainLength, int magnonNumber, NormCalculator& calc, int format = 0);
	void SaveMatrix(vector<vector<NExpansionSeries> > &mat, ofstream& ofs);
public:
	ScriptGenerator(string rootFolder);

	void SaveHamMatrix(int chainLength, int magnonNumber);
	void SaveStates(int chainLength, int magnonNumber);
	void SaveSpectrum(int chainLength, int magnonNumber, int N, f_type beta);
	void SaveNormToDataFile(int chainLength, int magnonNumber);
	void SaveNormToDataFile(int chainLength);
	void SaveNormToDataFile2(int chainLength, int magnonNumber);
	void SaveNormToDataFile2(int chainLength);
	void SaveOpMatrix(std::string op, int chainLength);
	void SaveOpMatrix(std::string op, int chainLength, int magnonNumber);
	void SaveOtoc(OtocParameters params, f_type tMin, f_type tMax, int steps, bool parts = true);
};
