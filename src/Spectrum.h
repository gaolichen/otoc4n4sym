#pragma once
#include <map>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "HamOperator.h"

using namespace std;

class Spectrum
{
private:
    DilatationOperator*op;
public:
    Spectrum();
    
    ~Spectrum();
    
    void Eigenvalues(int L, int M, int N, f_type beta, vector<var_t>& ret);
};
