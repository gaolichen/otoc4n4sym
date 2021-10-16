#include "Spectrum.h"

Spectrum::Spectrum()
{
    op = new DilatationOperator();
}

Spectrum::~Spectrum()
{
    delete op;
}

void Spectrum::Eigenvalues(int L, int M, int N, f_type beta, vector<var_t>& ret)
{
    vector<vector<NExpansionSeries> > res;
    op->ToMatrix(L, M, res);
    Matrix mat(res.size(), res.size());
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            mat(i, j) = res[i][j].ToNumerical(N, beta);
        }
    }
    
    Vector eigens = mat.eigenvalues();
    ret.resize(eigens.size());
    for (int i = 0; i < eigens.size(); i++) {
        ret[i] = Chop(eigens[i]);
    }
}
