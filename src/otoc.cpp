#include "otoc.h"
#include <fstream>
#include "NormCalculator.h"

#define IS_ZERO(x) (std::abs(x) < 1e-10)

std::vector<i64> StringToVector(std::string str) {
	std::istringstream iss(str);
	std::vector<i64> ret;
	i64 val;
	while (iss >> val) {
		ret.push_back(val);
	}

	return ret;
}

void LoadNormFromFile(int L, int M, f_type N, Matrix& res) {
	std::string file = "./otocdata/norm_L" + ToString(L) + "M" + ToString(M) + "_su.txt";
	std::ifstream fin(file.c_str());
	if (!fin.is_open()) {
		throw OtocException("cannot find file " + file + ". Please first generate norm matrix using otoc -n -su command.");
	}

	res.fill(.0L);
	
	f_type nn = (N == std::numeric_limits<f_type>::infinity() ? N : N * N);
	int size;
	fin >> size;
//	std::cout << "size=" << size << std::endl;
	std::string line;
	std::getline(fin, line);
	for (int i = 0; i < size; i++) {
		std::getline(fin, line);
		for (int j = i; j < size; j++) {
			std::getline(fin, line);
			vector<i64> val = StringToVector(line);
			if (val.size() == 1) {
				continue;
			}
			if (val.size() == 0) {
				std::cout << "error: val.size() should greater than 1" << std::endl;
			}
			f_type norm = .0L;
			if (nn == std::numeric_limits<f_type>::infinity()) {
				norm = val[1];
			}
			else {
				for (int k = val.size() - 1; k > 0; k--) {
					norm /= nn;
					norm += val[k];
				}
			}

			if (val[0]) {
				norm /= N;
			}

			res(i, j) = norm;
			if (i != j) {
				res(j, i) = norm;
			}
		}
	}

	fin.close();
}

void ComputeNorm(int L, int M, f_type N, Matrix& res) {
	StateCollection* inst = StateCollection::Inst();
	auto bIt = inst->Begin(L, M);
	auto eIt = inst->End(L, M);
	DfsNormCalculator calc;

	for (auto it1 = bIt; it1 != eIt; ++it1) {
		for (auto it2 = it1; it2 != eIt; ++it2) {
			PowerSeries<nType> poly = calc.Calculate(*it1, *it2);
			poly.ShiftOrder(-L);
			res(it1 - bIt, it2 - bIt) = poly.ToNumerical(N);
			if (it2 != it1) {
				res(it2 - bIt, it1 - bIt) = res(it1 - bIt, it2 - bIt);
			}
		}
	}
}

void FixedMOtoc::InitOperators() {
	// init W and V
	IOperator* opW = IOperator::Create(_params.W);
	IOperator* opV = IOperator::Create(_params.V);
	LOG("operator W=" << opW->ToString(), Verbos);
	LOG("operator V=" << opV->ToString(), Verbos);

	if (!opW->MagnonFixed() || !opV->MagnonFixed()) {
		throw OtocException("Operator " + _params.W + " and " + _params.V + " do not both preserve magnon number.");
	}

	int	matSize = StateCollection::Inst()->StateNumber(_params.L, _params.M);

	_W = new Matrix(matSize, matSize);
	_V = new Matrix(matSize, matSize);
	opW->ToMatrixN(_params.L, _params.M, _params.N, *_W);
	opV->ToMatrixN(_params.L, _params.M, _params.N, *_V);
	
	delete opW;
	delete opV;
}

void FixedMOtoc::InitHam() {
	DilatationOperator ham;
	StateCollection* inst = StateCollection::Inst();
	int size = inst->StateNumber(_params.L, _params.M);

	_H = new Matrix(size, size);
	_H->fill(0.0);
	ham.ToMatrixN(_params.L, _params.M, _params.N, _params.Beta, *_H);
}

void FixedMOtoc::InitNorm() {
	StateCollection* inst = StateCollection::Inst();
	int size = inst->StateNumber(_params.L, _params.M);

	_G = new Matrix(size, size);
	if (!_params.ComputeNorm) {
		LoadNormFromFile(_params.L, _params.M, _params.N, *_G);
	}
	else {
		ComputeNorm(_params.L, _params.M, _params.N, *_G);
	}
}

void FixedMOtoc::PostInit() {
	LOG("Diagonalizing dilitation matrices...", Verbos);
	Eigen::ComplexEigenSolver<Matrix> solver(*_H, true);
	_eigenValues = solver.eigenvalues().array();
	Matrix U = solver.eigenvectors();
//	std::cout << "U.imag()=" << (U - U.conjugate()).norm() << std::endl;
	Matrix Uinv = U.inverse();
	Array expHdiag = (_eigenValues * (-_params.B)).exp();

	LOG("Transform W, V, and norm matrix into basis of energy eigenstates...", Verbos);
	*_W = Uinv * (*_W) * U;
	*_V = Uinv * (*_V) * U;
	*_G = U.adjoint() * (*_G) * U;

	LOG("Computing exponential of Hamiltonian...", Verbos);
	_expHSigma = (_eigenValues * (-_params.B * _params.sigma)).exp();
	Array exp1 = (_eigenValues * (-_params.B) * (_params.alpha - _params.sigma)).exp();
	Array exp2 = (_eigenValues * (-_params.B) * (1 - _params.alpha - _params.sigma)).exp();
	Matrix gInv = _G->inverse();
	_expBH1 = (*_G) * exp1.matrix().asDiagonal();
	_expBH2 = exp2.matrix().asDiagonal() * gInv;

	if (_params.ExcludeZeroStates) {
		LOG("Restricting W and V in subspace of nonzero energy states...", Verbos);
		// restrict W and V to nonzero energy subspace.
		std::vector<int> zeros;
		for (int i = 0; i < _eigenValues.size(); i++) {
			if (IS_ZERO(_eigenValues(i).real())) {
				zeros.push_back(i);
//				expHdiag(i) = 0;
//				_expHSigma(i) = 0;
//				_W->row(i) = Vector::Zero(_W->cols());
//				_W->col(i) = Vector::Zero(_W->rows());

//				_V->row(i) = Vector::Zero(_V->cols());
//				_V->col(i) = Vector::Zero(_V->rows());
			}
		}

		for (int i = 0; i < zeros.size(); i++) {
			for (int j = 0; j < zeros.size(); j++) {
				(*_W)(zeros[i], zeros[j]) = 0;
				(*_V)(zeros[i], zeros[j]) = 0;
			}
		}
	}

	_normalization = expHdiag.sum();
	LOG("Tr(e^(-b*H)) = " << _normalization, Verbos);

	LOG("Computing <WW> and <VV>...", Verbos);
	Array expBHA = (_eigenValues * (-_params.B)).exp();
	Matrix expBH = (*_G) * expBHA.matrix().asDiagonal();
	_ww = (expBH * (*_W) * gInv * _W->adjoint()).trace() / _normalization;
	_vv = (expBH * (*_V) * gInv * _V->adjoint()).trace() / _normalization;

	_evW = (expBHA.matrix().asDiagonal() * (*_W)).trace() / _normalization;
	_evV = (expBHA.matrix().asDiagonal() * (*_V)).trace() / _normalization;

	Array tmp1 = (_eigenValues * (-_params.B) * (1 - _params.alpha + _params.sigma)).exp();

	_ww1 = (_expBH1 * (*_W) * tmp1.matrix().asDiagonal() * gInv * _W->adjoint()).trace() / _normalization;
	_vv2 = (_expBH1 * (*_V) * tmp1.matrix().asDiagonal() * gInv * _V->adjoint()).trace() / _normalization;

	Array tmp2 = (_eigenValues * (-_params.B) * (_params.alpha + _params.sigma)).exp();
	_ww2 = ((*_G) * tmp2.matrix().asDiagonal() * (*_W) * _expBH2 * _W->adjoint()).trace() / _normalization;
	_vv1 = ((*_G) * tmp2.matrix().asDiagonal() * (*_V) * _expBH2 * _V->adjoint()).trace() / _normalization;

	delete _H;
	delete _G;
	_G = _H = NULL;
}

void FixedMOtoc::ComputeWt(f_type t, Matrix& Wt) {
	var_t time(0.0, -t);
	Array exp1 = (_eigenValues * time).exp();
	if (_params.ExcludeZeroStates) {
		RemoveZeroStates(exp1, _eigenValues);
	}

	// since eigenvalues are real, eigenvalues * time is pure imaginary
	// so for its exponential, conjugate and inverse are the same.
	// but i think conjugate operation should be cheaper than inverse operation
	Array exp2 = exp1.conjugate();
//	std::cout << "diff norm=" << (exp1.inverse() - exp2).matrix().norm() << std::endl;

	Wt = exp1.matrix().asDiagonal() * (*_W) * exp2.matrix().asDiagonal();
}

void NonFixedMOtoc::InitOperators() {
	// init W and V
	IOperator* opW = IOperator::Create(_params.W);
	IOperator* opV = IOperator::Create(_params.V);
	LOG("operator W=" << opW->ToString(), Verbos);
	LOG("operator V=" << opV->ToString(), Verbos);

	_W = new BlockMatrix(_params.L);
	_V = new BlockMatrix(_params.L);

	opW->ToMatrixN(_params.L, _params.N, *_W);
	opV->ToMatrixN(_params.L, _params.N, *_V);
	delete opW;
	delete opV;
}

void NonFixedMOtoc::InitHam() {
	_H = new BlockMatrix(_params.L);
	DilatationOperator ham;
	StateCollection* inst = StateCollection::Inst();

	for (int i = 0; i <= _params.L; i++) {
		_H->fill(i, i);
		ham.ToMatrixN(_params.L, i, _params.N, _params.Beta, _H->block(i, i));
	}
}

void NonFixedMOtoc::InitNorm() {
	_G = new BlockMatrix(_params.L);
	for (int i = 0; i <= _params.L; i++) {
		_G->fill(i, i);
		if (!_params.ComputeNorm) {
			LoadNormFromFile(_params.L, i, _params.N, _G->block(i, i));
		}
		else {
			ComputeNorm(_params.L, i, _params.N, _G->block(i, i));
		}
	}
}

void NonFixedMOtoc::PostInit() {
	BlockMatrix U, Uinv;
	LOG("Diagonalizing dilatation matrices...", Verbos);
	_H->diagonalize(_eigenValues, U, Uinv);
	LOG("Transform W, V, and the norm matrix into basis of energy eigenstates...", Verbos);
	_W->mixWith(Uinv, U);
	_V->mixWith(Uinv, U);
	BlockMatrix Uadj = U.adjoint();
	_G->mixWith(Uadj, U);

	std::vector<int> zeros;
//	_noZeroProjector = Array::Constant(_eigenValues.size(), 1.0);
	for (int i = 0; i < _eigenValues.size(); i++) {
		if (IS_ZERO(_eigenValues(i).real())) {
//			_noZeroProjector(i) = 0.0;
			zeros.push_back(i);
		}
	}

	LOG("Computing exponential of Hamiltonian...", Verbos);
	Array expBHDiag = (_eigenValues * (-_params.B)).exp();
	Array expBHDiag1 = (_eigenValues * (-_params.B) * (_params.alpha - _params.sigma)).exp();
	Array expBHDiag2 = (_eigenValues * (_params.B) * (1.0 - _params.alpha - _params.sigma)).exp();

	_expBH1 = (*_G) * expBHDiag1;
	_expBH2 = (*_G) * expBHDiag2;
	_expBH2.inverseInPlace();
	_expHSigma = (_eigenValues * (-_params.B * _params.sigma)).exp();

	_normalization = Chop(expBHDiag.sum());
	if (!_params.ExcludeZeroStates) {
//		_normalization = Chop(expBHDiag.sum());
	}
	else {
		LOG("Restricting W and V in subspace of nonzero energy states...", Verbos);
		// restrict W to nonzero energy subspace. 
		_W->setZeroAt(zeros);
		_V->setZeroAt(zeros);
//		_W->mixWith(_noZeroProjector, _noZeroProjector);
//		_V->mixWith(_noZeroProjector, _noZeroProjector);
//		_expBH1.mixWith(_noZeroProjector, _noZeroProjector);
//		_expBH2.mixWith(_noZeroProjector, _noZeroProjector);
//		_normalization = Chop(expBHDiag.matrix().dot(_noZeroProjector.matrix()));
	}
	LOG("Tr(e^(-b*H))=" << _normalization, Verbos);

//	std::cout << "expHSigma=" << std::endl << _expHSigma << std::endl;
//	std::cout << "expBH1=" << std::endl << _expBH1 << std::endl;
//	std::cout << "expBH2=" << std::endl << _expBH2 << std::endl;

	LOG("Computing <WW> and <VV>...", Verbos);
	BlockMatrix expBH = (*_G) * expBHDiag;
	BlockMatrix gInv = _G->inverse();
	_ww = (expBH * (*_W) * gInv * _W->adjoint()).trace() / _normalization;
	_vv = (expBH * (*_V) * gInv * _V->adjoint()).trace() / _normalization;
	_evW = expBHDiag.matrix().dot(_W->diagonal().matrix()) / _normalization;
	_evV = expBHDiag.matrix().dot(_V->diagonal().matrix()) / _normalization;

	Array tmp1 = (_eigenValues * (_params.B) * (1 - _params.alpha + _params.sigma)).exp();
	BlockMatrix tmpExp1 = (*_G) * tmp1;
	tmpExp1.inverseInPlace();

	_ww1 = (_expBH1 * (*_W) * tmpExp1 * _W->adjoint()).trace() / _normalization;
	_vv2 = (_expBH1 * (*_V) * tmpExp1 * _V->adjoint()).trace() / _normalization;

	Array tmp2 = (_eigenValues * (-_params.B) * (_params.alpha + _params.sigma)).exp();
	_ww2 = ((*_G) * tmp2 * (*_W) * _expBH2 * _W->adjoint()).trace() / _normalization;
	_vv1 = ((*_G) * tmp2 * (*_V) * _expBH2 * _V->adjoint()).trace() / _normalization;

	delete _H;
	delete _G;
	_H = NULL;
	_G = NULL;
}

var_t NonFixedMOtoc::Vev(const BlockMatrix& mat) {
	//	BlockMatrix tmp = _expBH * mat * _Ginv * mat.adjoint() * (*_G);
	BlockMatrix tmp = _expBH1 * mat * _expBH2 * mat.adjoint();
	//	if (_params.ExcludeZeroStates) {
	//		tmp.mixWith(_noZeroProjector, _noZeroProjector);
	//	}
	return tmp.trace();
}


void NonFixedMOtoc::ComputeWt(f_type t, BlockMatrix& wt) {
	var_t time(0.0, -t);
	Array exp1 = (_eigenValues * time).exp();
	if (_params.ExcludeZeroStates) {
//		RemoveZeroStates(exp1, _eigenValues);
	}
	
	// since eigenvalues are real, eigenvalues * time is pure imaginary
	// so for its exponential, conjugate and inverse are the same.
	// but i think conjugate operation should be cheaper than inverse operation
	Array exp2 = exp1.conjugate();

	wt.fromMixing(exp1, *_W, exp2);
}

IOtoc* IOtoc::Create(OtocParameters params) {
	if (params.M >= 0) {
		return new FixedMOtoc(params);
	}
	else {
		return new NonFixedMOtoc(params);
	}
}

void RemoveZeroStates(Array& arr, Array eigenValues) {
	for (int i = 0; i < eigenValues.size(); i++) {
		if (std::abs(eigenValues(i).real()) < 1e-10) {
			arr(i) = 0.0;
		}
	}
}