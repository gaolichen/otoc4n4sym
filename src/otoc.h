#pragma once
#include "HamOperator.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::Matrix<f_type, Eigen::Dynamic, Eigen::Dynamic> MatrixXF;

struct OtocParameters {
	std::string W;
	std::string V;
	int L;
	int M = -1;
	f_type N;
	f_type Beta;
	f_type B;
	char ComputeNorm = 0;
	char ExcludeZeroStates = 0;
	f_type alpha = 1.0;
	f_type sigma = 0.0;
};

class IOtoc {
public:
	static IOtoc* Create(OtocParameters params);

	virtual var_t Compute(f_type t) = 0;

	virtual void ComputeParts(f_type t, var_t& toc, var_t& otoc) = 0;

	virtual void Init() = 0;

	virtual var_t AsymptoticWW() = 0;

	virtual var_t AsymptoticVV() = 0;

	virtual var_t WW1() const = 0;
	virtual var_t WW2() const = 0;
	virtual var_t VV1() const = 0;
	virtual var_t VV2() const = 0;

	virtual var_t ExpectationValueW() = 0;

	virtual var_t ExpectationValueV() = 0;
};

template<typename T> class OtocBase : public IOtoc{
protected:
	OtocParameters _params;

	T* _W = NULL;
	T* _V = NULL;
	T* _H = NULL;
	T* _G = NULL;

	// expBH1 = G * exp(-b(alpha-sigma)*H)
	T _expBH1;
	// expBH2 = exp(-b(1-alpha-sigma)*H)*G^-1
	T _expBH2;

	// expHSigma = exp(-b*sigma * H) 
	Array _expHSigma;

	Array _eigenValues;

	var_t _normalization = .0L;
	var_t _ww;
	var_t _vv;

	var_t _ww1;
	var_t _ww2;
	var_t _vv1;
	var_t _vv2;

	var_t _evW;
	var_t _evV;

	virtual void PostInit() = 0;

	virtual void InitHam() = 0;

	virtual void InitNorm() = 0;

	virtual void InitOperators() = 0;
public:
	OtocBase(OtocParameters params) {
		this->_params = params;
	}

	~OtocBase() {
		if (_H != NULL) {
			delete _H;
		}
		if (_G != NULL) {
			delete _G;
		}
		if (_W != NULL) {
			delete _W;
		}
		if (_V != NULL) {
			delete _V;
		}
	}

	virtual void Init() {
		LOG("Loading W and V operators...", Verbos);
		InitOperators();

		LOG("Generating Hamiltonian matrix...", Verbos);
		InitHam();

		LOG("Loading norm matrix...", Verbos);
		InitNorm();

		PostInit();
	}

	virtual var_t WW1() const {
		return _ww1;
	}

	virtual var_t WW2() const {
		return _ww2;
	}

	virtual var_t VV1() const {
		return _vv1;
	}

	virtual var_t VV2() const {
		return _vv2;
	}

	virtual var_t AsymptoticWW() {
		return _ww;
	}

	virtual var_t AsymptoticVV() {
		return _vv;
	}

	virtual var_t ExpectationValueW() {
		return _evW;
	}

	virtual var_t ExpectationValueV() {
		return _evV;
	}
};

class FixedMOtoc : public OtocBase<Matrix> {
protected:
	virtual void InitHam();

	virtual void InitNorm();

	virtual void InitOperators();

	virtual void PostInit();

	void ComputeWt(f_type t, Matrix& wt);
public:
	FixedMOtoc(OtocParameters params) : OtocBase(params) {
	}

	virtual var_t Compute(f_type t) {
		Matrix Wt;
		ComputeWt(t, Wt);

		Matrix mat2;
		if (_params.alpha == 0) {
			mat2 = Wt * (*_V) - (*_V) * Wt;
		}
		else {
			mat2 = Wt * _expHSigma.matrix().asDiagonal() * (*_V) - (*_V) * _expHSigma.matrix().asDiagonal() * Wt;
		}
		
		return (_expBH1 * mat2 * _expBH2 * mat2.adjoint()).trace() / _normalization;
	}

	virtual void ComputeParts(f_type t, var_t& toc, var_t& otoc) {
		Matrix Wt;
		ComputeWt(t, Wt);

		Matrix wv, vw;
		if (_params.alpha == 0) {
			wv = Wt * (*_V);
			vw = (*_V) * Wt;
		}
		else {
			wv = Wt * _expHSigma.matrix().asDiagonal() * (*_V);
			vw = (*_V) * _expHSigma.matrix().asDiagonal() * Wt;
		}

		Matrix mat = _expBH1 * wv * _expBH2;
		toc = (mat * wv.adjoint()).trace() / _normalization;
		otoc = (mat * vw.adjoint()).trace() / _normalization;
		
//		mat.noalias() = _expBH1 * vw * _expBH2;
//		var_t toc2 = (mat * vw.adjoint()).trace() / _normalization;
//		var_t otoc2 = (mat * wv.adjoint()).trace() / _normalization;
	}
};

class NonFixedMOtoc : public OtocBase<BlockMatrix> {
protected:
//	Array _noZeroProjector;

	virtual void PostInit();

	virtual void InitHam();

	virtual void InitNorm();

	virtual void InitOperators();

	var_t Vev(const BlockMatrix& mat);

	void ComputeWt(f_type t, BlockMatrix& wt);
public:
	NonFixedMOtoc(OtocParameters params) : OtocBase(params) {
	}

	virtual var_t Compute(f_type t) {
		BlockMatrix Wt;
		ComputeWt(t, Wt);

		BlockMatrix mat2;
		if (_params.sigma == 0) {
			mat2 = Wt * (*_V) - (*_V) * Wt;
		}
		else {
			BlockMatrix mat1;
			mat2.fromMixing(Wt, _expHSigma, *_V);
			mat1.fromMixing(*_V, _expHSigma, Wt);
			mat2 -= mat1;
		}

		return Vev(mat2) / _normalization;
	}

	virtual void ComputeParts(f_type t, var_t& toc, var_t& otoc) {
		BlockMatrix Wt;
		ComputeWt(t, Wt);

		BlockMatrix wv;
		BlockMatrix vw;
		if (_params.sigma == 0) {
			wv = Wt * (*_V);
			vw = (*_V) * Wt;
		}
		else {
			wv.fromMixing(Wt, _expHSigma, *_V);
			vw.fromMixing(*_V, _expHSigma, Wt);
		}

		BlockMatrix mat1 = _expBH1 * wv * _expBH2;
//		BlockMatrix mat2 = _expBH1 * vw * _expBH2;

		wv.adjointInPlace();
		vw.adjointInPlace();

		toc = (mat1 * wv).trace() / _normalization;
		otoc = (mat1 * vw).trace() / _normalization;
//		var_t toc2 = (mat2 * vw).trace() / _normalization;
//		var_t otoc2 = (mat2 * wv).trace() / _normalization;

//		if (abs(Chop(toc2 - toc)) > 1e-10 || abs(Chop(otoc2 - otoc).real()) > 1e-10 || abs(Chop(otoc2 + otoc).imag()) > 1e-10) {
//			std::cout << "toc2 = " << Chop(toc2) << ", otoc2=" << Chop(otoc2) << std::endl;
//		}
	}
};

void LoadNormFromFile(int L, int M, f_type N, Matrix& res);

void ComputeNorm(int L, int M, f_type N, Matrix& res);

void RemoveZeroStates(Array& arr, Array eigenValues);