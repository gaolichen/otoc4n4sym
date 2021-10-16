#pragma once
#include <map>
#include "MixState.h"
#include "SingleTrace.h"
#include "TraceState.h"
#include "StateCollection.h"
#include "BlockMatrix.h"
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;

extern f_type inverseTempeture;
extern f_type beta_nsym;

class IOperator
{
public:
	static IOperator* Create(std::string op);

	virtual bool PureMagnon() const = 0;

	virtual int MagnonNumberChange() const = 0;

	bool MagnonFixed() const {
		return PureMagnon() && MagnonNumberChange() == 0;
	}

	virtual void ToMatrixN(int L, f_type N, BlockMatrix& res) = 0;

	virtual void ToMatrixN(int L, int M, f_type N, Matrix& res) = 0;

	virtual string ToString() const = 0;
};

class ConstantOperator : public IOperator
{
private:
	var_t _val;
public:
	ConstantOperator(var_t val) {
		_val = val;
	}

	virtual bool PureMagnon() const {
		return true;
	};

	virtual int MagnonNumberChange() const {
		return 0;
	}

	var_t Value() const {
		return _val;
	}

	void MultiplyBy(var_t factor) {
		_val *= factor;
	}

	virtual string ToString() const {
		if (IsZero(_val.imag())) {
			return ::ToString(_val.real());
		}
		else {
			return ::ToString(_val);
		}
	}

	virtual void ToMatrixN(int L, f_type N, BlockMatrix& res) {
		res.init(L);
		for (int i = 0; i <= L; i++) {
			int size = StateCollection::Inst()->StateNumber(L, i);
			res.block(i, i) = Array::Constant(size, _val).matrix().asDiagonal();
		}
	}

	virtual void ToMatrixN(int L, int M, f_type N, Matrix& res) {
		int size = StateCollection::Inst()->StateNumber(L, M);
		res = Array::Constant(size, _val).matrix().asDiagonal();
	}
};


class NormalOrderedOperator : public IOperator
{
private:
	IOperator* _op;
public:
	NormalOrderedOperator(IOperator* op) {
		_op = op;
	}

	~NormalOrderedOperator() {
		if (_op != NULL) {
			delete _op;
		}
	}

	virtual int MagnonNumberChange() const {
		return _op->MagnonNumberChange();
	}

	bool PureMagnon() const {
		return _op->PureMagnon();
	}

	virtual string ToString() const {
		return "~" + _op->ToString();
	}

	virtual void ToMatrixN(int L, f_type N, BlockMatrix& res);

	virtual void ToMatrixN(int L, int M, f_type N, Matrix& res);
};

class ScalarMultiplyOperator : public IOperator
{
private:
	IOperator* _op;
	var_t _factor;
public:
	ScalarMultiplyOperator(IOperator* op, var_t factor) {
		_factor = factor;
		_op = op;
	}

	~ScalarMultiplyOperator() {
		if (_op != NULL) {
			delete _op;
		}
	}

	var_t Factor() const {
		return this->_factor;
	}

	virtual int MagnonNumberChange() const {
		return _op->MagnonNumberChange();
	}

	bool PureMagnon() const {
		return _op->PureMagnon();
	}

	void MultiplyBy(f_type val) {
		this->_factor *= val;
	}

	virtual string ToString() const {
		if (IsZero(_factor.imag())) {
			if (_factor.real() == 1) {
				return _op->ToString();
			}
			else if (_factor.real() == -1) {
				return "-" + _op->ToString();
			}
			return ::ToString(_factor.real()) + "*" + _op->ToString();
		}
		return ::ToString(_factor) + "*" + _op->ToString();
	}

	virtual void ToMatrixN(int L, f_type N, BlockMatrix& res) {
		_op->ToMatrixN(L, N, res);
		res *= _factor;
	}

	virtual void ToMatrixN(int L, int M, f_type N, Matrix& res) {
		_op->ToMatrixN(L, M, N, res);
		res *= _factor;
	}
};

class SumOperator : public IOperator
{
private:
	std::vector<IOperator*> _ops;
public:

	~SumOperator() {
		for (int i = 0; i < _ops.size(); i++) {
			delete _ops[i];
		}
	}

	virtual bool PureMagnon() const {
		for (int i = 0; i < _ops.size(); i++) {
			if (!_ops[i]->PureMagnon()) {
				return false;
			}
			if (_ops[i]->MagnonNumberChange() != _ops[0]->MagnonNumberChange()) {
				return false;
			}
		}
		return true;
	}

	virtual int MagnonNumberChange() const {
		if (_ops.size() > 0 && PureMagnon()) {
			return _ops[0]->MagnonNumberChange();
		}
		else {
			return 0;
		}
	}

	virtual void ToMatrixN(int L, f_type N, BlockMatrix& res) {
		res.init(L);
		for (int i = 0; i < _ops.size(); i++) {
			BlockMatrix tmp;
			_ops[i]->ToMatrixN(L, N, tmp);
			res += tmp;
		}
	}

	virtual void ToMatrixN(int L, int M, f_type N, Matrix& res) {
		int size = StateCollection::Inst()->StateNumber(L, M);
		res.resize(size, size);
		res.fill(0.0);
		for (int i = 0; i < _ops.size(); i++) {
			Matrix tmp;
			_ops[i]->ToMatrixN(L, M, N, tmp);
			res += tmp;
		}
	}

	void AddOperator(IOperator* op) {
		_ops.push_back(op);
	}

	virtual string ToString() const {
		if (_ops.size() == 0) {
			return "";
		}

		std::string ret = _ops[0]->ToString();
		for (int i = 1; i < _ops.size(); i++) {
			ScalarMultiplyOperator* op = dynamic_cast<ScalarMultiplyOperator*>(_ops[i]);
			ConstantOperator *cop = dynamic_cast<ConstantOperator*>(_ops[i]);
			if (op == NULL && cop == NULL) {
				ret += "+" + _ops[i]->ToString();
			}
			else {
				var_t val;
				if (op != NULL) {
					val = op->Factor();
				}
				else {
					val = cop->Value();
				}
				if (!IsZero(val.imag()) || val.real() > 0) {
					ret += "+" + _ops[i]->ToString();
				}
				else {
					ret += _ops[i]->ToString();
				}
			}
		}

		return ret;
	}
};

class ExponentialOperator : public IOperator
{
private:
	IOperator* _op = NULL;
public:
	ExponentialOperator(IOperator* op) {
		_op = op;
		if (!_op->PureMagnon()) {
			throw OtocException("Exponentiating the operator " + op->ToString() + " is not supported.");
		}
	}

	~ExponentialOperator() {
		if (_op != NULL) {
			delete _op;
		}
	}

	virtual int MagnonNumberChange() const {
		return 0;
	}

	bool PureMagnon() const {
		return _op->MagnonFixed();
	}

	virtual string ToString() const {
		return "exp(" + _op->ToString() + ")";
	}

	virtual void ToMatrixN(int L, f_type N, BlockMatrix& res) {
		BlockMatrix mat;
//		LOG("ExponentialOperator ToMatrixN 1", Verbos);
		_op->ToMatrixN(L, N, mat);
//		LOG("ExponentialOperator ToMatrixN 2", Verbos);
		res = mat.exp();
//		LOG("ExponentialOperator ToMatrixN 3", Verbos);
	}

	virtual void ToMatrixN(int L, int M, f_type N, Matrix& res) {
		Matrix mat;
		_op->ToMatrixN(L, M, N, mat);
		res = mat.exp();
	}
};

class ISimpleOperator : public IOperator
{
private:
	void ToMatrixN(int L, int M1, int M2, f_type N, BlockMatrix& res);
public:
	static ISimpleOperator* Create(std::string opStr);
	virtual bool PureMagnon() const {
		return true;
	}

	virtual void ApplyOn(const TraceState& state, MixState& res) = 0;

	virtual void ToMatrixN(int L, f_type N, BlockMatrix& res);

	virtual void ToMatrixN(int L, int M, f_type N, Matrix& res);
};

class TraceNumberOperator : public ISimpleOperator
{
public:
	virtual int MagnonNumberChange() const {
		return 0;
	}

	virtual string ToString() const {
		return "NTr";
	}

	virtual void ApplyOn(const TraceState& state, MixState& res);
};

class MagnonTraceNumberOperator : public ISimpleOperator
{
public:
	virtual int MagnonNumberChange() const {
		return 0;
	}

	virtual string ToString() const {
		return "NTrX";
	}

	virtual void ApplyOn(const TraceState& state, MixState& res);
};

class TraceOperator : public ISimpleOperator
{
private:
	map<SingleTrace, MixState> cache1;
	map<pair<SingleTrace, SingleTrace>, MixState> cache2;
	MixState& ApplyOnPrivate(const SingleTrace& single);
	MixState& ApplyOnPrivate(const SingleTrace& single1, const SingleTrace& single2);
protected:
	void AddState(TraceState& state, NExpansionSeries& coef, MixState& res);
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res) = 0;
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res) = 0;
public:
	TraceOperator();
	virtual void ApplyOn(const TraceState& state, MixState& res);
	friend ostream& operator<<(ostream& os, const TraceOperator& op);
	friend bool operator == (const TraceOperator& op1, const TraceOperator& op2);
};

class DilatationOperator : public TraceOperator
{
protected:
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res);
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res);
public:
    DilatationOperator();

	virtual int MagnonNumberChange() const {
		return 0;
	}

	void ToMatrix(int L, int M, vector<vector<NExpansionSeries> >& res);

	void ToMatrix(int L, vector<vector<NExpansionSeries> >& res);

	void ToMatrixN(int L, int M, f_type N, f_type beta, Matrix& res);

    virtual string ToString() const {
        return "D";
    }
};

class SingleTraceOperator : public TraceOperator
{
protected:
	int _annihilator;
	int _creator;
public:
	virtual int MagnonNumberChange() const {
		return BitCount(_creator) - BitCount(_annihilator);
	}
};

class TwoFieldOperator : public SingleTraceOperator
{
private:
protected:
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res);
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res) {
		// do nothing
	}
public:
	TwoFieldOperator(int creator, int annihilator);

	TwoFieldOperator(std::string op);

	virtual string ToString() const
	{
		std::string ret = "Tr(";
		ret += (_creator == 1 ? "X" : "Z");
		ret += (_annihilator == 1 ? "X" : "Z");
		ret += ")";
		return ret;
	}
};

class FourFieldOperator : public SingleTraceOperator
{
private:
protected:
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res);
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res);
public:
	FourFieldOperator(int creator, int annihilator);

	FourFieldOperator(std::string op);

	virtual string ToString() const;
};

class FourFieldOperatorB : public SingleTraceOperator
{
private:
protected:
	virtual void ApplyOnSingle(const SingleTrace& single, MixState& res);
	virtual void ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res);
public:
	FourFieldOperatorB(int creator, int annihilator);

	FourFieldOperatorB(std::string op);

	virtual string ToString() const;
};