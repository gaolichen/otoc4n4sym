#pragma once
#include "BitUtility.h"
#include "SingleTrace.h"
#include "StateId.h"
#include "TraceState.h"
#include "PowerSeries.h"
#include <map>

//typedef int nType;
typedef i64 nType;

class MixedTrace : public ITrace
{
private:
	int _maskLeft = 0;
	int _maskRight = 0;
public:
	MixedTrace();

	MixedTrace(int mask1, int mask2);

	MixedTrace(int trace1, int bit1, int trace2, int bit2);

	int LeftBitNumber() const;
	
	int RightBitNumber() const;

	int LeftTrace() const {
		return TRACE_BITS(_maskLeft);
	}

	int RightTrace() const {
		return TRACE_BITS(_maskRight);
	}

	bool Contractable() const {
		return LeftBitNumber() > 0 && RightBitNumber() > 0;
	}

	bool IsEmpty() const {
		return _maskLeft == 0 && _maskRight == 0;
	}

	int RightBit(int index) const;

	int LeftBit(int index) const;

	virtual int Bit(int index) const;

	virtual int BitNumber() const;

	i64 Hash() const;

//	void Split(int pos, SingleTrace& a, SingleTrace& b) const;
//	void Split(int pos1, int pos2, SingleTrace& a, SingleTrace& b, SingleTrace& c) const;

//	void Contract(ContractState& cs) const;
	bool operator< (const MixedTrace& other) const;

	std::string ToString() const;

	friend ostream& operator<<(ostream& os, const MixedTrace& mt);
};

class MixedMultiTrace
{
private:
	StateId _left;
	MixedTrace* _mixed = NULL;
	StateId _right;
public:
	MixedMultiTrace();

	MixedMultiTrace(const MixedMultiTrace& mmt);

	MixedMultiTrace(const StateId& left, const StateId& right, const MixedTrace& mixed);

	MixedMultiTrace(const TraceState& left, const TraceState& right);

	MixedMultiTrace(const StateId& left, const StateId& right);

	MixedMultiTrace(const TraceState& left, const TraceState& right, const MixedTrace& mixed);

	~MixedMultiTrace();

	bool operator< (const MixedMultiTrace& other) const;

	const StateId& Left() const {
		return _left;
	}

	const StateId& Right() const {
		return _right;
	}

	MixedTrace* Mixed() const {
		return _mixed;
	}

	template<class T> void Normalize(PowerSeries<T>& coef)
	{
		;
	}

	std::string ToString() const;

	friend ostream& operator<<(ostream& os, const MixedMultiTrace& mmt);

	pair<i64, i64> NormHash() const;
};
/*
class Contractable
{
private:
public:
	virtual bool SelfContract(ContractState& res) = 0;
	virtual bool Contract(const Contractable& other, ContractState& res) = 0;
};*/

class ContractState
{
private:
	std::map<MixedMultiTrace, PowerSeries<nType> > coefficients;
public:
	ContractState();

	void AddState(const MixedMultiTrace& mm, const PowerSeries<nType>& coefficient);

	void Clear() {
		coefficients.clear();
	}

	PowerSeries<nType> Get(MixedMultiTrace& mt) {
		auto it = coefficients.find(mt);
		if (it != coefficients.end()) {
			return it->second;
		}
		return PowerSeries<nType>();
	}

	void Extend(const TraceState& left, const TraceState& right);
	void Merge(const ContractState& state);
	void Merge(const ContractState& state, const PowerSeries<nType>& prefactor);

	std::map<MixedMultiTrace, PowerSeries<nType> >::const_iterator Begin() const;
	std::map<MixedMultiTrace, PowerSeries<nType> >::const_iterator End() const;

	int Size() const {
		return coefficients.size();
	}

	friend ostream& operator<<(ostream& os, const ContractState& cs);

};
