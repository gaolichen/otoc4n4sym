#include "MixedTrace.h"
#include "SingleTrace.h"
#include "StateCollection.h"
#define GET_BIT(mask, index) ((mask) & (1<<(index))) != 0 ? 1 : 0

MixedTrace::MixedTrace() {
}

MixedTrace::MixedTrace(int mask1, int mask2) {
	_maskLeft = mask1;
	_maskRight = mask2;
}

MixedTrace::MixedTrace(int trace1, int bit1, int trace2, int bit2) {
	_maskLeft = BuildMask(trace1, bit1);
	_maskRight = BuildMask(trace2, bit2);
}

bool MixedTrace::operator< (const MixedTrace& other) const {
	if (_maskLeft != other._maskLeft) {
		return _maskLeft < other._maskLeft;
	} else {
		return _maskRight < other._maskRight;
	}
}

int MixedTrace::LeftBitNumber() const {
	return BIT_LENGTH(_maskLeft);
}

int MixedTrace::RightBitNumber() const {
	return BIT_LENGTH(_maskRight);
}

int MixedTrace::RightBit(int index) const {
	return GET_BIT(_maskRight, index);
}

int MixedTrace::LeftBit(int index) const {
	return GET_BIT(_maskLeft, index);
}

int MixedTrace::Bit(int index) const {
	if (index >= RightBitNumber()) {
		return LeftBit(index - RightBitNumber());
	}
	else {
		return RightBit(index);
	}
}

int MixedTrace::BitNumber() const {
	return RightBitNumber() + LeftBitNumber();
}

i64 MixedTrace::Hash() const {
	return (((i64)_maskLeft) << 32) | (i64)_maskRight;
}

std::string MixedTrace::ToString() const {
	string ret;

	for (int i = 0; i < RightBitNumber(); i++)
	{
		if (RightBit(i) == 1)
		{
			ret = "X" + ret;
		}
		else
		{
			ret = "Z" + ret;
		}
	}

	for (int i = 0; i < LeftBitNumber(); i++)
	{
		if (LeftBit(i) == 1)
		{
			ret = "x" + ret;
		}
		else
		{
			ret = "z" + ret;
		}
	}


	return ret;
}

ostream& operator<<(ostream& os, const MixedTrace& mt) {
	os << "Tr(" << mt.ToString() << ")";

	return os;
}

/*
void Contract(ContractState& cs) const {
	if (LeftBitNumber() == 0 || RightBitNumber() == 0) {
		return;
	}

	int b = Bit(LeftBitNumber() - 1);
	int A = LeftTrace() - (b << (LeftBitNumber() - 1));
//	SingleTrace A(leftTrace, LeftBitNumber() - 1);
//	SingleTrace right(RightTrace(), RightBitNumber());

	for (int i = 0; i < RightBitNumber(); i++) {
		// for i = 1 it produce single trace of len 1, which vanishes for SU(N) gauge group
		if (i == 1) continue;
		if (Bit(i + LeftBitNumber()) != b) continue;
		// Tr(AzBZC) = Tr(AC)TrB - 1/N * Tr(ABC)
		// do contraction.
		int B = 0;
		if (i > 0) {
			B = RightTrace() & (1 << (i - 1));
		}
		int C = RightTrace() >> (i+1);

		// construct Tr(AC)
		PowerSeries<int> coef1(1);
		int len = LeftBitNumber() + RightBitNumber() - i - 2;
		if (len != 1) {
			if (len == 0) {
				coef1.ShiftOrder(1);
			}
			else {
				;
			}
		}
	}
}*/

MixedMultiTrace::MixedMultiTrace() {
}

MixedMultiTrace::MixedMultiTrace(const MixedMultiTrace& mmt) {
	_left = mmt._left;
	_right = mmt._right;
	if (mmt._mixed != NULL) {
		_mixed = new MixedTrace(*mmt._mixed);
	}
	else {
		_mixed = NULL;
	}
}

MixedMultiTrace::MixedMultiTrace(const StateId& left, const StateId& right, const MixedTrace& mixed) {
	_left = left;
	_right = right;
	_mixed = new MixedTrace(mixed);
}

MixedMultiTrace::MixedMultiTrace(const StateId& left, const StateId& right) {
	_left = left;
	_right = right;
	_mixed = NULL;
}

MixedMultiTrace::MixedMultiTrace(const TraceState& left, const TraceState& right) {
	StateCollection* sc = StateCollection::Inst();
	_left = sc->GetId(left);
	_right = sc->GetId(right);
	_mixed = NULL;
}

MixedMultiTrace::MixedMultiTrace(const TraceState& left, const TraceState& right, const MixedTrace& mixed) {
	StateCollection* sc = StateCollection::Inst();
	_left = sc->GetId(left);
	_right = sc->GetId(right);
	_mixed = new MixedTrace(mixed);
}

MixedMultiTrace::~MixedMultiTrace() {
	if (_mixed != NULL) {
		delete _mixed;
	}
}

pair<i64, i64> MixedMultiTrace::NormHash() const {
	i64 first = ((((i64)_left.BitNumber << 25) | _left.Index) << 32) | ((i64)_right.BitNumber << 25) | _right.Index;
	if (_mixed == NULL) {
		return make_pair(first, (i64)0);
	}
	else {
		return make_pair(first, _mixed->Hash());
	}
}

bool MixedMultiTrace::operator< (const MixedMultiTrace& other) const {
	if (_left < other._left) {
		return true;
	}
	else if (other._left < _left) {
		return false;
	}
	else if (_right < other._right) {
		return true;
	}
	else if (other._right < _right) {
		return false;
	}
	else {
		if ((_mixed == NULL) && (other._mixed == NULL)) {
			return false;
		}
		else if (_mixed == NULL) {
			return true;
		}
		else if (other._mixed == NULL) {
			return false;
		}
		else {
			return *_mixed < *(other._mixed);
		}
	}
}

std::string MixedMultiTrace::ToString() const {
	StateCollection* inst = StateCollection::Inst();
	const TraceState& left = inst->GetState(_left);
	const TraceState& right = inst->GetState(_right);
	string ret = "(" + left.ToString()+ "|";
	if (_mixed != NULL) {
		ret += _mixed->ToString();
	}
	ret += "|" + right.ToString() + ")";
	return ret;
}

ostream& operator<<(ostream& os, const MixedMultiTrace& mt) {
	StateCollection* inst = StateCollection::Inst();
	const TraceState& left = inst->GetState(mt._left);
	const TraceState& right = inst->GetState(mt._right);
	os << "(" << left << "|";
	if (mt._mixed != NULL) {
		os << *mt._mixed;
	}
	os << "|" << right << ")";
	return os;
}

ContractState::ContractState() {
}

void ContractState::AddState(const MixedMultiTrace& mm, const PowerSeries<nType>& coefficient) {
	if (coefficients.find(mm) == coefficients.end())
	{
		coefficients[mm] = coefficient;
	}
	else
	{
		coefficients[mm] += coefficient;
	}
}

void ContractState::Merge(const ContractState& state, const PowerSeries<nType>& prefactor) {
	for (auto it = state.Begin(); it != state.End(); ++it)
	{
		PowerSeries<nType> res = it->second;
		res *= prefactor;
		AddState(it->first, res);
	}
}

void ContractState::Merge(const ContractState& state) {
	for (auto it = state.Begin(); it != state.End(); ++it)
	{
		AddState(it->first, it->second);
	}
}

void ContractState::Extend(const TraceState& left, const TraceState& right) {
	NExpansionSeries coef;
	map<MixedMultiTrace, PowerSeries<nType> > newMix;
	TraceState empty;
	for (auto it = coefficients.begin(); it != coefficients.end(); ++it)
	{
		TraceState res1;
		if (it->first.Left().IsValid()) {
			const TraceState& traceL = StateCollection::Inst()->GetState(it->first.Left());
			TraceState::Merge(left, traceL, res1);
		}
		else {
			TraceState::Merge(left, empty, res1);
		}
		
		TraceState res2;
		if (it->first.Right().IsValid()) {
			const TraceState& traceR = StateCollection::Inst()->GetState(it->first.Right());
			TraceState::Merge(right, traceR, res2);
		}
		else {
			TraceState::Merge(right, empty, res2);
		}
		res1.Normalize(it->second);
		res2.Normalize(it->second);

		StateId idLeft = StateCollection::Inst()->GetId(res1);
		StateId idRight = StateCollection::Inst()->GetId(res2);
		if (!idLeft.IsValid() || !idRight.IsValid()) continue;

		if (it->first.Mixed() != NULL) {
			MixedMultiTrace mt(idLeft, idRight, *it->first.Mixed());
			newMix[mt] = it->second;
		}
		else {
			MixedMultiTrace mt(res1, res2);
			newMix[mt] = it->second;
		}
	}

	this->coefficients = newMix;
}


std::map<MixedMultiTrace, PowerSeries<nType> >::const_iterator ContractState::Begin() const
{
	return coefficients.begin();
}

std::map<MixedMultiTrace, PowerSeries<nType> >::const_iterator ContractState::End() const
{
	return coefficients.end();
}

ostream& operator<<(ostream& os, const ContractState& cs) {
	bool first = true;
	for (auto it = cs.Begin(); it != cs.End(); ++it) {
		if (first) {
			first = false;
		}
		else {
			os << " + ";
		}

		os << "(" <<  it->second << ")" << it->first;
	}
	return os;
}