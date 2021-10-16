#pragma warning(disable:4018)
#include "NormCalculator.h"
#include "MixedTrace.h"
#include "StateCollection.h"
#include <vector>
#include <algorithm>
#include <cstring>
using namespace std;

NormCalculator::NormCalculator()
{
	// ..
}

PowerSeries<nType> NormCalculator::Calculate(const TraceState& a, const TraceState& b)
{
	if (a.MagnonNumber() != b.MagnonNumber() || a.TotalBits() != b.TotalBits()) {
		return zero;
	}

	// exchange X and Z if magnonNumber is greater than 1/2 of total bits.
	if (a.MagnonNumber() + a.MagnonNumber() > a.TotalBits()) {
		TraceState aa = a.Reflect();
		TraceState bb = b.Reflect();
		if (aa < bb) {
			return DoCalculate(aa, bb);
		}
		else {
			return DoCalculate(bb, aa);
		}
	}
	else {
		if (a < b) {
			return DoCalculate(a, b);
		}
		else {
			return DoCalculate(b, a);
		}
	}
}

BruteForceCalculator::BruteForceCalculator()
{
	// ..
}

void BruteForceCalculator::Positions(const TraceState& state, vector<int>& bPos, vector<int>& fPos)
{
	int p = 0;
	for (int i = 0; i < state.TraceNumber(); i++)
	{
		const SingleTrace& single = state.Trace(i);
		for (int j = 0; j < single.BitNumber(); j++)
		{
			if (single.Bit(single.BitNumber() - j - 1) == 0)
			{
				bPos.push_back(p + j);
			}
			else
			{
				fPos.push_back(p + j);
			}
		}

		p += single.BitNumber();
	}
}

void BruteForceCalculator::IndexPosMap(const TraceState& state, vector<int>& index2Pos, vector<int>& pos2Index)
{
	int p = 0;
	for (int i = 0; i < state.TraceNumber(); i++)
	{
		int bits = state.Trace(i).BitNumber();
		for (int j = 0; j < bits; j++)
		{
			index2Pos[p + j] = p + (j + 1) % bits;
			pos2Index[p + (j + 1) % bits] = j + p;
		}

		p += bits;
	}
}

/// Calculate <0| A B |0> where A, B are two tracestates.
/// It bruteforcely go through all the possible contractions of A and B.
/// For each contraction, find how many pow of N it produce and parity of it (even or odd)
PowerSeries<nType> BruteForceCalculator::DoCalculate(const TraceState& a, const TraceState& b)
{
	vector<int> bPos1, bPos2, fPos1, fPos2;
	Positions(a, bPos1, fPos1);
	Positions(b, bPos2, fPos2);
	pair<TraceState, TraceState> key = make_pair(a, b);
	cIt = res.find(key);
	if (cIt != res.end())
	{
		return cIt->second;
	}

	PowerSeries<nType> ret;
	int bits = a.TotalBits();
	vector<int> pos2Tindex(bits);
	vector<int> dIndex2Pos(bits);
	vector<int> temp(bits);
	IndexPosMap(a, dIndex2Pos, temp);
	IndexPosMap(b, temp, pos2Tindex);
	vector<int> creator2Annihilator(bits);
	vector<int> annihilator2Creator(bits);

	sort(fPos1.begin(), fPos1.end());
	do
	{
		for (int i = 0; i < fPos1.size(); i++)
		{
			creator2Annihilator[fPos2[i]] = fPos1[i];
			annihilator2Creator[fPos1[i]] = fPos2[i];
		}

		sort(bPos1.begin(), bPos1.end());
		do
		{
			for (int i = 0; i < bPos1.size(); i++)
			{
				creator2Annihilator[bPos2[i]] = bPos1[i];
				annihilator2Creator[bPos1[i]] = bPos2[i];
			}

			memset(visited, 0x0, sizeof(visited));

			// find how many loops in this contraction configuration
			int loop = 0;
			for (int dIndex = 0; dIndex < bits; dIndex++)
			{
				if (visited[dIndex]) continue;
				loop++;
				int i = dIndex;
				int len = 0;
				do{
					visited[i] = true;
					len++;
					int j = creator2Annihilator[i];
					int k = dIndex2Pos[j];
					int x = annihilator2Creator[k];
					i = pos2Tindex[x];
				} while (!visited[i]);

				// For SU(N) case, for loop of length 1, it should vanish.
//				if (len == 1)
//				{
//					hasTrivialLoop = true;
//				}
			}
			ret += PowerSeries<nType>(loop, 1);
		} while(next_permutation(bPos1.begin(), bPos1.end()));
	} while(next_permutation(fPos1.begin(), fPos1.end()));

	res[key] = ret;

	return ret;
}

RecursiveNormCalculator::RecursiveNormCalculator() {
}

PowerSeries<nType> RecursiveNormCalculator::DoCalculate(const TraceState& left, const TraceState& right) {
	pair<TraceState, TraceState> key = make_pair(left, right);
	auto cIt = _cache.find(key);
	if (cIt != _cache.end())
	{
		return cIt->second;
	}

	StateCollection* sc = StateCollection::Inst();
	MixedMultiTrace mt(left, right);
	ContractState res1;
	ContractState res2;
	ContractState& prev = res1;
	ContractState& curr = res2;
	PowerSeries<nType> val(1);
	prev.AddState(mt, val);
	for (int i = 0; i < left.TotalBits(); i++) {
		for (auto it = prev.Begin(); it != prev.End(); it++) {
			DoContract(it->first, it->second, curr);
		}
//		std::cout << "i=" << i << ":" << curr << std::endl;
		std::swap(prev, curr);
		curr.Clear();
	}
//	return prev.Begin()->second;
	if (prev.Size() == 0) {
		PowerSeries<nType> res;
		_cache[key] = res;
	} else {
		PowerSeries<nType> res = prev.Begin()->second;
		_cache[key] = res;
	}
	
	return _cache[key];
}

void RecursiveNormCalculator::DoContract(const MixedMultiTrace& mt, const PowerSeries<nType>& factor, ContractState& res) {
	StateCollection* sc = StateCollection::Inst();
	const TraceState& right = sc->GetState(mt.Right());
	const TraceState& left = sc->GetState(mt.Left());

	if (mt.Mixed() != NULL) {
		if (mt.Mixed()->Contractable()) {
			ContractState state = Contract(*mt.Mixed());
			state.Extend(left, right);
			res.Merge(state, factor);
		}

		// copy all all traces to toAppend except the i-th trace.
		for (int i = 0; i < right.TraceNumber(); i++) {
			TraceState toAppend;
			right.CopyTo(toAppend, i);
			ContractState state = Contract(*mt.Mixed(), right.Trace(i));
			state.Extend(left, toAppend);
			res.Merge(state, factor);
		}
	}
	else {
		TraceState toAppendLeft;
		SingleTrace single = left.Trace(0);
		left.CopyTo(toAppendLeft, 0);
		for (int i = 0; i < right.TraceNumber(); i++) {
			TraceState toAppend;
			right.CopyTo(toAppend, i);
			ContractState state = Contract(single, right.Trace(i));
			state.Extend(toAppendLeft, toAppend);
			res.Merge(state, factor);
		}
	}
}

ContractState RecursiveNormCalculatorBase::Contract(const MixedTrace& mixed) {
	ContractState ret;
	if (mixed.LeftBitNumber() == 0 || mixed.RightBitNumber() == 0) {
		return ret;
	}

//	int b = mixed.Bit(mixed.LeftBitNumber() - 1);
	int b = mixed.LeftBit(0);
//	int A = mixed.LeftTrace() - (b << (mixed.LeftBitNumber() - 1));
	int A = (mixed.LeftTrace() >> 1);

	for (int i = 0; i < mixed.RightBitNumber(); i++) {
		if (mixed.RightBit(i) != b) continue;
		// Tr(AzBZC) = Tr(AC)TrB - 1/N * Tr(ABC)
		// do contraction.
//		int	B = mixed.RightTrace() & ((1 << i) - 1);
//		int C = mixed.RightTrace() >> (i + 1);
		int	B = mixed.RightTrace() >> (i + 1);
		int C = mixed.RightTrace() & ((1 << i) - 1);

		// construct Tr(AC)
//		ITrace* ac = CreateTrace(A, mixed.LeftBitNumber() - 1, C, mixed.RightBitNumber() - i - 1);
		ITrace* ac = CreateTrace(A, mixed.LeftBitNumber() - 1, C, i);
		// construct Tr(B)
//		SingleTrace* b = new SingleTrace(B, i);
		SingleTrace* b = new SingleTrace(B, mixed.RightBitNumber() - i - 1);
//		int BC = B + (C << i);
		int BC = C + (B << i);
		ITrace* abc = CreateTrace(A, mixed.LeftBitNumber() - 1, BC, mixed.RightBitNumber() - 1);
		PowerSeries<nType> coef1(1);
		AddStates(ac, b, coef1, ret);

		PowerSeries<nType> coef2(-1);
		coef2.ShiftOrder(-1);
		AddStates(abc, NULL, coef2, ret);

		delete ac;
		delete b;
		delete abc;
	}
	return ret;
}

ITrace* RecursiveNormCalculatorBase::CreateTrace(int left, int lLen, int right, int rLen) {
	if (lLen == 0) {
		return new SingleTrace(right, rLen, false);
	}
	else if (rLen == 0) {
		return new SingleTrace(BitReverse(left, lLen), lLen, true);
	}
	else {
		return new MixedTrace(left, lLen, right, rLen);
	}
}

ContractState RecursiveNormCalculatorBase::Contract(const SingleTrace& left, const SingleTrace& right) {
	// note that for left trace, bits are ordered (from left to right) as high bit to low bit
	// for right trace, bits are ordered (from left to right) as low bit to high bit
	// Tr(az)Tr(ZB) = Tr(aB) - 1/N Tr(a)Tr(B)
	ContractState ret;
	if (left.BitNumber() == 0) {
		return ret;
	}

//	int b = left.Bit(0);
	int b = left.Bit(left.BitNumber() - 1);
	// we need to reverse bits for A
//	int A = BitReverse(left.Trace() >> 1, left.BitNumber() - 1);
	int A = BitReverse(PickBits(left.Trace(), left.BitNumber() - 1), left.BitNumber() - 1);
	int rTrace = right.Trace();

	for (int i = 0; i < right.BitNumber(); i++) {
		rTrace = CyclicRotation(rTrace, right.BitNumber());
		if (right.Bit(i) != b) continue;
		// Tr(az)Tr(ZB) = Tr(aB) - 1/N Tr(a)Tr(B)
		// do contraction.
		int B = PickBits(rTrace, right.BitNumber() - 1);

		// construct Tr(a)
//		SingleTrace* a = new SingleTrace(left.Trace() >> 1, left.BitNumber() - 1, true);
		SingleTrace* a = new SingleTrace(PickBits(left.Trace(), left.BitNumber() - 1), left.BitNumber() - 1, true);
		// construct Tr(B)
		SingleTrace* b = new SingleTrace(B, right.BitNumber() - 1);

		// B cannot be empty, so we need to use reversed A here.
		// construct Tr(aB)
		ITrace* ab = CreateTrace(A, left.BitNumber() - 1, B, right.BitNumber() - 1);
		PowerSeries<nType> coef1(1);
		AddStates(ab, NULL, coef1, ret);

		PowerSeries<nType> coef2(-1);
		coef2.ShiftOrder(-1);
		AddStates(a, b, coef2, ret);

		delete ab;
		delete a;
		delete b;
	}

//	std::cout << "ret=" << ret << std::endl;

	return ret;
}

ContractState RecursiveNormCalculatorBase::Contract(const MixedTrace& mixed, const SingleTrace& right) {
	ContractState ret;
	if (mixed.LeftBitNumber() == 0) {
		return ret;
	}

//	int b = mixed.Bit(mixed.LeftBitNumber() - 1);
	int b = mixed.LeftBit(0);
//	int A = mixed.LeftTrace() - (b << (mixed.LeftBitNumber() - 1));
	int A = mixed.LeftTrace() >> 1;
	int rTrace = right.Trace();

	for (int i = 0; i < right.BitNumber(); i++) {
		rTrace = CyclicRotation(rTrace, right.BitNumber());
		if (right.Bit(i) != b) continue;
		// Tr(azB)Tr(ZC) = Tr(aCB) - 1/N * Tr(aB)Tr(C)
		// do contraction.
		int	B = mixed.RightTrace();
		int C = PickBits(rTrace, right.BitNumber() - 1);
//		int CB = C | (B << (right.BitNumber() - 1));
		int CB = B | (C << (mixed.RightBitNumber()));

		// construct Tr(aCB)
		ITrace* acb = CreateTrace(A, mixed.LeftBitNumber() - 1, CB, mixed.RightBitNumber() + right.BitNumber() - 1);
		// construct Tr(C)
		SingleTrace* c = new SingleTrace(C, right.BitNumber() - 1);
		// construct Tr(aB)
		ITrace* ab = CreateTrace(A, mixed.LeftBitNumber() - 1, B, mixed.RightBitNumber());
		PowerSeries<nType> coef1(1);
		AddStates(acb, NULL, coef1, ret);

		PowerSeries<nType> coef2(-1);
		coef2.ShiftOrder(-1);
		AddStates(ab, c, coef2, ret);

		delete acb;
		delete ab;
		delete c;
	}
	return ret;
}

void RecursiveNormCalculatorBase::AddStates(ITrace* trace, SingleTrace* single, PowerSeries<nType>& coef, ContractState& res) {
	TraceState ts;
	if (single != NULL) {
		ts.AddTrace(*single);
	}
	ts.Normalize(coef);
	if (coef.IsZero()) {
		return;
	}

	StateId emptyId(0, 0);
	StateCollection* inst = StateCollection::Inst();

	MixedTrace* mt = dynamic_cast<MixedTrace*>(trace);
	if (mt != NULL) {
		StateId rightId = inst->GetId(ts);
		MixedMultiTrace mmt(emptyId, rightId, *mt);
		res.AddState(mmt, coef);
	}
	else {
		// both are SingleTrace
		SingleTrace* st = dynamic_cast<SingleTrace*>(trace);
		if (st->IsAnnihilator()) {
			TraceState left;
			left.AddTrace(*st);
			left.Normalize(coef);
			if (coef.IsZero()) {
				return;
			}
//			std::cout << "left=" << left << std::endl;
			MixedMultiTrace mmt(left, ts);
			res.AddState(mmt, coef);
		}
		else {
			ts.AddTrace(*st);
			ts.Normalize(coef);
			if (coef.IsZero()) {
				return;
			}
			StateId rightId = inst->GetId(ts);
			MixedMultiTrace mmt(emptyId, rightId);
			res.AddState(mmt, coef);
		}
	}
}
/*
void RecursiveNormCalculator::AddStates(ITrace* trace, PowerSeries<int>& prefactor, ContractState& res) {
	MixedTrace* mt = dynamic_cast<MixedTrace*>(trace);
	StateId invalidId(-1, -1);
	if (mt != NULL) {
		MixedMultiTrace mmt(invalidId, invalidId, *mt);
		res.AddState(mmt, prefactor);
	}
	else {
		SingleTrace* st = dynamic_cast<SingleTrace*>(trace);
		TraceState ts;
		ts.AddTrace(*st);
		ts.Normalize(prefactor);
		if (prefactor.IsZero()) {
			return;
		}
		StateCollection* inst = StateCollection::Inst();
		StateId id = inst->GetId(ts);
		if (st->IsAnnihilator()) {
			MixedMultiTrace mmt(id, invalidId);
			res.AddState(mmt, prefactor);
		}
		else {
			MixedMultiTrace mmt(invalidId, id);
			res.AddState(mmt, prefactor);
		}
	}
}*/


DfsNormCalculator::DfsNormCalculator() {
	_one.Set(0, 1);
}

PowerSeries<nType> DfsNormCalculator::DoCalculate(const TraceState& left, const TraceState& right) {
	MixedMultiTrace mmt(left, right);
	return Dfs(mmt);
}

PowerSeries<nType>& DfsNormCalculator::Dfs(const MixedMultiTrace& mt){
	if (mt.Mixed() == NULL && mt.Left().BitNumber == 0 && mt.Right().BitNumber == 0) {
		return _one;
	}

	pair<i64, i64> key = mt.NormHash();
	auto cIt = _cache.find(key);
	if (cIt != _cache.end())
	{
		return cIt->second;
	}

	PowerSeries<nType>& ret = _cache[key];

	StateCollection* sc = StateCollection::Inst();
	const TraceState& right = sc->GetState(mt.Right());
	const TraceState& left = sc->GetState(mt.Left());
	ContractState res;
	if (mt.Mixed() != NULL) {
		if (mt.Mixed()->Contractable()) {
			ContractState state = Contract(*mt.Mixed());
			state.Extend(left, right);
			res.Merge(state);
		}

		// copy all all traces to toAppend except the i-th trace.
		for (int i = 0; i < right.TraceNumber(); i++) {
			TraceState toAppend;
			right.CopyTo(toAppend, i);
			ContractState state = Contract(*mt.Mixed(), right.Trace(i));
			state.Extend(left, toAppend);
			res.Merge(state);
		}
	}
	else if (left.TraceNumber() > 0) {
		TraceState toAppendLeft;
		SingleTrace single = left.Trace(0);
		left.CopyTo(toAppendLeft, 0);
		for (int i = 0; i < right.TraceNumber(); i++) {
			TraceState toAppend;
			right.CopyTo(toAppend, i);
			ContractState state = Contract(single, right.Trace(i));
			state.Extend(toAppendLeft, toAppend);
			res.Merge(state);
		}
	}

	for (auto it = res.Begin(); it != res.End(); ++it) {
		PowerSeries<nType> tmp = it->second;
		tmp *= Dfs(it->first);
		ret += tmp;
	}

	return ret;
}