#pragma once
#include "TraceState.h"
#include "PowerSeries.h"
#include "MixedTrace.h"
#include <vector>
#include <map>
using namespace std;

class NormCalculator
{
private:
protected:
	PowerSeries<nType> zero;
	virtual PowerSeries<nType> DoCalculate(const TraceState&, const TraceState&) = 0;
public:
	NormCalculator();
	virtual PowerSeries<nType> Calculate(const TraceState&, const TraceState&);
};

class BruteForceCalculator : public NormCalculator
{
private:
	bool visited[1<<5];
	map< pair<TraceState, TraceState>, PowerSeries<nType> >::iterator cIt;
	map< pair<TraceState, TraceState>, PowerSeries<nType> > res;
	static void Positions(const TraceState&, vector<int>& bPos, vector<int>& fPos);
	static void IndexPosMap(const TraceState&, vector<int>& index2Pos, vector<int>& pos2Index);
protected:
	virtual PowerSeries<nType> DoCalculate(const TraceState&, const TraceState&);
public:
	BruteForceCalculator();
};

class RecursiveNormCalculatorBase : public NormCalculator
{
protected:
	static ITrace* CreateTrace(int left, int lLen, int right, int rLen);
	//	void AddStates(ITrace* trace, PowerSeries<int>& prefactor, ContractState& res);
	void AddStates(ITrace* trace, SingleTrace* single, PowerSeries<nType>& prefactor, ContractState& res);
	ContractState Contract(const SingleTrace& left, const SingleTrace& right);
	ContractState Contract(const MixedTrace& mixed, const SingleTrace& right);
	ContractState Contract(const MixedTrace& mixed);
};

class RecursiveNormCalculator : public RecursiveNormCalculatorBase
{
private:
	map< pair<TraceState, TraceState>, PowerSeries<nType> > _cache;
	void DoContract(const MixedMultiTrace& mt, const PowerSeries<nType>& factor, ContractState& res);
protected:
	virtual PowerSeries<nType> DoCalculate(const TraceState&, const TraceState&);
public:
	RecursiveNormCalculator();
};

class DfsNormCalculator : public RecursiveNormCalculatorBase
{
private:
	PowerSeries<nType> _one;
//	map<MixedMultiTrace, PowerSeries<nType> > _cache;
	map<pair<i64, i64>, PowerSeries<nType> > _cache;
protected:
	virtual PowerSeries<nType> DoCalculate(const TraceState&, const TraceState&);
	PowerSeries<nType>& Dfs(const MixedMultiTrace& mmt);
public:
	DfsNormCalculator();
};

