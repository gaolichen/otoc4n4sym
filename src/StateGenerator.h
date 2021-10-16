#pragma once
#include "SingleTrace.h"
#include "TraceState.h"
#include "BitUtility.h"
#include "StateCollection.h"
#include <vector>

using namespace std;

class StateGenerator
{
private:
	int maxBit2Generate;
	bool* myFlags;
	bool visited[30][30][2];

	vector<vector<snum> > stateNumbers;
	vector<snum> singleTraceNumbers;
	vector<vector<SingleTrace> > singles;
	vector<vector<vector<TraceState> > > all;

	void InitSingleTraceNumber();
	void FindSingleStates(int n);
	snum StateNumbers(int bit, int remain);
	void GeneratSingleStates();
	static void DoPickBoson(vector<SingleTrace>& allstates, int index, int remain, vector<SingleTrace>& curr, vector<TraceState>& ret);
	static vector<TraceState> PickBosonFromSingleState(vector<SingleTrace> &states, int number);
	static TraceState CombineStates(TraceState& a, TraceState& b);

//	void BuildSingleOperatorStates(int remBit, int currBits, vector<int>& res);
public:
	StateGenerator(int maxBitToGenerate = MAX_BIT_TO_GENERATE);
	~StateGenerator();
	
	void GenerateAllStates();
	snum StateNumber(int n);
	snum SingleStateNumber(int n);
	SingleTrace& SingleTraceState(int n, int index);
	TraceState State(int n, int index);
	void InitStateCollection(StateCollection* collection);
//	void GenerateSingleOperatorStates(int bits, vector<TraceState> &res);
};
