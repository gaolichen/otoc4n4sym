#pragma once
#include <map>
#include <vector>
#include "TraceState.h"
#include "StateId.h"
using namespace std;

class StateCollection
{
private:
	map<TraceState, StateId> state2Id;
	vector<vector<TraceState> > stateList;
	map<int, int> begins;
	map<int, int> ends;
	StateCollection();
	static int GetKey(int L, int M);
	static StateCollection* inst;
public:
	static StateCollection* Inst();
	void Init(int bits, vector<TraceState>& states);
	StateId GetId(const TraceState& state) const;
	const TraceState& GetState(const StateId& id) const;
	int StateNumber(int bits) const;
	int StateNumber(int chainLength, int magnonNumber);
	const TraceState& GetState(int bits, int index) const;
	
	// the begin interator of states with given chain length and magnon number.
	vector<TraceState>::const_iterator Begin(int chainLength, int magnonNumber) const;
	
	// the begin interator of states with given chain length and magnon number.
	vector<TraceState>::const_iterator End(int chainLength, int magnonNumber) const;
};
