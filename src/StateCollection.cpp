#pragma warning(disable:4018)
#include "StateCollection.h"

StateCollection* StateCollection::inst = NULL;

StateCollection::StateCollection()
{
	// .. do nothing.
	stateList.resize(MAX_BIT_TO_GENERATE + 1, vector<TraceState>());
	TraceState empty;
	stateList[0].push_back(empty);
	this->state2Id[empty] = StateId(0, 0);
}

StateCollection* StateCollection::Inst()
{
	if (inst == NULL)
	{
		inst = new StateCollection();
	}

	return inst;
}

void StateCollection::Init(int bits, vector<TraceState>& states)
{
	vector<TraceState>& vt = stateList[bits];
	int magnonNumber = 0;
	begins[GetKey(bits, 0)] = 0;
	for (int i = 0; i < states.size(); i++)
	{
	    if (states[i].MagnonNumber() > magnonNumber) {
		    ends[GetKey(bits, magnonNumber)] = vt.size();
		    magnonNumber = states[i].MagnonNumber();
		    begins[GetKey(bits, magnonNumber)] = vt.size();
		}
		vt.push_back(states[i]);
		this->state2Id[states[i]] = StateId(bits, i);
	}
	ends[GetKey(bits, magnonNumber)] = vt.size();
}

StateId StateCollection::GetId(const TraceState& state) const
{
	map<TraceState, StateId>::const_iterator it;
	it = state2Id.find(state); 
	if (it == state2Id.end())
	{
		return StateId(-1, -1);
	}

	return it->second;
}

const TraceState& StateCollection::GetState(const StateId& id) const
{
	return stateList[id.BitNumber][id.Index];
}

int StateCollection::StateNumber(int bits) const
{
	return stateList[bits].size();
}

int StateCollection::StateNumber(int chainLength, int magnonNumber)
{
    return End(chainLength, magnonNumber) - Begin(chainLength, magnonNumber);
}

const TraceState& StateCollection::GetState(int bits, int index) const
{
    return stateList[bits][index];
}

// the begin interator of states with given chain length and magnon number.
vector<TraceState>::const_iterator StateCollection::Begin(int chainLength, int magnonNumber) const
{
    int key = GetKey(chainLength, magnonNumber);
    vector<TraceState>::const_iterator bit = stateList[chainLength].begin();
    vector<TraceState>::const_iterator eit = stateList[chainLength].end();
    auto it = begins.find(key);
    if (it != begins.end()) {
        return bit + it->second;
    } else {
        return eit;
    }
}

// the end interator of states with given chain length and magnon number.
vector<TraceState>::const_iterator StateCollection::End(int chainLength, int magnonNumber) const
{
    int key = GetKey(chainLength, magnonNumber);
    vector<TraceState>::const_iterator bit = stateList[chainLength].begin();
    vector<TraceState>::const_iterator eit = stateList[chainLength].end();
    auto it = ends.find(key);
    if (it != ends.end()) {
        return bit + it->second;
    } else {
        return eit;
    }
}

int StateCollection::GetKey(int L, int M)
{
    return (MAX_BIT_TO_GENERATE + 1) * L + M;
}
