#pragma warning(disable:4018)
#include "TraceState.h"
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

TraceState::TraceState()
{
}

TraceState::TraceState(vector<SingleTrace>& traces)
{
	vs = traces;
}

TraceState::TraceState(string traces)
{
	int pos, next = -1;
	do
	{
		pos = next + 1;
		next = traces.find_first_of(',', pos);
		if (next != string::npos)
			vs.push_back(SingleTrace(traces.substr(pos, next - pos)));
		else vs.push_back(SingleTrace(traces.substr(pos, next)));
	} while (next != string::npos);
}

TraceState TraceState::Reflect() const {
	TraceState ret;
	
	for (int i = 0; i < vs.size(); i++) {
		ret.AddTrace(vs[i].Reflect());
	}
	ret.QuickNormalize();

	return ret;
}

int TraceState::MagnonNumber() const
{
	int ret = 0;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += vs[i].MagnonNumber();
	}

	return ret;
}

int TraceState::TotalBits() const
{
	int ret = 0;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += vs[i].BitNumber();
	}

	return ret;
}

void TraceState::AddTrace(SingleTrace trace)
{
	vs.push_back(trace);
}

void TraceState::AddTrace(vector<SingleTrace> traces)
{
    vs.insert(vs.end(), traces.begin(), traces.end());
}


vector<SingleTrace>& TraceState::Traces()
{
	return vs;
}

const SingleTrace& TraceState::Trace(int n) const
{
	return vs[n];
}

void TraceState::QuickNormalize()
{
	if (vs.size() <= 1) return;
	sort(vs.begin(), vs.end());
}

bool TraceState::operator< (const TraceState& other) const
{
	int f1 = MagnonNumber();
	int f2 = other.MagnonNumber();
	if (f1 != f2)
	{
		return f1 < f2;
	}

	if (vs.size() != other.vs.size())
	{
		return vs.size() < other.vs.size();
	}

	for (int i = 0; i < vs.size(); i++)
	{
		if (vs[i] < other.vs[i]) return true;
		else if (other.vs[i] < vs[i]) return false;
	}

	return false;
}

int TraceState::TraceNumber() const
{
	return vs.size();
}

string TraceState::ToString() const
{
	string ret;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += "(" + vs[i].ToString() + ")";
	}

	return ret;
}

ostream& operator<<(ostream& os, const TraceState& ts)
{
	for (int i = 0; i < ts.vs.size(); i++)
	{
		os << ts.vs[i];
	}

	return os;
}

string TraceState::ToLaTeX() const
{
	string ret;
	for (int i = 0; i < vs.size(); i++)
	{
		ret += vs[i].ToLaTeX();
	}

	return ret + "\\left|0\\right\\rangle";
}


void TraceState::Merge(const TraceState& a, const TraceState& b, TraceState& res)
{
	res.vs = a.vs;
	res.vs.insert(res.vs.end(), b.vs.begin(), b.vs.end());
}

void TraceState::CopyTo(TraceState& state, int exclude) const
{
	state.vs.insert(state.vs.end(), vs.begin(), vs.begin() + exclude);
	state.vs.insert(state.vs.end(), vs.begin() + exclude + 1, vs.end());
}

void TraceState::CopyTo(TraceState& state, int exclude1, int exclude2) const
{
	state.vs.insert(state.vs.end(), vs.begin(), vs.begin() + exclude1);
	state.vs.insert(state.vs.end(), vs.begin() + exclude1 + 1, vs.begin() + exclude2);
	state.vs.insert(state.vs.end(), vs.begin() + exclude2 + 1, vs.end());
}
