#pragma once
#include <vector>
#include <iostream>
#include <string>
#include "SingleTrace.h"
//#include "Coefficient.h"
#include "PowerSeries.h"
using namespace std;

class TraceState
{
private:
	vector<SingleTrace> vs;
public:
	TraceState();
	TraceState(string);
	TraceState(vector<SingleTrace> &);

	// exchange X and Z
	TraceState Reflect() const;

	int TraceNumber() const;
	int MagnonNumber() const;
	int TotalBits() const;
	void AddTrace(SingleTrace trace);
	void AddTrace(vector<SingleTrace> trace);
	void QuickNormalize();
//	void Normalize(NExpansionSeries& coef);
	vector<SingleTrace>& Traces();
	const SingleTrace& Trace(int n) const;
	string ToString() const;
	string ToLaTeX() const;
	void CopyTo(TraceState& state, int exclude) const;
	void CopyTo(TraceState& state, int exclude1, int exclude2) const;
	bool operator< (const TraceState&) const;
	friend ostream& operator<<(ostream& os, const TraceState& ms);
	static void Merge(const TraceState& a, const TraceState& b, TraceState& res);

	template<class T> void Normalize(PowerSeries<T>& coef)
	{
		for (int i = vs.size() - 1; i >= 0; i--)
		{
			if (vs[i].BitNumber() == 0)
			{
				coef.ShiftOrder(1);
				vs.erase(vs.begin() + i);
			}
			else if (vs[i].BitNumber() == 1)
			{
				coef.Clear();
				return;
			}
			else {
				vs[i].Normalize();
			}
		}

		for (int i = 0; i < vs.size(); i++)
		{
			int m = i;
			for (int j = i + 1; j < vs.size(); j++)
			{
				if (vs[j] < vs[m])
				{
					m = j;
				}
			}

			if (m != i)
			{
				// swap vs[i], vs[m]
				swap(vs[i], vs[m]);
			}
		}
	}

};
