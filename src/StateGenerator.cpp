#pragma warning(disable:4018)
#include "StateGenerator.h"
#include "BitUtility.h"
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;
//#define TEST_STATENUMBER

StateGenerator::StateGenerator(int maxBitToGenerate)
{
	this->maxBit2Generate = maxBitToGenerate;
	vector<snum> numbers = vector<snum>(MAX_BIT_TO_COUNT + 1, (snum)(-1));
	stateNumbers.resize(MAX_BIT_TO_COUNT + 1, numbers);
	myFlags = new bool[1 << maxBit2Generate];
	InitSingleTraceNumber();
}

StateGenerator::~StateGenerator()
{
	delete[] myFlags;
}

void StateGenerator::InitSingleTraceNumber()
{
	singleTraceNumbers.resize(MAX_BIT_TO_COUNT + 1, 0);
	for (int m = 2; m <= MAX_BIT_TO_COUNT; m++)
	{
#ifdef TEST_STATENUMBER
		singleTraceNumbers[m] = ((((i64)1) << (m - 1)) + m - 1) / m;
#else
		snum res = 0;
		for (int i = 1; i <= m; i++)
		{
			int toShift = Gcd(i, m);
			res += ((i64)1 << toShift);
		}

		singleTraceNumbers[m] = res / m;
#endif
	}
}

void StateGenerator::GeneratSingleStates()
{
	singles.resize(maxBit2Generate + 1);
	for (int i = 2; i <= maxBit2Generate; i++)
	{
		FindSingleStates(i);
	}
}

void StateGenerator::FindSingleStates(int bit)
{
	memset(myFlags, 0, sizeof(bool) * (1<<bit));

	for (int i = 0; i < (1<<bit); i++)
	{
		if (myFlags[i]) continue;
		int n = i;
		for (int j = 1; j < bit; j++)
		{
			n = CyclicRotation(n, bit);
			if (n == i)
			{
			    break;
			}
			else myFlags[n] = true;
		}

		singles[bit].push_back(SingleTrace(i, bit));
	}
}

snum StateGenerator::SingleStateNumber(int n)
{
	return singleTraceNumbers[n];
}

snum StateGenerator::StateNumber(int n)
{
	return StateNumbers(n, n);
}

snum StateGenerator::StateNumbers(int bit, int remain)
{
	if (bit == 0)
	{
		if (remain == 0) return 1;
		else return 0;
	}
	if (remain == 0)
	{
		return 1;
	}

	snum ret = stateNumbers[bit][remain];
	if (ret >= 0) return ret;
	int a = remain / bit;

	ret = 0;
	for (int i = 0; i <= a; i++)
	{
		//int n1 = BinomialCoefficient(singleFermions[bit].size(), i);
		snum n2 = BinomialCoefficient(singleTraceNumbers[bit] - 1 + i, (i64)i);
		ret += n2 * StateNumbers(bit - 1, remain - i * bit);
	}
	stateNumbers[bit][remain] = ret;

	return ret;
}

void StateGenerator::GenerateAllStates()
{
	GeneratSingleStates();

	memset(visited, 0, sizeof(visited));
	all.resize(maxBit2Generate + 1, vector<vector<TraceState> >(maxBit2Generate + 1, vector<TraceState>()));
	all[0][0] = vector<TraceState>(1, TraceState());

	vector<vector<TraceState> > bv;

	for (int i = 1; i <= maxBit2Generate; i++)
	{
		// find TraceStates for i-bit.
		bv.clear();

		// maximum number of i-bit states can be in one state.
		int n = maxBit2Generate / i;
		// bv[k]: all states cosnsiting of k number of i-bit single trace.
		for (int k = 0; k <= n; k++)
		{
			bv.push_back(PickBosonFromSingleState(singles[i], k));
		}

		for (int j = 0; j <= maxBit2Generate; j++)
		{
			int a = j / i;
			for (int h = 0; h <= a; h++)
			{
				for (int y = 0; y < bv[h].size(); y++)
				{
					vector<TraceState>& vts2 = all[i - 1][j - h * i];
					for (int z = 0; z < vts2.size(); z++)
					{
						all[i][j].push_back(CombineStates(vts2[z], bv[h][y]));
					}
				}
			}
		}
	}

	for (int i = 1; i <= maxBit2Generate; i++)
	{
		for (int j = 0; j < all[i][i].size(); j++)
		{
			all[i][i][j].QuickNormalize();
		}

		sort(all[i][i].begin(), all[i][i].end());
	}
}

void StateGenerator::DoPickBoson(vector<SingleTrace>& allstates, int index, int remain,
	vector<SingleTrace>& curr, vector<TraceState>& ret)
{
	if (remain == 0)
	{
		ret.push_back(TraceState(curr));
		return;
	}
	for (int i = index; i < allstates.size(); i++) {
		curr.push_back(allstates[i]);
		DoPickBoson(allstates, i, remain - 1, curr, ret);
		curr.pop_back();
	}

/*	if (index == allstates.size())
	{
		if (remain == 0)
		{
			ret.push_back(TraceState(curr));
		}
		else
		{
			//.. do nothing.
		}

		return;
	}

	for (int i = 0; i <= remain; i++)
	{
		DoPickBoson(allstates, index + 1, remain - i, curr, ret);
		curr.push_back(allstates[index]);
	}

	for (int i = 0; i <= remain; i++)
	{
		curr.pop_back();
	}*/
}

vector<TraceState> StateGenerator::PickBosonFromSingleState(vector<SingleTrace> &states, int number)
{
	vector<TraceState> ret;
	vector<SingleTrace> curr;
	DoPickBoson(states, 0, number, curr, ret);

	return ret;
}

TraceState StateGenerator::CombineStates(TraceState& a, TraceState& b)
{
	TraceState ret;
	for (int i = 0; i < a.TraceNumber(); i++)
	{
		ret.AddTrace(a.Traces()[i]);
	}

	for (int i = 0; i < b.TraceNumber(); i++)
	{
		ret.AddTrace(b.Traces()[i]);
	}

	return ret;
}

SingleTrace& StateGenerator::SingleTraceState(int n, int index) {
    return singles[n][index];
}

TraceState StateGenerator::State(int n, int index)
{
	return all[n][n][index];
}

void StateGenerator::InitStateCollection(StateCollection* collection)
{
	for (int i = 1; i <= maxBit2Generate; i++)
	{
		collection->Init(i, all[i][i]);
	}
}

