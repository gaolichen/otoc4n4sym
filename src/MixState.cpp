#include "MixState.h"
#include <sstream>
using namespace std;

MixState::MixState()
{
}

void MixState::AddState(const StateId& stateId, const NExpansionSeries& coefficient)
{
	if (coefficients.find(stateId) == coefficients.end())
	{
		coefficients[stateId] = coefficient;
	}
	else
	{
		coefficients[stateId] += coefficient;
	}
}

void MixState::Merge(const MixState& state, int prefactor)
{
	map<StateId, NExpansionSeries>::const_iterator it;
	ExpBetaSeries eb = prefactor;
	for (it = state.coefficients.begin(); it != state.coefficients.end(); ++it)
	{
		if (prefactor != 1)
		{
			AddState(it->first, it->second * eb);
		}
		else
		{
			AddState(it->first, it->second);
		}
	}
}

void MixState::Merge(const MixState& state)
{
	Merge(state, 1);
}

MixState& MixState::operator*= (int n)
{
	map<StateId, NExpansionSeries>::iterator it;
	for (it = coefficients.begin(); it != coefficients.end(); ++it)
	{
		it->second *= n;
	}

	return *this;
}

ostream& operator<<(ostream& os, const MixState& ms)
{
	for (map<StateId, NExpansionSeries>::const_iterator it = ms.coefficients.begin(); it != ms.coefficients.end(); ++it)
	{
		os << it->second << "|" << it->first << ">";
		os << ",";
	}

	return os;
}

string MixState::ToString()
{
	ostringstream oss;
	bool isFirst = true;
	for (map<StateId, NExpansionSeries>::const_iterator it = coefficients.begin(); it != coefficients.end(); ++it)
	{
		if (it->second.IsZero()) continue;
		if (!isFirst)
		{
			oss << "+";
		}
		oss << it->second  << StateCollection::Inst()->GetState(it->first);
		isFirst = false;
	}
	
	// if all terms are zero.
	if (isFirst)
	{
		oss << "0";
	}

	return oss.str();
}

void MixState::Extend(TraceState& state, int parity)
{
	NExpansionSeries coef;
	map<StateId, NExpansionSeries> newMix;
	for (map<StateId, NExpansionSeries>::iterator it = coefficients.begin(); it != coefficients.end(); ++it)
	{
		const TraceState& ts = StateCollection::Inst()->GetState(it->first);
		TraceState res;
		TraceState::Merge(ts, state, res);
		// Here, coef can only be either 1 or -1.
		res.Normalize(it->second);
		StateId id = StateCollection::Inst()->GetId(res);
		if (!id.IsValid()) continue;

        newMix[id] = it->second;
	}

	this->coefficients = newMix;
}

map<StateId, NExpansionSeries>::const_iterator MixState::Begin() const
{
	return coefficients.begin();
}

map<StateId, NExpansionSeries>::const_iterator MixState::End() const
{
	return coefficients.end();
}
