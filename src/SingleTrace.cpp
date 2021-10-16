#pragma warning(disable:4018)
#include <string>
#include <iostream>
#include "BitUtility.h"
#include "SingleTrace.h"

using namespace std;

SingleTrace::SingleTrace()
{
	SingleTrace(0, 0);
}

SingleTrace::SingleTrace(int trace, int bitNumber, bool isAnnihilator)
{
	mask = BuildMask(trace, bitNumber, isAnnihilator);
}

SingleTrace::SingleTrace(string trace)
{
	int n = 0;
	for (int i = 0; i < trace.length(); i++)
	{
		if (trace[i] == 'X')
		{
			// TODO: Why the bits are in reversed order?
//			n |= (1 << i);
			n |= (1 << (trace.length() - i - 1));
		}
	}

	mask = BuildMask(n, trace.length());
}

SingleTrace SingleTrace::Reflect() const {
	int len = this->BitNumber();
	int trace = this->Trace() ^ ((1 << len) - 1);
	SingleTrace ret = SingleTrace(trace, len, this->IsAnnihilator());
	ret.Normalize();
	return ret;
}

int SingleTrace::BitNumber() const
{
	return BIT_LENGTH(mask);
}

int SingleTrace::Bit(int index) const
{
	if ((mask & (1 << index)) != 0)
	{
		return 1;
	}

	return 0;
}

int SingleTrace::Trace() const
{
	return TRACE_BITS(mask);
}

bool SingleTrace::IsAnnihilator() const {
	return ANNIHILATOR_FLAG(mask) != 0;
}

int SingleTrace::MagnonNumber() const
{
	return BitCount(Trace());
}

string SingleTrace::ToString() const
{
	string ret;
	int trace = Trace();
	for (int i = 0; i < BitNumber(); i++)
	{
		if ((trace & (1<<i)) != 0)
		{
			// TODO: Why the order is reversed here?
			// The reason is that we want to the single string ordered in alphabeta
			// by simply order the bits increasely.
			ret = "X" + ret;
//			ret += "X";
		}
		else
		{
			ret = "Z" + ret;
//			ret += "Z";
		}
	}

	return ret;
}

ostream& operator<<(ostream& os, const SingleTrace& single)
{
	os << "Tr(" << single.ToString() << ")";

	return os;
}

string SingleTrace::ToLaTeX() const
{
	string ret;
	int trace = Trace();
	for (int i = 0; i < BitNumber(); i++)
	{
		if ((trace & (1<<i)) != 0)
		{
			ret = "\\bar{X}" + ret;
		}
		else
		{
			ret = "\\bar{Z}" + ret;
		}
	}

	return "\\mathrm{Tr}" + ret;
}

bool SingleTrace::operator< (const SingleTrace& other) const
{
	if (this->BitNumber() != other.BitNumber()) {
		return this->BitNumber() < other.BitNumber();
	}
	return this->Trace() < other.Trace();
}

// a + Bit(pos) + b = original trace
// pos count from right to left: len, len-1, ... pos+1, pos, pos-1, ..., 1, 0.
// a = Tr{len, len-1, ... pos+1}, b = Tr{pos-1, ..., 1, 0}
void SingleTrace::Split(int pos, SingleTrace& a, SingleTrace& b) const
{
	b.mask = BuildMask(PickBits(mask, pos), pos, this->IsAnnihilator());
	int remain = (mask >> (pos + 1));
	a.mask = BuildMask(PickBits(remain, BitNumber() - pos - 1), BitNumber() - pos - 1, this->IsAnnihilator());
}

// a + {pos1} + b + {pos2} + c = original trace
// pos count from right to left and pos1 > pos2: len, len-1, ... pos1,..., pos2,...1, 0.
// a = Tr{len, len-1, ... pos1+1}, b = Tr{pos1-1, ..., pos2+1}, c= {po2-1, ..., 0}
void SingleTrace::Split(int pos1, int pos2, SingleTrace& a, SingleTrace& b, SingleTrace& c) const
{
	SingleTrace d(0, 0);
	Split(pos1, a, d);
	d.Split(pos2, b, c);
}

SingleTrace SingleTrace::Merge(SingleTrace& a, SingleTrace& b)
{
	int trace = (a.Trace() << b.BitNumber()) | b.Trace();
	return SingleTrace(trace, a.BitNumber() + b.BitNumber(), a.IsAnnihilator());	
}

SingleTrace SingleTrace::Merge(SingleTrace& a, SingleTrace& b, SingleTrace& c)
{
	int trace = (a.Trace() << b.BitNumber()) | b.Trace();
	trace = (trace << c.BitNumber()) | c.Trace();
	return SingleTrace(trace, a.BitNumber() + b.BitNumber() + c.BitNumber(), a.IsAnnihilator());
}

SingleTrace SingleTrace::Merge(SingleTrace& a, SingleTrace& b, SingleTrace& c, SingleTrace& d)
{
	int trace = (a.Trace() << b.BitNumber()) | b.Trace();
	trace = (trace << c.BitNumber()) | c.Trace();
	trace = (trace << d.BitNumber()) | d.Trace();
	return SingleTrace(trace, a.BitNumber() + b.BitNumber() + c.BitNumber() + d.BitNumber(), a.IsAnnihilator());
}

void SingleTrace::Normalize()
{
	int trace = this->Trace();
	int n = this->BitNumber();
	int min = trace;
	for (int i = 1; i < n; i++)
	{
		trace = CyclicRotation(trace, n);
		if (trace < min)
		{
			min = trace; 
		}
	}

	this->mask = BuildMask(min, BitNumber(), this->IsAnnihilator());
}
