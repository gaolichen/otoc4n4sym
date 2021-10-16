#pragma once
#include <string>
#include <vector>
#include "BitUtility.h"
using namespace std;

class ITrace
{
public:
	virtual int Bit(int index) const = 0;
	virtual int BitNumber() const = 0;
};

class SingleTrace : public ITrace
{
private:
	int mask;
public:
	SingleTrace();
	SingleTrace(string trace);
	SingleTrace(int trace, int bitNumber, bool isAnnihilator = false);

	// exchange X and Z to obtain a new trace
	SingleTrace Reflect() const;

	virtual int BitNumber() const;
	virtual int Bit(int index) const;
	bool IsAnnihilator() const;
	int Trace() const;
	int MagnonNumber() const;

	// normalize the trace to it's minimum form
	void Normalize();
	void Split(int pos, SingleTrace& a, SingleTrace& b) const;
	void Split(int pos1, int pos2, SingleTrace& a, SingleTrace& b, SingleTrace& c) const;
	string ToString() const;
	string ToLaTeX() const;

	bool operator< (const SingleTrace& other) const;
	friend ostream& operator<<(ostream& os, const SingleTrace& single);

	static SingleTrace Merge(SingleTrace& a, SingleTrace& b);
	static SingleTrace Merge(SingleTrace& a, SingleTrace& b, SingleTrace& c);
	static SingleTrace Merge(SingleTrace& a, SingleTrace& b, SingleTrace& c, SingleTrace& d);
};
