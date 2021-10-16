#pragma once
#include <iostream>
using namespace std;

class StateId
{
public:
	int BitNumber;
	int Index;
	bool IsValid() const;
	StateId();
	StateId(int bitNumber, int index);
};

bool operator< (const StateId& a, const StateId& b);
ostream& operator<<(ostream& os, const StateId& ms);
