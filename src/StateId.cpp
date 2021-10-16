#include"StateId.h"

StateId::StateId()
	: BitNumber(0), Index(0)
{
}

StateId::StateId(int bitNumber, int index)
	: BitNumber(bitNumber), Index(index)
{
	// ...
}

bool StateId::IsValid() const
{
	return BitNumber >= 0 && Index >= 0;
}

bool operator< (const StateId& a, const StateId& b)
{
	if (a.BitNumber != b.BitNumber)
	{
		return a.BitNumber < b.BitNumber;
	}

	return a.Index < b.Index;
}

ostream& operator<<(ostream& os, const StateId& id)
{
//	os << (id.Index / 2) + 1;
    os << id.Index;
	return os;
}
