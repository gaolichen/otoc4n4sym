#pragma warning(disable:4018)
#include "BitUtility.h"
#include <time.h>
using namespace std;

std::ostream* logOut = &std::cout;
LogLevel currentLogLevel = Verbos;

std::ostream& operator<< (std::ostream& out, LogLevel level) {
	if (level == Debug) {
		out << "Debug";
	}
	if (level == Verbos) {
//		out << "Verbos";
	}
	if (level == Info) {
//		out << "Info";
	}
	else if (level == Warning) {
		out << "Warning";
	}
	else if (level == Error) {
		out << "Error";
	}
	else {
		// do nothing.
	}

	return out;
}

int BitCount(int i)
{
	 i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

int CyclicRotation(int n, int bitNumber)
{
	return ((n & 1) << (bitNumber -1)) | (n >> 1);
}

int BitReverse(int n, int bitNumber) {
	int ret = 0;
	while (bitNumber > 0) {
		ret <<= 1;
		if ((n & 1) == 1) {
			ret ^= 1;
		}
		n >>= 1;
		bitNumber--;
	}
	return ret;
}

snum BinomialCoefficient(snum a,i64 b)
{
	if (b > a) return 0;
	if (b + b > a) b = a - b;
	snum ret = 1;
	for (int i = 0; i < b; i++)
	{
		ret *= (a - i);
		ret /= (i+1);
	}

	return ret;
}

int PickBits(int n, int bitNumber)
{
	return n & ((1 << bitNumber) - 1);
}

bool IsBitSet(int n, int bit)
{
	return (n & (1 << bit)) != 0;
}

vector<int> Num2Digit(i64 n, int maxBit)
{
	vector<int> ret;
	ret.reserve(maxBit);
	for (int i = 0; i < maxBit; i++)
	{
		if (IsBitSet(n, i)) ret.push_back(i);
	}

	return ret;
}

i64 Digit2Num(int* v, int size)
{
	i64 ret = 0;
	for (int i = 0; i < size; i++)
	{
		ret |= (1<< v[i]);
	}

	return ret;
}

var_t ToComplex(string str) {
	f_type re, im;
	if (str.find('(') == 0 && str.find(')') == str.length() - 1) {
		int pos = str.find(',');
		if (!ToValue(str.substr(1, pos - 1), re) || !ToValue(str.substr(pos + 1, str.length() - pos - 2), im)) {
			throw std::invalid_argument(str + " is not a valid complex number.");
		}
		else {
			return var_t(re, im);
		}
	}
	else {
		if (!ToValue(str, re)) {
			throw std::invalid_argument(str + " is not a valid complex number.");
		}
		return var_t(re, 0);
	}
}


i64 Factorial(int n)
{
	i64 ret = 1;
	while (n > 0) { ret *= n; n--; }
	return ret;
}

int SymmetryFactor(vector<i64>& mode, int s)
{
	i64 ret = Factorial(s);
	int r = 1;
	for (int j = 1; j < mode.size(); j++)
	{
		if (mode[j] == mode[j-1]) r++;
		else
		{
			ret /= Factorial(r);
			r = 1;
		}
	}

	ret /= Factorial(r);

	return (int)ret;
}

int InverseNumber(const vector<int>& v)
{
	int ret = 0;
	for (int i = 1; i < v.size(); i++)
	{
		for (int j = 0; j < i; j++)
			if (v[i] < v[j]) ret++;
	}

	return ret;
}

int BuildMask(int trace, int bitNumber, bool annihilator)
{
	if (annihilator) {
		return ANNIHILATOR_MASK(trace, bitNumber);
	} else {
		return CREATOR_MASK(trace, bitNumber);
	}
}

int Gcd(int a, int b)
{
	if (b == 0) return a;
	return Gcd(b, a % b);
}

string CombinePath(string path1, string path2)
{
    if (path1.length() > 0) {
	    return path1 + "/" + path2;
	} else {
	    return path2;
	}
}

string Bits2String(int bits, int bitNumber)
{
	string ret(bitNumber, 'Z');
	for (int i = 0; i < bitNumber; i++)
	{
		if ((bits & (1 << i)) != 0)
		{
			ret[bitNumber - 1 - i] = char('X');
		}
	}

	return ret;
}

string ToUpper(string s)
{
	string ret = s;
	for (int i = 0; i < s.length(); i++)
	{
		if (s[i] >= 'a' && s[i] <= 'z')
		{
			ret[i] = char((s[i] - 'a') + 'A');
		}
	}

	return ret;
}

string ToLower(string s)
{
	string ret = s;
	for (int i = 0; i < s.length(); i++)
	{
		if (s[i] >= 'A' && s[i] <= 'Z')
		{
			ret[i] = char((s[i] - 'A') + 'a');
		}
	}

	return ret;
}

f_type Chop(f_type v, f_type eps)
{
    if (std::abs(v) < eps) return 0;
    return v;
}

var_t Chop(var_t v, f_type eps)
{
    return var_t(Chop(v.real(), eps), Chop(v.imag(), eps));
}

var_t ExpI(f_type phase)
{
    return std::polar((f_type)1.0, phase);
}

bool IsZero(f_type value, f_type eps)
{
	return std::abs(value) < eps;
}

Stopwatch::Stopwatch()
{
	start = 0;
}

void Stopwatch::Start()
{
	start = clock();
}

string Stopwatch::Now()
{
    time_t t = time(0);   // get time now
    struct tm * now = localtime(&t);
    char buf[80];
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", now);
    return buf;
}

/// return how much time in second since last Start() function call
double Stopwatch::Stop()
{
	double ret = Elapsed();
	start = 0;
	return ret;
}

double Stopwatch::Elapsed()
{
	return (clock() - start) / (double)CLOCKS_PER_SEC;
}
