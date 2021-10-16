#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <sstream>
#include <complex>
#include <stdexcept>
#include <regex>
#include <time.h>
#include <Eigen/Dense>
using namespace std;

//#define HAM_PARAMETER
#define EPS 1e-10

#ifdef _WIN32
typedef __int64 i64;
#else
typedef long long i64;
#endif

//typedef double snum;
typedef i64 snum;

// float type.
//typedef long double f_type;
typedef double f_type;
typedef std::complex<f_type> var_t;
typedef Eigen::Matrix<var_t, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Array<var_t, Eigen::Dynamic, 1> Array;
typedef Eigen::Matrix<var_t, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorI;

#define MAX_TRACE_BITS					25
#define FULL_TRACE_BITS					((1 << MAX_TRACE_BITS)-1)
#define TRACE_BITS(n)					((n) & FULL_TRACE_BITS)
#define BIT_LENGTH(n)					(((n) >> MAX_TRACE_BITS) & 31)
#define ANNIHILATOR_FLAG(n)				(((n) >> 30) & 1)
#define CREATOR_MASK(trace, len)		((trace) | ((len) << MAX_TRACE_BITS))
#define ANNIHILATOR_MASK(trace, len)	((trace) | (((len) | 32) << MAX_TRACE_BITS))
#define DEBUG 0
#if DEBUG
#define MAX_BIT_TO_COUNT				30
#define MAX_BIT_TO_GENERATE				10
#else
#define MAX_BIT_TO_COUNT				61
#define MAX_BIT_TO_GENERATE				20
#endif

enum LogLevel { Debug = 0, Verbos = 1, Info = 2, Warning = 3, Error = 4, Off = 100 };

std::ostream& operator<< (std::ostream& out, LogLevel level);

extern std::ostream* logOut;
extern LogLevel currentLogLevel;

#define LOG(message, level) \
if (logOut && level >= currentLogLevel) { \
    (*logOut) << Stopwatch::Now() << ' ' << level << ' ' << message << std::endl; \
}

#define SET_LOG_LEVEL(level) \
currentLogLevel = level; \
if (level != Off) { \
    logOut = &std::cout; \
} else { \
    logOut = NULL; \
}



int BitCount(int);
int BitReverse(int n, int bitNumber);
int CyclicRotation(int, int);
snum BinomialCoefficient(snum a, i64 b);
int PickBits(int n, int bitNumber);
bool IsBitSet(int n, int bit);
int InverseNumber(const vector<int>& v);
int BuildMask(int trace, int bitNumber, bool annihilator = false);
int Gcd(int a, int b);
string CombinePath(string path1, string path2);
string Bits2String(int bits, int bitNumber);
string ToUpper(string s);
string ToLower(string s);
vector<int> Num2Digit(i64 n, int maxBit);
i64 Digit2Num(int* v, int size);
int SymmetryFactor(vector<i64>& mode, int s);
f_type Chop(f_type v, f_type eps = EPS);
var_t Chop(var_t v, f_type eps = EPS);
var_t ExpI(f_type phase);
bool IsZero(f_type value, f_type eps = EPS);

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  out << '[';
  if (!v.empty()) {
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b";
  }
  //out << "\b\b]";
  out << ']';
  return out;
}

template<class T> string ToString(T a)
{
	ostringstream oss;
	oss << a;
	return oss.str();
};

template<class T> bool ToValue(string str, T& res)
{
	T a;
	istringstream iss(str);
	if (iss >> a) {
		res = a;
		return true;
	}
	else {
		return false;
	}
};

var_t ToComplex(string str);

template<class T> T ToExpression(string str, T defaultValue)
{
    T a;
	istringstream iss(str);
	if (iss >> a) {
	    return a;
	} else {
	    return defaultValue;
	}
};


class Stopwatch
{
private:
	clock_t start;
public:
	Stopwatch();
	static string Now();
	void Start();
	double Stop();
	double Elapsed();
};


class OtocException : public std::exception
{
private:
	std::string _message;
public:
	OtocException(std::string message) {
		this->_message = message;
	}

	virtual const char* what() const noexcept {
		return _message.c_str();
	}
};
