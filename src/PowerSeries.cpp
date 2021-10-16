#include "PowerSeries.h"
using namespace std;


ExpBetaSeries::ExpBetaSeries()
{
    zero = 0;
}

ExpBetaSeries::ExpBetaSeries(int val)
{
    this->list.push_back(make_pair(0, val));
}

ExpBetaSeries::ExpBetaSeries(int order, int val)
{
    this->list.push_back(make_pair(order, val));
}

ExpBetaSeries ExpBetaSeries::operator+ (const ExpBetaSeries& a) const
{
    ExpBetaSeries ret = *this;
    ret += a;
    return ret;
}

ExpBetaSeries ExpBetaSeries::operator* (int n) const
{
    ExpBetaSeries ret = *this;
    ret *= n;
    return ret;
}

ExpBetaSeries ExpBetaSeries::operator* (const ExpBetaSeries& a) const
{
	ExpBetaSeries ret = *this;
    ret *= a;
    return ret;
}

bool ExpBetaSeries::operator == (int a) const
{
    for (int i = 0; i < this->list.size(); i++) {
        if (this->list[i].first != 0) {
            if (this->list[i].second != 0) {
                return false;
            }
        } else {
            if (this->list[i].second != a) {
                return false;
            }
        }
    }

    return true;
}

bool ExpBetaSeries::operator < (const ExpBetaSeries& other) const
{
    if (this->list.size() != other.list.size()) {
        return this->list.size() < other.list.size();
    }
    for (int i = 0; i < this->list.size(); i++) {
        if (this->list[i].first != other.list[i].first) {
            return this->list[i].first < other.list[i].first;
        } else if (this->list[i].second != other.list[i].second) {
            return this->list[i].second < other.list[i].second;
        }
    }
    return false;
}

int ExpBetaSeries::TermNumber() const
{
    int ret = 0;
    for (int i = 0; i < this->list.size(); i++) {
        if (this->list[i].second != 0) {
            ret++;
        }
    }
    return ret;
}

bool ExpBetaSeries::FirstTermPositive() const
{
    for (int i = 0; i < this->list.size(); i++) {
        if (this->list[i].second != 0) {
            return this->list[i].second > 0;
        }
    }
    // it's zero
    return false;
}

string ExpBetaSeries::ExpToString(int exponent)
{
    if (exponent == 1) {
//        return "exp(i*b)";
        return "x";
    } else if (exponent == -1) {
//        return "exp(-i*b)";
        return "/x";
    } else if (exponent == 0) {
        return "";
    } else if (exponent > 1) {
//        return "exp(" + ToString(exponent) + "i*b)";
        return "x^" + ToString(exponent);
    } else {
//        return "exp(-" + ToString(-exponent) + "i*b)";
        return "/x^" + ToString(-exponent);
    }
}

complex<f_type> ExpBetaSeries::ToNumerical(f_type beta) const
{
    complex<f_type> ret = 0.0L;
    for (int i = 0; i < list.size(); i++) {
        ret += ExpI(beta * list[i].first) * (f_type)list[i].second;
    }
    
    return ret;
}
/*
ostream& operator<<(ostream& os, const ExpBetaSeries& coef)
{
    if (coef.IsZero())
    {
	    os << "0";
	    return os;
    }

    bool isFirst = true;
    for (int i = 0; i < coef.list.size(); i++)
    {
	    if (coef.list[i].second == 0) continue;
	    if (isFirst)
	    {
		    os << coef.list[i].second;
		    isFirst = false;
	    }
	    else
	    {
		    if (coef.list[i].second > 0)
		    {
			    os << "+";
		    }
		    os << coef.list[i].second;
	    }
	    
	    os << ExpBetaSeries::ExpToString(coef.list[i].first);
    }

    return os;
}*/

string NExpansionSeries::ExpToString(int exponent)
{
	if (exponent == 1) {
        return "N";
    } else if (exponent == -1) {
        return "/N";
    } else if (exponent == 0) {
        return "";
    } else if (exponent > 1) {
        return "N^" + ToString(exponent);
    } else {
        return "/N^" + ToString(-exponent);
    }
}
    
NExpansionSeries::NExpansionSeries()
{
}

NExpansionSeries::NExpansionSeries(int order, ExpBetaSeries val)
{
//    cout << "NExpansionSeries val=" << val << endl;
    this->list.push_back(make_pair(order, val));
}

NExpansionSeries NExpansionSeries::operator+ (const NExpansionSeries& a) const
{
    NExpansionSeries ret = *this;
    ret += a;
    return ret;
}

NExpansionSeries NExpansionSeries::operator* (const ExpBetaSeries& n) const
{
    NExpansionSeries ret = *this;
    ret *= n;
    return ret;
}

NExpansionSeries NExpansionSeries::operator* (const NExpansionSeries& val) const
{
    NExpansionSeries ret = *this;
    ret *= val;
    return ret;
}

complex<f_type> NExpansionSeries::ToNumerical(f_type N, f_type beta) const
{
    complex<f_type> ret = 0.0L;
    for (int i = 0; i < list.size(); i++) {
        ret += pow((f_type)N, list[i].first) * list[i].second.ToNumerical(beta);
    }
    
    return ret;
}

void NExpansionSeries::OutputExp(ostream& os, const ExpBetaSeries& val, bool isFirst)
{
    if (isFirst == false && (val.TermNumber() > 1 || val.FirstTermPositive())) {
        os << "+";
    }
    if (val.TermNumber() > 1) {
        os << "(";
    }
    os << val;
    if (val.TermNumber() > 1) {
        os << ")";
    }
}

ostream& operator<<(ostream& os, const NExpansionSeries& coef)
{
    if (coef.IsZero())
    {
	    os << "0";
	    return os;
    }

    bool isFirst = true;
    for (int i = 0; i < coef.list.size(); i++)
    {
	    if (coef.list[i].second == 0) continue;
	    NExpansionSeries::OutputExp(os, coef.list[i].second, isFirst);
	    isFirst = false;		    
	    os << NExpansionSeries::ExpToString(coef.list[i].first);
    }

    return os;
}


