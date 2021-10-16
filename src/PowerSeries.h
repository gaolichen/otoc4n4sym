#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include "BitUtility.h"
using namespace std;

template<typename T> class PowerSeries
{
protected:
    T zero = 0;
    vector<pair<int, T> > list;
    
    void Add(const pair<int, T>& item)
    {
        auto it = list.begin();
        while (it != list.end() && it->first < item.first) {
            it++;
        }
        if (it == list.end() || it->first != item.first) {
            list.insert(it, item);
        } else {
            it->second += item.second;
        }
    }
public:
    PowerSeries()
    {
    }

    PowerSeries(T val)
    {
        this->list.push_back(make_pair(0, val));
    }

    PowerSeries(int order, T val)
    {
        this->list.push_back(make_pair(order, val));
    }
    
    PowerSeries(const PowerSeries& ps)
    {
//        cout << "copy constractor ps.list.size()=" << ps.list.size() << endl;
        for (int i = 0; i < ps.list.size(); i++) {
            this->list.push_back(make_pair(ps.list[i].first, ps.list[i].second));
        }
    }

    void Set(int order, T val) {
        int i = 0;
        for (; i < list.size(); i++) {
            if (list[i].first == order) {
                list[i].second = val;
                return;
            }
            else if (list[i].first > order) {
                break;
            }
        }
        list.insert(list.begin() + i, make_pair(order, val));
    }

    int HighestOrder() {
        for (int i = list.size() - 1; i >= 0; i--) {
            if (!(list[i].second == 0)) {
                return list[i].first;
            }
        }
        return 0;
    }

    const T& CoefAt(int pos) const {
        return list[pos].second;
    }

    const int OrderAt(int pos) const {
        return list[pos].first;
    }

    int Size() const {
        return list.size();
    }

    const T& Coef(int order) const
    {
        for (int i = 0; i < list.size(); i++) {
            if (list[i].first == order) {
                return list[i].second;
            } else if (list[i].first > order) {
                return zero;
            }
        }
        return zero;
    }
        
    PowerSeries<T>& operator+= (const PowerSeries<T>& a)
    {
        for (int i = 0; i < a.list.size(); i++) {
            Add(a.list[i]);
        }
        return *this;
    }
	
    PowerSeries<T>& operator*= (const T& n)
	{
	    for (int i = 0; i < list.size(); i++) {
	        list[i].second *= n;
	    }
        return *this;
	}
		
	PowerSeries<T>& operator*= (const PowerSeries<T>& ps)
    {
	    if (this->IsZero()) return *this;
	    if (ps.IsZero())
	    {
		    this->list.resize(0);
		    return *this;
	    }
	    
	    vector<pair<int, T> > res;
	    res.reserve(this->list.size() * ps.list.size());
	    for (int i = 0; i < this->list.size(); i++) {
	        for (int j = 0; j < ps.list.size(); j++) {
	            T prod = this->list[i].second * ps.list[j].second;
	            pair<int, T> p = make_pair(this->list[i].first + ps.list[j].first, prod);
	            res.push_back(p);
	        }
	    }
	    
	    sort(res.begin(), res.end());
	    int j = -1;
	    for (int i = 0; i < res.size(); i++) {
	        if (j == -1 || res[i].first != this->list[j].first) {
	            j++;
                if (this->list.size() > j) {
                    this->list[j] = res[i];
                }
                else {
                    this->list.push_back(res[i]);
                }	            
	        } else {
	            this->list[j].second += res[i].second;
	        }
	    }
	    this->list.resize(j + 1);
	    
	    return *this;
    }
    		
	void ShiftOrder(int order)
	{
	    for (int i = 0; i < list.size(); i++) {
	        list[i].first += order;
	    }
	}
	
	void Clear()
	{
	    list.resize(0);
	}
	
	bool IsZero() const
	{
	    for (int i = 0; i < list.size(); i++) {
	        if (!(list[i].second == 0)) {
	            return false;
	        }
	    }
	    return true;
	}

    string VarToString(int exponent) const
    {
        if (exponent == 1) {
            return "x";
        }
        else if (exponent == -1) {
            return "/x";
        }
        else if (exponent == 0) {
            return "";
        }
        else if (exponent > 1) {
            return "x^" + ToString(exponent);
        }
        else {
            return "/x^" + ToString(-exponent);
        }
    }

    f_type ToNumerical(f_type N) const {
        f_type ret = .0;
        for (int i = 0; i < list.size(); i++) {
            ret += pow(N, list[i].first) * list[i].second;
        }

        return ret;
    }
};

template<typename T>
//ostream& operator<<(ostream& os, const PowerSeries<T>& ps);
ostream& operator<<(ostream& os, const PowerSeries<T>& coef) {
    if (coef.IsZero())
    {
        os << "0";
        return os;
    }

    bool isFirst = true;
    for (int i = 0; i < coef.Size(); i++)
    {
        if (coef.CoefAt(i) == 0) continue;
        T val = coef.CoefAt(i);
        if (isFirst)
        {
            if (val < 0)
            {
                os << "-";
            }

            if (abs(val) != 1 || coef.OrderAt(i) < 1) {
                os << abs(val);
            }

            isFirst = false;
        }
        else
        {
            if (val > 0)
            {
                os << "+";
            }
            else {
                os << "-";
            }
            if (abs(val) != 1 || coef.OrderAt(i) < 1) {
                os << abs(val);
            }
        }

        os << coef.VarToString(coef.OrderAt(i));
    }

    return os;
}

//ostream& operator<<(ostream& os, const PowerSeries<i64>& ps);

class ExpBetaSeries : public PowerSeries<int>
{
private:
	static string ExpToString(int exponent);
public:
    ExpBetaSeries();
    
    ExpBetaSeries(int val);
    
    ExpBetaSeries(int order, int val);
        
    ExpBetaSeries operator+ (const ExpBetaSeries& a) const;
    
    ExpBetaSeries operator* (int n) const;
	
	ExpBetaSeries operator* (const ExpBetaSeries& a) const;
    
    bool operator == (int a) const;
    
    bool operator < (const ExpBetaSeries&) const;
	
	int TermNumber() const;
	
	bool FirstTermPositive() const;
    	    
//    friend ostream& operator<<(ostream& os, const ExpBetaSeries& coef);
    
    // to numerical value. replay variable x with e^(i * beta)
    complex<f_type> ToNumerical(f_type beta) const;
};

class NExpansionSeries : public PowerSeries<ExpBetaSeries>
{
private:
    static string ExpToString(int exponent);
    static void OutputExp(ostream& os, const ExpBetaSeries& val, bool isFirst);
public:
    NExpansionSeries();

    NExpansionSeries(int order, ExpBetaSeries val);

    NExpansionSeries operator+ (const NExpansionSeries& a) const;

	NExpansionSeries operator* (const ExpBetaSeries& n) const;
	
	NExpansionSeries operator* (const NExpansionSeries& val) const;
	
	// to numerical value. replay variable x with e^(i * beta) and N with integer value.
    complex<f_type> ToNumerical(f_type N, f_type beta) const;
    
    friend ostream& operator<<(ostream& os, const NExpansionSeries& coef);
};
