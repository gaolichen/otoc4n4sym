#pragma warning(disable:4018)
#include "Hamiltonian.h"

Hamiltonian::Hamiltonian()
{
	//this->inverted = false;
	//Init(0);
}

//Hamiltonian::Hamiltonian(int xi, bool invert)
//{
//	this->inverted = invert;
//	Init(xi);
//}
//
//void Hamiltonian::Init(int xi)
//{
//	int a = 1;
//	if (inverted)
//	{
//		a = -1;
//	}
//
//	realOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::AA));
//	rePrefactors.push_back(2 * a + 2 * xi);
//	realOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::BB));
//	rePrefactors.push_back(-2 * a + 2 * xi);
//	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::BA));
//	rePrefactors.push_back(2 * a + 2 * xi);
//	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::BA));
//	rePrefactors.push_back(2 * a);
//	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::AB));
//	rePrefactors.push_back(2 * a);
//	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::AB));
//	rePrefactors.push_back(2 * xi - 2 * a);
//
//	realOps.push_back(new BitNumberHamOperator());
//	rePrefactors.push_back(-2 * xi);
//
//	imaginaryOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::AA));
//	imPrefactor.push_back(-2 * a);
//	imaginaryOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::BB));
//	imPrefactor.push_back(2 * a);
//
//	//if (aaaa < 0)
//	//{
//	//	realOps.push_back(new HamOperatorB(HamOperator::AA, HamOperator::AA));
//	//	rePrefactors.push_back(-aaaa);
//	//	realOps.push_back(new HamOperatorB(HamOperator::BA, HamOperator::BA));
//	//	rePrefactors.push_back(-aaaa);
//	//	realOps.push_back(new HamOperatorB(HamOperator::AB, HamOperator::BA));
//	//	rePrefactors.push_back(aaaa);
//	//}
//}

Hamiltonian::~Hamiltonian()
{
	for (int i = 0; i < realOps.size(); i++)
	{
		delete realOps[i];
	}

	for (int i = 0; i < imaginaryOps.size(); i++)
	{
		delete imaginaryOps[i];
	}
}

void Hamiltonian::AddReadOp(HamOperator *op, int prefactor)
{
	for (int i = 0; i < realOps.size(); i++)
	{
		if (*realOps[i] == *op)
		{
			rePrefactors[i] += prefactor;
			return;
		}
	}

	this->realOps.push_back(op);
	this->rePrefactors.push_back(prefactor);
}

void Hamiltonian::AddImaginaryOp(HamOperator *op, int prefactor)
{
	for (int i = 0; i < realOps.size(); i++)
	{
		if (*imaginaryOps[i] == *op)
		{
			imPrefactor[i] += prefactor;
			return;
		}
	}

	this->imaginaryOps.push_back(op);
	this->imPrefactor.push_back(prefactor);
}

void Hamiltonian::Apply(const TraceState& state, MixState& real, MixState& imaginary)
{
	for (int i = 0; i < realOps.size(); i++)
	{
		MixState res;
		realOps[i]->ApplyOn(state, res);
		real.Merge(res, rePrefactors[i]);
	}

	for (int i = 0; i < imaginaryOps.size(); i++)
	{
		MixState res;
		imaginaryOps[i]->ApplyOn(state, res);
		imaginary.Merge(res, imPrefactor[i]);
	}
}

void Hamiltonian::Matrix(int bits, StateType type, vector<vector<Coefficient> >& rem, vector<vector<Coefficient> >& imm)
{
	int n = StateCollection::Inst()->StateNumber(bits);
	map<StateId, Coefficient>::const_iterator it;
	rem.resize(n, vector<Coefficient>());
	imm.resize(n, vector<Coefficient>());
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			rem[i].push_back(Coefficient(0, 0));
			imm[i].push_back(Coefficient(0, 0));
		}
	}

	for (int i = 0; i < n; i++)
	{
		const TraceState& state = StateCollection::Inst()->GetState(bits, i, type);
		MixState re, im;
		Apply(state, re, im);
		for (it = re.Begin(); it != re.End(); ++it)
		{
			rem[it->first.Index / 2][i] = it->second;
		}

		for (it = im.Begin(); it != im.End(); ++it)
		{
			imm[it->first.Index / 2][i] = it->second;
		}
	}
}

ostream& operator<<(ostream& os, const Hamiltonian& ham)
{
	bool isFirst = true;
	for (int i = 0; i < ham.realOps.size(); i++)
	{
		if (ham.rePrefactors[i] == 0) continue;
		if (!isFirst && ham.rePrefactors[i] > 0)
		{
			os << " + ";
		}
		os << ham.rePrefactors[i] << *ham.realOps[i];
		isFirst = false;
	}

	os << " + i*{";

	for (int i = 0; i < ham.imaginaryOps.size(); i++)
	{
		if (ham.imPrefactor[i] == 0) continue;
		if (!isFirst && ham.imPrefactor[i] > 0)
		{
			os << " + ";
		}
		os << ham.imPrefactor[i] << *ham.realOps[i];
		isFirst = false;
	}

	os << "}";

	return os;
}

HamOperator* Hamiltonian::RealOp(int index)
{
	return realOps[index];
}

int Hamiltonian::RealOpSize()
{
	return realOps.size();
}

H0Hamiltonian::H0Hamiltonian()
{
	filePrefix = "h0Ham";
	realOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::AA));
	rePrefactors.push_back(2);
	realOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::BB));
	rePrefactors.push_back(-2);
	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::BA));
	rePrefactors.push_back(2);
	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::BA));
	rePrefactors.push_back(2);
	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::AB));
	rePrefactors.push_back(2);
	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::AB));
	rePrefactors.push_back(-2);

	imaginaryOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::AA));
	imPrefactor.push_back(-2);
	imaginaryOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::BB));
	imPrefactor.push_back(2);
}

DeltaHamiltonian::DeltaHamiltonian()
{
	filePrefix = "deltaHam";

	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::BA));
	rePrefactors.push_back(2);

	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::AB));
	rePrefactors.push_back(2);

	realOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::AA));
	rePrefactors.push_back(2);

	realOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::BB));
	rePrefactors.push_back(2);

	realOps.push_back(new BitNumberHamOperator());
	rePrefactors.push_back(-2);
}

ZeroHamiltonian::ZeroHamiltonian()
{
	filePrefix = "zeroHam";
	realOps.push_back(new HamOperatorA(HamOperator::AA, HamOperator::AA));
	rePrefactors.push_back(1);
	realOps.push_back(new HamOperatorA(HamOperator::BB, HamOperator::BB));
	rePrefactors.push_back(1);
	realOps.push_back(new HamOperatorA(HamOperator::AB, HamOperator::BA));
	rePrefactors.push_back(1);
	realOps.push_back(new HamOperatorA(HamOperator::BA, HamOperator::AB));
	rePrefactors.push_back(1);
	realOps.push_back(new BitNumberHamOperator());
	rePrefactors.push_back(-1);
	realOps.push_back(new HamOperatorB(HamOperator::AA, HamOperator::AA));
	rePrefactors.push_back(-1);
	realOps.push_back(new HamOperatorB(HamOperator::BA, HamOperator::BA));
	rePrefactors.push_back(-1);
	realOps.push_back(new HamOperatorB(HamOperator::AB, HamOperator::BA));
	rePrefactors.push_back(1);
}
