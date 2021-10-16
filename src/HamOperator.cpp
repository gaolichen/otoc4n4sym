#include "HamOperator.h"
#include "SingleTrace.h"
#include <iostream>
#include "StateId.h"

using namespace std;

f_type inverseTempeture = 0.0;
f_type beta_nsym = 0.0;

void TraceNumberOperator::ApplyOn(const TraceState& state, MixState& res) {
	NExpansionSeries coef(0, state.TraceNumber());
	res.AddState(StateCollection::Inst()->GetId(state), coef);
}

void MagnonTraceNumberOperator::ApplyOn(const TraceState& state, MixState& res) {
	int n = 0;
	for (int i = 0; i < state.TraceNumber(); i++) {
		if (state.Trace(i).MagnonNumber() > 0) {
			n++;
		}
	}
	NExpansionSeries coef(0, n);
	res.AddState(StateCollection::Inst()->GetId(state), coef);
}

void NormalOrderedOperator::ToMatrixN(int L, f_type N, BlockMatrix& res) {
	_op->ToMatrixN(L, N, res);
	BlockMatrix ham(L);
	DilatationOperator op;
	for (int m = 0; m <= L; m++) {
		ham.fill(m, m);
		op.ToMatrixN(L, m, N, beta_nsym, ham.block(m, m));
		ham.block(m, m) *= -inverseTempeture;
	}
	BlockMatrix expBH = ham.exp();
	var_t vev = (expBH * res).trace() / expBH.trace();
	int size = StateCollection::Inst()->StateNumber(L);
	res -= BlockMatrix(L, Array::Constant(size, vev));
}

void NormalOrderedOperator::ToMatrixN(int L, int M, f_type N, Matrix& res) {
	_op->ToMatrixN(L, M, N, res);
	DilatationOperator op;
	int size = StateCollection::Inst()->StateNumber(L, M);
	Matrix ham(size, size);
	ham.fill(0.0);
	op.ToMatrixN(L, M, N, beta_nsym, ham);
	Matrix expBH = (ham * (-inverseTempeture)).exp();
	var_t vev = (expBH * res).trace() / expBH.trace();
	res -= Array::Constant(res.rows(), vev).matrix().asDiagonal();
}

TraceOperator::TraceOperator()
{
	// ..
}

void TraceOperator::ApplyOn(const TraceState& state, MixState& res)
{
	for (int i = 0; i < state.TraceNumber(); i++)
	{
		TraceState toAppend;
		// copy all all traces to toAppend except the i-th trace.
		state.CopyTo(toAppend, i);
		MixState ms = ApplyOnPrivate(state.Trace(i));
//		cout << "toAppend=" << toAppend << endl;
//		cout << "mixstate=" << ms << endl;
		ms.Extend(toAppend, 0);
//		cout << "single trace after extend=" << ms << endl;
		res.Merge(ms);
		for (int j = i + 1; j < state.TraceNumber(); j++)
		{
			TraceState toAdd;
		    // copy all all traces to toAdd except the i-th and j-th traces.
			state.CopyTo(toAdd, i, j);
			MixState ms2 = ApplyOnPrivate(state.Trace(i), state.Trace(j));
			ms2.Extend(toAdd, 0);
//		    cout << "double trace after extend=" << ms << endl;
			res.Merge(ms2);
		}
	}
}

MixState& TraceOperator::ApplyOnPrivate(const SingleTrace& single)
{
	map<SingleTrace, MixState>::iterator it = cache1.find(single);
	if (it != cache1.end())
	{
		return it->second;
	}

	MixState res;
	ApplyOnSingle(single, res);

	cache1[single] = res;
	return cache1[single];
}

MixState& TraceOperator::ApplyOnPrivate(const SingleTrace& single1, const SingleTrace& single2)
{
	pair<SingleTrace, SingleTrace> key = make_pair(single1, single2);

	map<pair<SingleTrace, SingleTrace>, MixState>::iterator it = cache2.find(key);
	if (it != cache2.end())
	{
		return it->second;
	}

	MixState res;
	ApplyOnTwoSingle(single1, single2, res);

	cache2[key] = res;
	return cache2[key];
}

void TraceOperator::AddState(TraceState& state, NExpansionSeries& coef, MixState& res)
{
//    cout << "TraceOperator::AddState before state=" << state << ", coef=" << coef << endl;
	state.Normalize(coef);
//	cout << "TraceOperator::AddState normalized state=" << state << ", coef=" << coef << endl;
	if (coef.IsZero())
	{
	    return;
	}
//	cout << "TraceOperator::AddState nonzero" << endl;

	StateId id = StateCollection::Inst()->GetId(state);
//	cout << "TraceOperator::AddState state id=" << id << endl;
	if (!id.IsValid())
	{
		return;
	}
//	cout << "after state=" << state << ", id=" << id << endl;

	res.AddState(id, coef);
}

void DilatationOperator::ToMatrix(int L, int M, vector<vector<NExpansionSeries> >& res)
{
    StateCollection* inst = StateCollection::Inst();
    auto beginIt = inst->Begin(L, M);
    int beginIndex = inst->GetId(*beginIt).Index;
	int n = inst->StateNumber(L, M);
	map<StateId, NExpansionSeries>::const_iterator it;
	res.resize(n, vector<NExpansionSeries>(n));
	for (auto it = inst->Begin(L, M); it != inst->End(L, M); it++)
	{
		MixState ms;
		ApplyOn(*it, ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			res[it2->first.Index - beginIndex][it - beginIt] = it2->second;
		}
	}
}

void DilatationOperator::ToMatrixN(int L, int M, f_type N, f_type beta, Matrix& res) {
	StateCollection* inst = StateCollection::Inst();
	auto beginIt = inst->Begin(L, M);
	int beginIndex = inst->GetId(*beginIt).Index;
	int n = inst->StateNumber(L, M);
	for (auto it = inst->Begin(L, M); it != inst->End(L, M); it++)
	{
		MixState ms;
		ApplyOn(*it, ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			res(it2->first.Index - beginIndex, it - beginIt) = it2->second.ToNumerical(N, beta);
		}
	}
}

void DilatationOperator::ToMatrix(int L, vector<vector<NExpansionSeries> >& res)
{
	StateCollection* inst = StateCollection::Inst();
	int n = inst->StateNumber(L);
	res.resize(n, vector<NExpansionSeries>(n));
	for (int i = 0; i < n; i++)
	{
		MixState ms;
		ApplyOn(inst->GetState(L, i), ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			res[it2->first.Index][i] = it2->second;
		}
	}
}


ostream& operator<<(ostream& os, const TraceOperator& op)
{
	os << op.ToString();
	return os;
}

bool operator == (const TraceOperator& op1, const TraceOperator& op2)
{
	return op1.ToString() == op2.ToString();
}

DilatationOperator::DilatationOperator()
{
}

void DilatationOperator::ApplyOnSingle(const SingleTrace& single, MixState& res)
{
	int n = single.BitNumber();
	SingleTrace xz("XZ");
	SingleTrace zx("ZX");
	SingleTrace A, B, s1, s2;
//	cout << "ApplyOnSingle, single=" << single << endl;
	for (int i = 0; i < n; i++)
	{
	    if (single.Bit(i) == 1) continue;
	    // find a "Z"
		for (int jj = i + 1; jj < n + i; jj++)
		{
		    int j = jj % n;
		    if (single.Bit(j) == 0) continue;
		    // find a "X"
			if (j > i)
			{
				single.Split(j, i, s2, A, s1);
				B = SingleTrace::Merge(s1, s2);
			} else {
				single.Split(i, j, s2, B, s1);
				A = SingleTrace::Merge(s1, s2);
			}
			
			// first term
			// 2/N * e^(-ib)Tr(A)Tr([X,Z]_{b}B)
			vector<SingleTrace> traces1 {A,  SingleTrace::Merge(xz, B)};
			vector<SingleTrace> traces2 {A,  SingleTrace::Merge(zx, B)};
			NExpansionSeries coef1(-1, 2); // 2/N
			NExpansionSeries coef2(-1, ExpBetaSeries(-2, -2)); // -2*e^(-2ib)/N
			TraceState state1(traces1);
			TraceState state2(traces2);
			AddState(state1, coef1, res);
			AddState(state2, coef2, res);
			
			// second term
			// -2/N * e^(ib)Tr([X,Z]_{b}A)Tr(B)
		    vector<SingleTrace> traces3 {SingleTrace::Merge(xz, A), B};
			vector<SingleTrace> traces4 {SingleTrace::Merge(zx, A), B};
			NExpansionSeries coef3(-1, ExpBetaSeries(2, -2)); // -2*e^(2ib)/N
			NExpansionSeries coef4(-1, 2); // 2/N
			
			TraceState state3(traces3);
			TraceState state4(traces4);
			AddState(state3, coef3, res);
			AddState(state4, coef4, res);
			
			
			// third term
			// 2(e^(ib) - e^(-ib))/N^2 * Tr([X,Z]_{b}{A,B})
			vector<SingleTrace> traces5 {SingleTrace::Merge(xz, A, B)};
			vector<SingleTrace> traces6 {SingleTrace::Merge(xz, B, A)};
			// 2(e^(2ib) - 1)/N^2
			NExpansionSeries coef5(-2, ExpBetaSeries(2, 2) + ExpBetaSeries(0, -2));
			
			TraceState state5(traces5);
			TraceState state6(traces6);
			AddState(state5, coef5, res);
			AddState(state6, coef5, res);
			
		    vector<SingleTrace> traces7 {SingleTrace::Merge(zx, A, B)};
			vector<SingleTrace> traces8 {SingleTrace::Merge(zx, B, A)};
			// 2(e^(-2ib) - 1)/N^2
			NExpansionSeries coef7(-2, ExpBetaSeries(-2, 2) + ExpBetaSeries(0, -2));
			
			TraceState state7(traces7);
			TraceState state8(traces8);
			AddState(state7, coef7, res);
			AddState(state8, coef7, res);

            // fourth term
			// 2(e^(ib) - e^(-ib))/N^2 * Tr([X,Z]_{b})Tr(A)Tr(B)
			vector<SingleTrace> traces9 {xz, A, B};
			vector<SingleTrace> traces10 {zx, A, B};
			// 2(e^(2ib) - 1)/N^2
			NExpansionSeries coef9(-2, ExpBetaSeries(2, 2) + ExpBetaSeries(0, -2));
			// 2(e^(-2ib)-1)
			NExpansionSeries coef10(-2, ExpBetaSeries(-2, 2) + ExpBetaSeries(0, -2));
			
			TraceState state9(traces9);
			TraceState state10(traces10);
			AddState(state9, coef9, res);
			AddState(state10, coef10, res);
			
			// fifth term
			// -4(e^(ib) - e^(-ib))/N^3 * Tr([X,Z]_{b})Tr(AB)
			vector<SingleTrace> traces11 {xz, SingleTrace::Merge(A, B)};
			vector<SingleTrace> traces12 {zx, SingleTrace::Merge(A, B)};
			// -4(e^(2ib) - 1)/N^3
			NExpansionSeries coef11(-3, ExpBetaSeries(2, -4) + ExpBetaSeries(0, 4));
			// 4(1 - e^(-2ib))/N^3
			NExpansionSeries coef12(-3, ExpBetaSeries(0, 4) + ExpBetaSeries(-2, -4));
			
			TraceState state11(traces11);
			TraceState state12(traces12);
			AddState(state11, coef11, res);
			AddState(state12, coef12, res);
		}
	}
}

void DilatationOperator::ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res)
{
//    cout << "ApplyOnTwoSingle" << endl;
    int n1 = single1.BitNumber();
    int n2 = single2.BitNumber();
    SingleTrace A, B;
    SingleTrace xz("XZ");
    SingleTrace zx("ZX");
    SingleTrace s1, s2;
    for (int i = 0; i < n1; i++) {
        single1.Split(i, s1, s2);
        if (single1.Bit(i) == 1) {
            // find a "X"
            A = SingleTrace::Merge(s2, s1);
        } else {
            // find a "Z"
            B = SingleTrace::Merge(s2, s1);
        }
        
        for (int j = 0; j < n2; j++) {
            if (single1.Bit(i) == single2.Bit(j)) continue;
            single2.Split(j, s1, s2);
            if (single2.Bit(j) == 1) {
                A = SingleTrace::Merge(s2, s1);
            } else {
                B = SingleTrace::Merge(s2, s1);
            }
        
            // first term
            // 2/N * e^(-ib) Tr([X,Z]_{b}BA)
            vector<SingleTrace> traces1 {SingleTrace::Merge(xz, B, A)};
		    vector<SingleTrace> traces2 {SingleTrace::Merge(zx, B, A)};
		    NExpansionSeries coef1(-1, 2); // 2/N
		    NExpansionSeries coef2(-1, ExpBetaSeries(-2, -2)); // -2*e^(-2ib)/N
		    
		    TraceState state1(traces1);
		    TraceState state2(traces2);
		    AddState(state1, coef1, res);
		    AddState(state2, coef2, res);

            
            // second term
            // -2/N * e^(ib) Tr([X,Z]_{b}AB)
            vector<SingleTrace> traces3 {SingleTrace::Merge(xz, A, B)};
		    vector<SingleTrace> traces4 {SingleTrace::Merge(zx, A, B)};
		    NExpansionSeries coef3(-1, ExpBetaSeries(2, -2)); // -2/N * e^(2ib)
		    NExpansionSeries coef4(-1, 2); // 2/N
		    
		    TraceState state3(traces3);
		    TraceState state4(traces4);
		    AddState(state3, coef3, res);
		    AddState(state4, coef4, res);
            
            // third term
            // 2(e^(ib) - e^(-ib))/N^2 * Tr(A)Tr([X,Z]_{b}B)
            vector<SingleTrace> traces5 {A, SingleTrace::Merge(xz, B)};
		    vector<SingleTrace> traces6 {A, SingleTrace::Merge(zx, B)};
		    NExpansionSeries coef5(-2, ExpBetaSeries(2, 2) + ExpBetaSeries(-2)); // 2(e^(2ib) - 1)/N^2
		    NExpansionSeries coef6(-2, ExpBetaSeries(-2, 2) + ExpBetaSeries(-2)); // 2(-1+e^(-2ib))/N^2
		    
		    TraceState state5(traces5);
		    TraceState state6(traces6);
		    AddState(state5, coef5, res);
		    AddState(state6, coef6, res);
            
            // fourth term
            // 2(e^(ib) - e^(-ib))/N^2 * Tr([X,Z]_{b}A)Tr(B)
            vector<SingleTrace> traces5a {SingleTrace::Merge(xz, A), B};
		    vector<SingleTrace> traces6a {SingleTrace::Merge(zx, A), B};
		    NExpansionSeries coef5a(-2, ExpBetaSeries(2, 2) + ExpBetaSeries(-2)); // 2(e^(2ib) - 1)/N^2
		    NExpansionSeries coef6a(-2, ExpBetaSeries(-2, 2) + ExpBetaSeries(-2)); // 2(-1+e^(-2ib))/N^2

		    
		    TraceState state5a(traces5a);
		    TraceState state6a(traces6a);
		    AddState(state5a, coef5a, res);
		    AddState(state6a, coef6a, res);
            
            // fifth term
            // 2(e^(ib) - e^(-ib))/N^2 * Tr([X,Z]_{b})Tr(AB)
            vector<SingleTrace> traces5b {xz, SingleTrace::Merge(A, B)};
		    vector<SingleTrace> traces6b {zx, SingleTrace::Merge(A, B)};
		    NExpansionSeries coef5b(-2, ExpBetaSeries(2, 2) + ExpBetaSeries(-2)); // 2(e^(2ib) - 1)/N^2
		    NExpansionSeries coef6b(-2, ExpBetaSeries(-2, 2) + ExpBetaSeries(-2)); // 2(-1+e^(-2ib))/N^2

		
		    TraceState state5b(traces5b);
		    TraceState state6b(traces6b);
		    AddState(state5b, coef5b, res);
		    AddState(state6b, coef6b, res);
            
            // sixth term
            // -4(e^(ib) - e^(-ib))/N^3 * Tr([X,Z]_{b})Tr(A)Tr(B)
            vector<SingleTrace> traces7 {xz, A, B};
		    vector<SingleTrace> traces8 {zx, A, B};
		    NExpansionSeries coef7(-3, ExpBetaSeries(2, -4) + 4); // -4(e^(2ib) - 1)/N^3
		    NExpansionSeries coef8(-3, ExpBetaSeries(-2, -4) + 4); // 4(1-e^(-2ib))/N^3
		    
		    TraceState state7(traces7);
		    TraceState state8(traces8);
		    AddState(state7, coef7, res);
		    AddState(state8, coef8, res);
		}
    }
}

TwoFieldOperator::TwoFieldOperator(int creator, int annihilator) {
	_creator = creator & 1;
	_annihilator = annihilator & 1;
}

TwoFieldOperator::TwoFieldOperator(std::string op) {
	_creator = 0;
	_annihilator = 0;
	if (op.length() != 2 ||
		std::string("XZ").find(op[0]) == std::string::npos ||
		std::string("XZ").find(op[1]) == std::string::npos) {
		throw OtocException("Argument Error: op=" + op + ": must be 2-character string consisting of either X or Z.");
	}

	if (op[0] == 'X') {
		_creator = 1;
	}
	if (op[1] == 'X') {
		_annihilator = 1;
	}
}

void TwoFieldOperator::ApplyOnSingle(const SingleTrace& single, MixState& res) {
	int n = single.BitNumber();
	SingleTrace b, c;
	SingleTrace a(this->_creator & 1, 1);
	int a1 = (this->_annihilator & 1);
	for (int i = 0; i < n; i++)
	{
		if (single.Bit(i) == a1)
		{
			NExpansionSeries coef(0, ExpBetaSeries(1));
			single.Split(i, b, c);
			TraceState state;
			state.AddTrace(SingleTrace::Merge(a, c, b));
			AddState(state, coef, res);

			NExpansionSeries coef2(-1, ExpBetaSeries(-1));
			TraceState state2;
			state2.AddTrace(SingleTrace::Merge(c, b));
			state2.AddTrace(a);
			AddState(state2, coef2, res);
		}
	}
}

FourFieldOperator::FourFieldOperator(int creator, int annihilator) {
	this->_creator = creator;
	this->_annihilator = annihilator;
}

int OpStringToInteger(std::string op) {
	int ret = 0;
	for (int i = 0; i < op.length(); i++) {
		ret <<= 1;
		if (op[i] == 'X') {
			ret += 1;
		}
	}

	return ret;
}

FourFieldOperator::FourFieldOperator(std::string op) {
	if (op.length() != 4 ||
		std::string("XZ").find(op[0]) == std::string::npos ||
		std::string("XZ").find(op[1]) == std::string::npos ||
		std::string("XZ").find(op[2]) == std::string::npos ||
		std::string("XZ").find(op[3]) == std::string::npos) {
		throw OtocException("Argument Error: op=" + op + ": must be 4-character string consisting of either X or Z.");
	}

	int n = OpStringToInteger(op);
	this->_creator = (n>>2) & 3;
	this->_annihilator = n & 3;
}

string FourFieldOperator::ToString() const
{
	int n = (this->_creator << 2) | this->_annihilator;
	return "Tr(" + Bits2String(n, 4) + ")";
}
/*
void ISimpleOperator::ToMatrixN(int L, f_type N, Matrix& res) {
	StateCollection* inst = StateCollection::Inst();
	int n = inst->StateNumber(L);
	for (int i = 0; i < n; i++)
	{
		MixState ms;
		this->ApplyOn(inst->GetState(L, i), ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			res(it2->first.Index, i) = it2->second.ToNumerical(N, .0L);
		}
	}
}*/

void ISimpleOperator::ToMatrixN(int L, f_type N, BlockMatrix& res) {
	res.init(L);
	for (int m1 = 0; m1 <= L; m1++) {
		int m2 = m1 + MagnonNumberChange();
		if (m2 < 0 || m2 > L) {
			continue;
		}
		ToMatrixN(L, m1, m2, N, res);
	}
}

void ISimpleOperator::ToMatrixN(int L, int M1, int M2, f_type N, BlockMatrix& res) {
	StateCollection* inst = StateCollection::Inst();
	auto bIt1 = inst->Begin(L, M1);
	auto bIt2 = inst->Begin(L, M2);
	int beginIndex = inst->GetId(*bIt2).Index;

	res.fill(M2, M1);
	for (auto it = inst->Begin(L, M1); it != inst->End(L, M1); it++)
	{
		MixState ms;
		ApplyOn(*it, ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			(res.block(M2, M1))(it2->first.Index - beginIndex, it - bIt1) = it2->second.ToNumerical(N, .0);
		}
	}
}

void ISimpleOperator::ToMatrixN(int L, int M, f_type N, Matrix& res) {
	if (MagnonNumberChange() != 0) {
		throw OtocException("The operator " + this->ToString() + " does not preserve magnon number.");
	}

	StateCollection* inst = StateCollection::Inst();
	res.resize(inst->StateNumber(L, M), inst->StateNumber(L, M));
	res.fill(0.0);
	auto beginIt = inst->Begin(L, M);
	int beginIndex = inst->GetId(*beginIt).Index;
	int n = inst->StateNumber(L, M);
	for (auto it = inst->Begin(L, M); it != inst->End(L, M); it++)
	{
		MixState ms;
		ApplyOn(*it, ms);
		for (auto it2 = ms.Begin(); it2 != ms.End(); ++it2)
		{
			res(it2->first.Index - beginIndex, it - beginIt) = it2->second.ToNumerical(N, .0);
		}
	}
}

void FourFieldOperator::ApplyOnSingle(const SingleTrace& single, MixState& res) {
	int n = single.BitNumber();
	SingleTrace b, c, d;
	SingleTrace a(this->_creator, 2);
	int a2 = (this->_annihilator & 1);
	int a1 = this->_annihilator / 2;
	// position of a1 and a2: Tr(...a1a2)

	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			// Tr(...j...i...)
			if (single.Bit(i) == a1 && single.Bit(j) == a2)
			{
				// i == a1, j == a2
				// 1/N*Tr(Aa1a2) Tr(BjCiD) = 1/N Tr(Aa1a2) Tr(jCiDB)
				// = 1/N*Tr(ADB)Tr(C)-1/N^2*(Tr(ACDB) + Tr(ADBC)) + 1/N^3*Tr(A)Tr(CDB)
				single.Split(j, i, b, c, d);

				NExpansionSeries coef(-1, ExpBetaSeries(1));
				TraceState state;
				state.AddTrace(SingleTrace::Merge(a, d, b));
				state.AddTrace(c);
				AddState(state, coef, res);

				NExpansionSeries coef2(-2, ExpBetaSeries(-1));
				TraceState state2;
				state2.AddTrace(SingleTrace::Merge(a, c, d, b));
				AddState(state2, coef2, res);

				NExpansionSeries coef3(-2, ExpBetaSeries(-1));
				TraceState state3;
				state3.AddTrace(SingleTrace::Merge(a, d, b, c));
				AddState(state3, coef3, res);

				NExpansionSeries coef4(-3, ExpBetaSeries(1));
				TraceState state4;
				state4.AddTrace(SingleTrace::Merge(c, d, b));
				state4.AddTrace(a);
				AddState(state4, coef4, res);
			}

			if (single.Bit(i) == a2 && single.Bit(j) == a1)
			{
				// i == a2, j == a1
				// 1/N*Tr(Aa1a2) Tr(BjCiD) = 1/N Tr(Aa1a2) Tr(iDBjC)
				// = 1/N*Tr(AC)Tr(DB)-1/N^2*(Tr(ACDB) + Tr(ADBC)) + 1/N^3*Tr(A)Tr(CDB)

				NExpansionSeries coef(-1, ExpBetaSeries(1));
				single.Split(j, i, b, c, d);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(a, c));
				state.AddTrace(SingleTrace::Merge(d, b));
				AddState(state, coef, res);

				NExpansionSeries coef2(-2, ExpBetaSeries(-1));
				TraceState state2;
				state2.AddTrace(SingleTrace::Merge(a, c, d, b));
				AddState(state2, coef2, res);

				NExpansionSeries coef3(-2, ExpBetaSeries(-1));
				TraceState state3;
				state3.AddTrace(SingleTrace::Merge(a, d, b, c));
				AddState(state3, coef3, res);

				NExpansionSeries coef4(-3, ExpBetaSeries(1));
				TraceState state4;
				state4.AddTrace(SingleTrace::Merge(c, d, b));
				state4.AddTrace(a);
				AddState(state4, coef4, res);
			}
		}
	}
}

void FourFieldOperator::ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res) {
	int a1 = _annihilator / 2;
	int a2 = (_annihilator & 1);
	// position of a1 and a2: Tr(...a1a2)

	SingleTrace a(_creator, 2);
	SingleTrace b, c;

	for (int i = 0; i < single1.BitNumber(); i++)
	{
		if (single1.Bit(i) != a1 && single1.Bit(i) != a2)
		{
			continue;
		}

		SingleTrace st1, st2;
		single1.Split(i, st1, st2);
		b = SingleTrace::Merge(st2, st1);
		for (int j = 0; j < single2.BitNumber(); j++)
		{
			SingleTrace st3, st4;
			single2.Split(j, st3, st4);
			c = SingleTrace::Merge(st4, st3);
			if (single1.Bit(i) == a2 && single2.Bit(j) == a1)
			{
				// i == a2, j == a1
				// 1/N*Tr(Aa1a2)Tr(iB)Tr(jC) = 
				// 1/N*Tr(ACB) - 1/N^2*(Tr(AB)Tr(C) + Tr(AC)Tr(B)) + 1/N^3*Tr(A)Tr(B)Tr(C)
				NExpansionSeries coef(-1, ExpBetaSeries(1));
				TraceState state;
				state.AddTrace(SingleTrace::Merge(a, c, b));
				AddState(state, coef, res);

				NExpansionSeries coef2(-2, ExpBetaSeries(-1));
				TraceState state2;
				state2.AddTrace(SingleTrace::Merge(a, b));
				state2.AddTrace(c);
				AddState(state2, coef2, res);

				NExpansionSeries coef3(-2, ExpBetaSeries(-1));
				TraceState state3;
				state3.AddTrace(SingleTrace::Merge(a, c));
				state3.AddTrace(b);
				AddState(state3, coef3, res);

				NExpansionSeries coef4(-3, ExpBetaSeries(1));
				TraceState state4;
				state4.AddTrace(a);
				state4.AddTrace(b);
				state4.AddTrace(c);
				AddState(state4, coef4, res);
			}

			if (single1.Bit(i) == a1 && single2.Bit(j) == a2)
			{
				// i == a1, j == a2
				// 1/N*Tr(Aa1a2)Tr(iB)Tr(jC) = 
				// 1/N*Tr(ABC) - 1/N^2*(Tr(AB)Tr(C) + Tr(AC)Tr(B)) + 1/N^3*Tr(A)Tr(B)Tr(C)

				NExpansionSeries coef(-1, ExpBetaSeries(1));
				TraceState state;
				state.AddTrace(SingleTrace::Merge(a, b, c));
				AddState(state, coef, res);

				NExpansionSeries coef2(-2, ExpBetaSeries(-1));
				TraceState state2;
				state2.AddTrace(SingleTrace::Merge(a, b));
				state2.AddTrace(c);
				AddState(state2, coef2, res);

				NExpansionSeries coef3(-2, ExpBetaSeries(-1));
				TraceState state3;
				state3.AddTrace(SingleTrace::Merge(a, c));
				state3.AddTrace(b);
				AddState(state3, coef3, res);

				NExpansionSeries coef4(-3, ExpBetaSeries(1));
				TraceState state4;
				state4.AddTrace(a);
				state4.AddTrace(b);
				state4.AddTrace(c);
				AddState(state4, coef4, res);
			}
		}
	}
}

ISimpleOperator* ISimpleOperator::Create(std::string opStr) {
	if (opStr == "NTr") {
		return new TraceNumberOperator();
	}
	else if (opStr == "NTrX") {
		return new MagnonTraceNumberOperator();
	}
	else if (opStr.length() == 2) {
		return new TwoFieldOperator(opStr);
	}
	else if (opStr.length() == 4) {
		return new FourFieldOperator(opStr);
	}
	else {
		throw OtocException("Unknown operator " + opStr);
	}
}

IOperator* IOperator::Create(std::string op) {
	std::vector<int> pos;
	int parenthsis = 0;
	int opType = -1;
	int add = 1;
	int multiply = 0;
	for (int i = 0; i < op.length(); i++) {
		if (op[i] == '(') {
			parenthsis++;
			continue;
		}
		else if (op[i] == ')') {
			parenthsis--;
			continue;
		}
		else if (parenthsis != 0) {
			continue;
		}

		if (op[i] == '+' || op[i] == '-') {
			if (opType < add) {
				opType = add;
				pos.clear();
			}
			pos.push_back(i);
		}
		else if (op[i] == '*' && opType <= multiply) {
			opType = multiply;
			pos.push_back(i);
		}
	}

	if (opType < 0) {
		if (op.find("exp(") == 0 && op.rfind(')') == op.length() - 1) {
			return new ExponentialOperator(Create(op.substr(4, op.length() - 5)));
		}
		else if (op.find("~") == 0) {
			return new NormalOrderedOperator(Create(op.substr(1, op.length() - 1)));
		}
		else if (op.length() > 2 && op[0] == '(' && op[op.length() - 1] == ')') {
			if (op.find(',') != std::string::npos && op.find('X') == std::string::npos && op.find('Z') == std::string::npos) {
				// this is const operator.
				var_t val;
				try {
					val = ToComplex(op);
				}
				catch (std::invalid_argument& iae) {
					throw OtocException(iae.what());
				}

				return new ConstantOperator(val);
			}
			else {
				return Create(op.substr(1, op.length() - 2));
			}
		}
		else {
			f_type val;
			if (ToValue(op, val)) {
				return new ConstantOperator(val);
			}
			else {
				return ISimpleOperator::Create(op);
			}
		}
	}
	else if (opType == multiply) {
		if (pos.size() > 1) {
			throw OtocException("Invalid operator expression " + op + " Cannot contain more than 1 multiplication at one level.");
		}
		std::string part1 = op.substr(0, pos[0]);
		var_t factor;
		try {
			factor = ToComplex(part1);
		}
		catch (std::invalid_argument& iae) {
			throw OtocException(iae.what());
		}
		return new ScalarMultiplyOperator(Create(op.substr(pos[0] + 1)), factor);
	}
	else if (opType == add) {
		SumOperator* ret = new SumOperator();
		if (pos[0] > 0) {
			ret->AddOperator(Create(op.substr(0, pos[0])));
		}

		pos.push_back(op.length());
		for (size_t i = 0; i < pos.size() - 1; i++) {
			std::string part = op.substr(pos[i] + 1, pos[i + 1] - pos[i] - 1);
			IOperator* inst = Create(part);
			if (op[pos[i]] == '+') {
				ret->AddOperator(inst);
			}
			else if (op[pos[i]] == '-') {
				ScalarMultiplyOperator* sop = dynamic_cast<ScalarMultiplyOperator*>(inst);
				ConstantOperator* cop = dynamic_cast<ConstantOperator*>(inst);
				if (sop != NULL) {
					sop->MultiplyBy(-1.0);
					ret->AddOperator(sop);
				}
				else if (cop != NULL) {
					cop->MultiplyBy(-1.0);
					ret->AddOperator(cop);
				}
				else {
					ret->AddOperator(new ScalarMultiplyOperator(inst, -1.0));
				}
			}
		}

		return ret;
	}
}

FourFieldOperatorB::FourFieldOperatorB(int creator, int annihilator)
{
	this->_creator = creator;
	this->_annihilator = annihilator;
}

FourFieldOperatorB::FourFieldOperatorB(std::string op)
{
	if (op.length() != 4 ||
		std::string("XZ").find(op[0]) == std::string::npos ||
		std::string("XZ").find(op[1]) == std::string::npos ||
		std::string("XZ").find(op[2]) == std::string::npos ||
		std::string("XZ").find(op[3]) == std::string::npos) {
		throw OtocException("Argument Error: op=" + op + ": must be 4-character string consisting of either X or Z.");
	}

	std::swap(op[1], op[2]);

	int n = OpStringToInteger(op);
	this->_creator = (n >> 2) & 3;
	this->_annihilator = n & 3;
}

void FourFieldOperatorB::ApplyOnSingle(const SingleTrace& single, MixState& res)
{
	int n = single.BitNumber();
	SingleTrace c, d, e;
	SingleTrace a(_creator / 2, 1);
	SingleTrace b(_creator & 1, 1);
	int a2 = (this->_annihilator & 1);
	int a1 = this->_annihilator / 2;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (single.Bit(i) == a1 && single.Bit(j) == a2)
			{
				NExpansionSeries coef(0, ExpBetaSeries(1));
				single.Split(j, i, c, d, e);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(a, e, c));
				state.AddTrace(SingleTrace::Merge(b, d));
				AddState(state, coef, res);
			}

			if (single.Bit(i) == a2 && single.Bit(j) == a1)
			{
				NExpansionSeries coef(0, ExpBetaSeries(1));
				single.Split(j, i, c, d, e);
				TraceState state2;
				state2.AddTrace(SingleTrace::Merge(a, d));
				state2.AddTrace(SingleTrace::Merge(b, e, c));
				AddState(state2, coef, res);
			}
		}
	}

	if (a1 != b.Bit(0))
	{
		return;
	}

	// for normal ordered operator, we do not need it.

/*	for (int i = 0; i < n; i++)
	{
		if (single.Bit(i) == a2)
		{
			NExpansionSeries coef(0, ExpBetaSeries(1));
			single.Split(i, c, d);
			TraceState state;
			state.AddTrace(SingleTrace::Merge(a, d, c));
			AddState(state, coef, res);
		}
	}*/
}

void FourFieldOperatorB::ApplyOnTwoSingle(const SingleTrace& single1, const SingleTrace& single2, MixState& res)
{
	int a1 = _annihilator / 2;
	int a2 = (_annihilator & 1);
	SingleTrace a(_creator / 2, 1);
	SingleTrace b(_creator & 1, 1);
	SingleTrace c, d, e, f, g, h;
	for (int i = 0; i < single1.BitNumber(); i++)
	{
		if (single1.Bit(i) != a1 && single1.Bit(i) != a2)
		{
			continue;
		}

		single1.Split(i, c, d);
		for (int j = 0; j < single2.BitNumber(); j++)
		{
			if (single1.Bit(i) == a2 && single2.Bit(j) == a1)
			{
				//..
				NExpansionSeries coef(0, ExpBetaSeries(1));
				single2.Split(j, e, f);
				g = SingleTrace::Merge(a, f, e);
				h = SingleTrace::Merge(b, d, c);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(g, h));

				AddState(state, coef, res);
			}

			if (single1.Bit(i) == a1 && single2.Bit(j) == a2)
			{
				//..
				NExpansionSeries coef(0, ExpBetaSeries(1));
				single2.Split(j, e, f);
				g = SingleTrace::Merge(a, d, c);
				h = SingleTrace::Merge(b, f, e);
				TraceState state;
				state.AddTrace(SingleTrace::Merge(g, h));
				AddState(state, coef, res);
			}
		}
	}
}

string FourFieldOperatorB::ToString() const
{
	string a = ToLower(Bits2String(this->_creator, 2));
	string b = Bits2String(this->_annihilator, 2);
	string ret = "Tr(";
	ret += a[0];
	ret += b[0];
	ret += a[1];
	ret += b[1] + ")";

	return ret;
}