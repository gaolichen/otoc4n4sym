#pragma once
#include "BitUtility.h"
#include "StateCollection.h"

class BlockMatrix {
private:
	int _L;
	std::vector<std::vector<Matrix> > _mats;

	static std::vector<int> buildIndex(int L);
public:
	BlockMatrix() {
		// default constructor
		_L = 0;
	}

	BlockMatrix(int L);

	BlockMatrix(const BlockMatrix& bm);

	BlockMatrix(int L, const Array& diag);

	~BlockMatrix();

	void init(int L);

	void fill(int i, int j);

	Matrix& block(int i, int j);

	bool empty(int i, int j) const {
		return _mats[i][j].rows() == 0;
	}

	int size() const {
		return this->_L + 1;
	}

	bool isLadder(int* ladderSize) const;

	var_t trace() const;

	Array diagonal() const;

	BlockMatrix exp() const;

	BlockMatrix conjugate() const;

	BlockMatrix adjoint() const;

	BlockMatrix& adjointInPlace();

	BlockMatrix inverse() const;

	BlockMatrix& inverseInPlace();

	BlockMatrix& operator += (const BlockMatrix& other);

	BlockMatrix& operator -= (const BlockMatrix& other);

	BlockMatrix operator - (const BlockMatrix& other) const;

	BlockMatrix operator * (const BlockMatrix& other) const;
	void multiplyByLadder(const BlockMatrix& ladder, int ladderSize);

	BlockMatrix operator * (const var_t& val) const;
	BlockMatrix& operator *= (const var_t& val);

	BlockMatrix& operator *= (const Array& diag);
	BlockMatrix operator * (const Array& diag);

	// only diagonalize diagonal block matrix
	void diagonalize(Array& eigenValues, BlockMatrix& eigenVectors, BlockMatrix& eigenVectorsInv) const;

	// *this = DiagonalMatrix(left) * center * DiagonalMatrix(right)
	void fromMixing(const Array& left, const BlockMatrix& center, const Array& right);

	// *this = DiagonalMatrix(left) * this * DiagonalMatrix(right)
	void mixWith(const Array& left, const Array& right);

	// *this = left * this * right
	void mixWith(const BlockMatrix& left, const BlockMatrix& right);

	// both the left and right are diagonal block matrices
	void fromMixing(const BlockMatrix& left, const Array& diag, const BlockMatrix& right);

	void setZeroAt(const std::vector<int>& zeroPos);

	Matrix matrix() const;

	friend ostream& operator << (ostream& os, const BlockMatrix& coef);
};