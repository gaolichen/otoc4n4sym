#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include "BlockMatrix.h"

BlockMatrix::BlockMatrix(int L) {
	init(L);
}

BlockMatrix::BlockMatrix(const BlockMatrix& bm) {
	init(bm._L);
	for (int i = 0; i < bm.size(); i++) {
		for (int j = 0; j < bm.size(); j++) {
			_mats[i][j] = bm._mats[i][j];
		}
	}
}

BlockMatrix::BlockMatrix(int L, const Array& diag) {
	init(L);
	std::vector<int> index = buildIndex(_L);
	for (int i = 0; i < size(); i++) {
		_mats[i][i] = diag.block(index[i], 0, StateCollection::Inst()->StateNumber(_L, i), 1).matrix().asDiagonal();
	}
}

BlockMatrix::~BlockMatrix() {
/*	for (size_t i = 0; i < _mats.size(); i++) {
		for (size_t j = 0; j < _mats.size(); j++) {
			if (_mats[i][j] != NULL) {
				delete _mats[i][j];
			}
		}
	}*/
}

void BlockMatrix::init(int L) {
	this->_L = L;
	_mats.resize((size_t)L + 1, std::vector<Matrix>(L + 1));
}

std::vector<int> BlockMatrix::buildIndex(int L) {
	StateCollection* inst = StateCollection::Inst();
	std::vector<int> index(L + 1);
	index[0] = 0;
	for (int i = 1; i <= L; i++) {
		index[i] = inst->StateNumber(L, i - 1) + index[i - 1];
	}

	return index;
}

bool BlockMatrix::isLadder(int* ladderSize) const {
	int diff = 0;
	bool first = true;
	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (empty(i, j)) continue;
			if (first) {
				diff = i - j;
				first = false;
			}
			else {
				if (i - j != diff) {
					return false;
				}
			}
		}
	}

	*ladderSize = diff;
	return true;
}

void BlockMatrix::fromMixing(const Array& left, const BlockMatrix& center, const Array& right) {
	init(center._L);

	std::vector<int> index = buildIndex(_L);
	for (int i = 0; i < center.size(); i++) {
		for (int j = 0; j < center.size(); j++) {
			if (center.empty(i, j)) continue;
			_mats[i][j].resize(center._mats[i][j].rows(), center._mats[i][j].cols());
			_mats[i][j].fill(0.0);

			for (int x = 0; x < _mats[i][j].rows(); x++) {
				for (int y = 0; y < _mats[i][j].cols(); y++) {
					_mats[i][j](x, y) = left(x + index[i]) * center._mats[i][j](x, y) * right(y + index[j]);
				}
//				_mats[i][j].row(x) = center._mats[i][j].row(x) * left(x + index[i]);
			}
//			for (int y = 0; y < _mats[i][j].cols(); y++) {
//				_mats[i][j].col(y) *= right(y + index[j]);
//			}
		}
	}
}

void BlockMatrix::mixWith(const Array& left, const Array& right) {
	std::vector<int> index = buildIndex(_L);
	for (int i = 0; i < this->size(); i++) {
		for (int j = 0; j < this->size(); j++) {
			if (empty(i, j)) continue;

			for (int x = 0; x < _mats[i][j].rows(); x++) {
				_mats[i][j].row(x) *= left(x + index[i]);
			}
			for (int y = 0; y < _mats[i][j].cols(); y++) {
				_mats[i][j].col(y) *= right(y + index[j]);
			}
		}
	}
}

void BlockMatrix::mixWith(const BlockMatrix& left, const BlockMatrix& right) {
	for (int i = 0; i < this->size(); i++) {
		for (int j = 0; j < this->size(); j++) {
			if (empty(i, j)) continue;
			_mats[i][j] = left._mats[i][i] * _mats[i][j] * right._mats[j][j];
		}
	}
}

void BlockMatrix::fromMixing(const BlockMatrix& left, const Array& diag, const BlockMatrix& right) {
	if (left._L != right._L) {
		throw OtocException("Dimension of left and right matrices are different.");
	}

	init(left._L);
	std::vector<int> index = buildIndex(_L);
	for (int i = 0; i < left.size(); i++) {
		for (int k = 0; k < left.size(); k++) {
			if (left.empty(i, k)) continue;
			Matrix tmp = Matrix::Zero(left._mats[i][k].rows(), left._mats[i][k].cols());
			for (int h = 0; h < tmp.cols(); h++) {
				tmp.col(h) = left._mats[i][k].col(h) * diag(h + index[k]);
			}

			for (int j = 0; j < right.size(); j++) {
				if (right.empty(k, j)) continue;
				if (this->empty(i, j)) {
					this->fill(i, j);
				}
				this->_mats[i][j].noalias() += tmp * right._mats[k][j];
			}
		}
	}
}

void BlockMatrix::setZeroAt(const std::vector<int>& zeroPos) {
	std::vector<int> index = buildIndex(_L);
	StateCollection* inst = StateCollection::Inst();
	int x = 0;
	for (int i = 0; i < zeroPos.size(); i++) {
		if (zeroPos[i] >= index[x] + inst->StateNumber(_L, x)) x++;
//		std::cout << "i=" << i << ", x=" << x << ", state#=" << inst->StateNumber(_L, x) << ", index[x]=" << index[x] << std::endl;
		int y = 0;
		for (int j = 0; j < zeroPos.size(); j++) {
			if (zeroPos[j] >= index[y] + inst->StateNumber(_L, y)) y++;
//			std::cout << "j=" << j << ", y=" << y << ", state#=" << inst->StateNumber(_L, y) << ", index[y]=" << index[y] << std::endl;
			if (empty(x, y)) {
				continue;
			}
			_mats[x][y](zeroPos[i] - index[x], zeroPos[j] - index[y]) = 0;
		}
	}
}

void BlockMatrix::diagonalize(Array& eigenValues, BlockMatrix& eigenVectors, BlockMatrix& eigenVectorsInv) const {
	eigenValues.resize(StateCollection::Inst()->StateNumber(_L));
	eigenVectors.init(_L);
	eigenVectorsInv.init(_L);

	std::vector<int> index = buildIndex(_L);
	Eigen::ComplexEigenSolver<Matrix> solver;
	for (int i = 0; i < size(); i++) {
		int n = _mats[i][i].rows();
		solver.compute(_mats[i][i], true);
		eigenValues.block(index[i], 0, n, 1) = solver.eigenvalues();
		eigenVectors._mats[i][i] = solver.eigenvectors();
		eigenVectorsInv._mats[i][i] = eigenVectors._mats[i][i].inverse();
	}
}

Matrix BlockMatrix::matrix() const {
	std::vector<int> index = buildIndex(_L);
	int size = StateCollection::Inst()->StateNumber(_L);
	Matrix ret = Matrix::Zero(size, size);
	for (int i = 0; i < this->size(); i++) {
		for (int j = 0; j < this->size(); j++) {
			if (this->empty(i, j)) {
				continue;
			}

			ret.block(index[i], index[j], this->_mats[i][j].rows(), this->_mats[i][j].cols()) = this->_mats[i][j];
		}
	}

	return ret;
}

void BlockMatrix::fill(int i, int j) {
	StateCollection* inst = StateCollection::Inst();
	_mats[i][j] = Matrix::Zero(inst->StateNumber(_L, i), inst->StateNumber(_L, j));
//	_mats[i][j]->fill(0.0);
}

Matrix& BlockMatrix::block(int i, int j) {
	return _mats[i][j];
}

BlockMatrix& BlockMatrix::operator += (const BlockMatrix& other) {
	if (this->_L != other._L) {
		throw OtocException("Matrix with different dimensions cannot add.");
	}

	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (other.empty(i, j)) {
				continue;
			}
			if (this->empty(i, j)) {
				_mats[i][j] = other._mats[i][j];
			}
			else {
				_mats[i][j] += other._mats[i][j];
			}
		}
	}

	return *this;
}

BlockMatrix& BlockMatrix::operator -= (const BlockMatrix& other) {
	if (this->_L != other._L) {
		throw OtocException("Matrix with different dimensions cannot add.");
	}

	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (other.empty(i, j)) {
				continue;
			}
			if (this->empty(i, j)) {
				_mats[i][j] = -other._mats[i][j];
			}
			else {
				_mats[i][j] -= other._mats[i][j];
			}
		}
	}

	return *this;
}

BlockMatrix BlockMatrix::operator - (const BlockMatrix& other) const {
	if (this->_L != other._L) {
		throw OtocException("Matrix with different dimensions cannot multiply.");
	}

	BlockMatrix res(this->size() - 1);

	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (this->empty(i, j) && other.empty(i, j)) {
				continue;
			}
			if (other.empty(i, j)) {
				res._mats[i][j] = this->_mats[i][j];
			}
			else if (this->empty(i, j)) {
				res._mats[i][j] = other._mats[i][j];
			}
			else {
				res._mats[i][j] = this->_mats[i][j] - other._mats[i][j];
			}
		}
	}

	return res;
}

BlockMatrix BlockMatrix::operator * (const BlockMatrix& other) const {
	if (this->_L != other._L) {
		throw OtocException("Matrix with different dimensions cannot multiply.");
	}

	BlockMatrix res(this->_L);

	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			for (int k = 0; k < size(); k++) {
				if (this->empty(i, k) || other.empty(k, j)) {
					continue;
				}
				if (res.empty(i, j)) {
					res.fill(i, j);
				}

				res._mats[i][j].noalias() += this->_mats[i][k] * other._mats[k][j];
			}
		}
	}

	return res;
}

void BlockMatrix::multiplyByLadder(const BlockMatrix& ladder, int ladderSize) {
	for (int i = 0; i < size(); i++) {
		if (ladderSize >= 0) {
			for (int j = 0; j < size(); j++) {
				int k = j + ladderSize;
				if (k >= size()) {
					_mats[i][j] = Matrix();
				}
				else {
					_mats[i][j] = _mats[i][k] * ladder._mats[k][j];
				}
			}
		}
		else {
			for (int j = size() - 1; j >= 0; j--) {
				int k = j + ladderSize;
				if (k < 0) {
					_mats[i][j] = Matrix();
				}
				else {
					_mats[i][j] = _mats[i][k] * ladder._mats[k][j];
				}
			}
		}
	}
}

BlockMatrix BlockMatrix::operator * (const var_t& val) const {
	BlockMatrix res(this->size() - 1);

	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (this->empty(i, j)) continue;
			res._mats[i][j] = _mats[i][j] * val;
		}
	}

	return res;
}

BlockMatrix& BlockMatrix::operator *= (const var_t& val) {
	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (this->empty(i, j)) continue;
			_mats[i][j] *= val;
		}
	}

	return *this;
}

BlockMatrix& BlockMatrix::operator *= (const Array& diag) {
	std::vector<int> index = buildIndex(_L);
	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (empty(i, j)) continue;
			for (int x = 0; x < _mats[i][j].cols(); x++) {
				_mats[i][j].col(x) *= diag(index[j] + x);
			}
		}
	}

	return *this;
}

BlockMatrix BlockMatrix::operator * (const Array& diag) {
	BlockMatrix ret(_L);
	std::vector<int> index = buildIndex(_L);
	for (int i = 0; i < size(); i++) {
		for (int j = 0; j < size(); j++) {
			if (empty(i, j)) continue;
			ret.fill(i, j);
			for (int x = 0; x < _mats[i][j].cols(); x++) {
				ret._mats[i][j].col(x) = _mats[i][j].col(x) * diag(index[j] + x);
			}
		}
	}

	return ret;
}

var_t BlockMatrix::trace() const {
	var_t ret(0.0);
	for (int i = 0; i < size(); i++) {
		if (!empty(i, i)) {
			ret += _mats[i][i].trace();
		}
	}

	return ret;
}

BlockMatrix BlockMatrix::exp() const {
	int ladderSize;
	if (!isLadder(&ladderSize)) {
		throw OtocException("Only ladder-like matrix can be exponentiated.");
	}

	BlockMatrix res(this->size() - 1);
	StateCollection* inst = StateCollection::Inst();
	if (ladderSize == 0) {
		for (int i = 0; i < res.size(); i++) {
			if (this->empty(i, i)) {
				res._mats[i][i] = Matrix::Identity(inst->StateNumber(this->_L, i), inst->StateNumber(this->_L, i));
			}
			else {
				res._mats[i][i] = _mats[i][i].exp();
			}
		}
	}
	else {                                
		Array id = Array::Constant(inst->StateNumber(_L), 1);
		BlockMatrix mat(_L, id);
		res += mat;
		for (int i = 1; i * std::abs(ladderSize) <= _L; i++) {
			mat.multiplyByLadder(*this, ladderSize);
			mat *= 1 / (f_type)i;
			res += mat;
		}
	}

	return res;
}

BlockMatrix BlockMatrix::conjugate() const {
	BlockMatrix res(this->size() - 1);
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res.size(); j++) {
			if (this->empty(i, j)) continue;
			res._mats[i][j] = _mats[i][j].conjugate();
		}
	}

	return res;
}

BlockMatrix BlockMatrix::adjoint() const {
	BlockMatrix res(this->size() - 1);
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res.size(); j++) {
			if (this->empty(i, j)) continue;
			res._mats[j][i] = _mats[i][j].adjoint();
		}
	}

	return res;
}

BlockMatrix& BlockMatrix::adjointInPlace() {
	for (int i = 0; i < size(); i++) {
		for (int j = i; j < size(); j++) {
			if (this->empty(i, j) && this->empty(j, i)) continue;
			if (i == j) {
				_mats[i][j].adjointInPlace();
			} else if (this->empty(i, j)) {
				_mats[i][j] = _mats[j][i].adjoint();
				_mats[j][i] = Matrix();
			}
			else if (this->empty(j, i)){
				_mats[j][i] = _mats[i][j].adjoint();
				_mats[i][j] = Matrix();
			}
			else {
				Matrix tmp = _mats[i][j];
				_mats[i][j] = _mats[j][i].adjoint();
				_mats[j][i] = tmp.adjoint();
/*				for (int x = 0; x < _mats[i][j].rows(); x++) {
					for (int y = 0; y < _mats[i][j].cols(); y++) {
						var_t tmp = _mats[i][j](x, y);
						_mats[i][j](x, y) = std::conj(_mats[j][i](y, x));
						_mats[j][i](y, x) = std::conj(tmp);
					}
				}*/
			}
		}
	}

	return *this;
}

BlockMatrix BlockMatrix::inverse() const {
	BlockMatrix res(this->size() - 1);
	for (int i = 0; i < res.size(); i++) {
		res._mats[i][i] = _mats[i][i].inverse();
	}

	return res;
}

BlockMatrix& BlockMatrix::inverseInPlace() {
	for (int i = 0; i < size(); i++) {
		_mats[i][i] = _mats[i][i].inverse();
	}

	return *this;
}

Array BlockMatrix::diagonal() const {
	std::vector<int> index = buildIndex(_L);
	StateCollection* inst = StateCollection::Inst();
	Array ret = Array::Zero(inst->StateNumber(this->_L));
	for (int m = 0; m <= this->_L; m++) {
		if (empty(m, m)) continue;
		ret.block(index[m], 0, _mats[m][m].rows(), 1) = _mats[m][m].diagonal().array();
	}

	return ret;
}

ostream& operator << (ostream& os, const BlockMatrix& bm) {
	for (int i = 0; i < bm.size(); i++) {
		for (int j = 0; j < bm.size(); j++) {
			if (bm.empty(i, j)) continue;
			os << "block(" << i << ", " << j << ")=" << std::endl;
			os << bm._mats[i][j] << std::endl;
		}
	}

	return os;
}