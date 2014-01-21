#include "Matrix.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>

#define ROUND_CUTOFF 0.0001
//MAKE USE BOOST
//Constructors: 
Matrix::Matrix() 
  :numRows_(1), numCols_(1), boostMatrix_(zero_matrix(1,1))
{
}

Matrix::Matrix(const Matrix& copy)
  :numRows_(copy.rows()), numCols_(copy.cols()), boostMatrix_(copy.boostMatrix_)
{
}

Matrix::Matrix(const size_t rows, const size_t cols) 
  :numRows_(rows), numCols_(cols), boostMatrix_(zero_matrix(rows,cols))
{
}

Matrix::Matrix(const size_t dim)
  :numRows_(dim), numCols_(dim), boostMatrix_(identity_matrix(dim))
{
}

//Assignment Operator:
Matrix& Matrix::operator=(const Matrix& copy) {
  this->numCols_ = copy.numCols_;
  this->numRows_ = copy.numRows_;
  delete[] this->data_;
  size_t size = this->numCols_ * this->numRows_;
  this->data_ = new scalar[size];
  for (size_t i = 0; i < size; i++)
    this->data_[i] = copy.data_[i];
  return *this;
}

//Math: 
// Linear Algebra: 
//  Transpose:
Matrix Matrix::transpose() const{
  Matrix mT(this->numCols_, this->numRows_);
  for (size_t n = 0; n < this->numRows_; ++n) {
    mT.setCol(n, this->getRow(n));
  }
  return mT;
}

void Matrix::transposeIP() {
  Matrix mT = this->transpose();
  (*this) = mT;
}

//  Trace:
scalar Matrix::trace() const {
  size_t bound = std::min(this->numRows_, this->numCols_);
  scalar sum = 0;
  for (size_t i = 0; i < bound; ++i)
     sum += (*this)(i,i);
  return sum;
}

//  Row Echelon Form:
Matrix Matrix::rowEchelonForm() const {
  Matrix m = (*this);
  m.rowEchelonFormRIP(0,0); 
  return m;
}

void Matrix::rowEchelonFormIP() {
  this->rowEchelonFormRIP(0,0);
}

//  Gram Schmidt:
void Matrix::gramSchmidtColsIP() {
  for (size_t j = 0; j < this->numCols_; ++j) {
    vec jCol = this->getCol(j);
    vec copy = jCol;
    for (size_t jj = 0; jj < j; ++jj)
      jCol -= ::proj(this->getCol(jj),copy);
    this->setCol(j,jCol);
  }
  for (size_t j = 0; j < this->numCols_; ++j) {
    vec column = this->getCol(j);
    if (!::isZero(column))
      ::normalize(column);
    this->setCol(j, column);
  }
  this->roundZero();
}

Matrix Matrix::gramSchmidtCols() const {
  Matrix temp = (*this);
  temp.gramSchmidtColsIP();
  return temp;
}

//  Inversion:
Matrix Matrix::inverse() const {
  if (this->orthogonal(false)) {
    return this->orthMatInv();
  } else {
    std::cout << "I don't know how to invert this type of matrix.";
    std::cout << " It's not orthogonal." << std::endl;
    this->prettyPrint();
    return Matrix(this->numRows_, this->numCols_);
  }
}

Matrix Matrix::orthMatInv() const {
  Matrix mT = this->transpose();
  Matrix temp = mT * (*this);
  Matrix scalars(temp.rows(), temp.cols());
  for (size_t i = 0; i < temp.rows(); ++i)
    scalars(i,i) = 1.0/temp(i,i);
  return scalars * mT;
}

//  Orthogonality: 
bool Matrix::orthogonal(bool debug = false) const {
  for (size_t i = 0; i < this->numCols_; ++i) {
    for (size_t j = i+1; j < this->numCols_; ++j) {
      if (!::orthogonal(this->getCol(i), this->getCol(j))) {
        if (debug) {
          std::cout << "MATRIX NOT ORTHOGONAL" << std::endl;
          this->prettyPrint();
        }
        return false;
      }
    }
  }
  return true;
}

// Arithmetic: 
//  Multiplication:
Matrix  Matrix::operator* (const Matrix& rhs) const {
  if (this->numCols_ != rhs.numRows_) {
    std::cout << "You can't multiply a " << this->numRows_ << "x";
    std::cout << this->numCols_ << " matrix by a " << rhs.numRows_ << "x";
    std::cout << rhs.numCols_ << std::endl;
    //std::cout << "Here's left hand operand: " << std::endl;
    //this->prettyPrint();
    //std::cout << "Here's right hand operand: " << std::endl;
    //rhs.prettyPrint();
    throw;
  }
  Matrix result = Matrix(this->numRows_, rhs.numCols_);
  for (size_t row = 0; row < result.numRows_; row++) {
    for (size_t col = 0; col < result.numRows_; col++) {
      scalar sum = 0;
      for (size_t step = 0; step < this->numCols_; step++)
        sum += (*this)(row,step) * rhs(step,col);
      result(row, col) = sum;
    }
  }
  result.roundZero();
  return result;
}

Matrix& Matrix::operator*=(const Matrix& rhs) {
  Matrix result = *this * rhs;
  *this = result;
  return *this;
}

// Addition:
Matrix  Matrix::operator+ (const Matrix& rhs) const {
  Matrix copy = *this;
  copy += rhs;
  return copy;
}

Matrix  Matrix::operator- (const Matrix& rhs) const {
  return *this + rhs * (-1.0);
}

Matrix& Matrix::operator+=(const Matrix& rhs) {
  //TODO: Error checking. 
  size_t thisSize, rhsSize;
  thisSize = this->numRows_ * this->numCols_;
  rhsSize = rhs.numRows_ * rhs.numCols_;
  if ((this->numRows_ != rhs.numRows_) || (this->numCols_ != rhs.numCols_)) {
    std::cout << "You can't add matricies that aren't of the same size! ";
    std::cout << "Lhs size: " << this->numRows_ << " x " << this->numCols_;
    std::cout << " | Rhs size: " << rhs.numRows_ << " x " << rhs.numCols_;
    std::cout << "." << std::endl;
    throw;
  }
  for (size_t i = 0; i < thisSize; i++)
    this->data_[i] += rhs.data_[i];
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
  return *this += rhs * (-1.0);
}

// Scalar Multiplication: 
Matrix  Matrix::operator* (const scalar rhs) const {
  Matrix copy = *this;
  copy *= rhs;
  return copy;
}

Matrix& Matrix::operator*=(const scalar rhs) {
  size_t size = this->numRows_ * this->numCols_;
  for (size_t i = 0; i < size; i++)
    this->data_[i] *= rhs;
  return *this;
}

//  Vector Multiplication:
vec Matrix::operator*(const vec& v) const {
  vec sum;
  sum = v[0] * this->getCol(0);
  for (size_t i = 1; i < v.size(); ++i)
    sum += v[i]*this->getCol(i);
  return sum;
}

//Num-nonzero Entires:
size_t Matrix::numNonzeroEntries() const {
  size_t count = 0;
  for (size_t row = 0; row < this->numRows_; ++row)
    for (size_t col = 0; col < this->numCols_; ++col) 
      if ((*this)(row,col) != 0)
        ++count;
  return count;
}

//Element Access: 
const scalar Matrix::operator()(const size_t row, const size_t col) const {
  return this->data_[row*this->numCols_ + col];
}

scalar& Matrix::operator()(const size_t row, const size_t col) {
  return this->data_[row*this->numCols_ + col];
}

vec Matrix::getCol(const size_t colNum) const {
  vec column;
  for (size_t row = 0; row < this->numRows_; ++row)
    column.push_back((*this)(row,colNum));
  return column;
}

vec Matrix::getRow(const size_t rowNum) const {
  vec row;
  for (size_t col = 0; col < this->numCols_; ++col)
    row.push_back((*this)(rowNum, col));
  return row;
}

void Matrix::setCol(const size_t colNum, const vec& newCol) {
  for (size_t row = 0; row < this->numRows_; ++row)
    (*this)(row,colNum) = newCol[row];
}

void Matrix::setRow(const size_t rowNum, const vec& newRow) {
  for (size_t col = 0; col < this->numCols_; ++col)
    (*this)(rowNum,col) = newRow[col];
}

//Getters:
size_t Matrix::rows() const {
  return this->numRows_;
}

size_t Matrix::cols() const {
  return this->numCols_;
}

void Matrix::prettyPrintTo(std::ostream& out) const {
  std::stringstream ss;
  std::vector<size_t> colWidth;
  size_t numRows = this->numRows_;
  size_t numCols = this->numCols_;
  for (size_t j=0; j < numCols; ++j) {
    ss.str("");
    ss << (*this)(0,j);
    colWidth.push_back(ss.str().length());
  }
  for (size_t i=1; i < numRows; ++i) {
    for (size_t j=0; j < numCols; ++j) {
      ss.str("");
      ss << (*this)(i,j);
      if (ss.str().length() > colWidth[j])
        colWidth[j] = ss.str().length();
    }
  }
  std::stringstream row;
  std::stringstream colEntry;
  for (size_t i=0; i < numRows; ++i) {
    colEntry << (*this)(i,0);
    size_t colLength = colEntry.str().length();
    size_t maxColLength = colWidth[0];
    for (size_t k=colLength; k < maxColLength; ++k) 
      colEntry << " ";
    row << colEntry.str();
    colEntry.str("");
    for (size_t j=1; j < numCols; ++j) {
      row << " ";
      colEntry <<  (*this)(i,j);
      colLength = colEntry.str().length();
      maxColLength = colWidth[j];
      for (size_t k=colLength; k < maxColLength; ++k) 
        colEntry << " ";
      row << colEntry.str();
      colEntry.str("");
    }
    if (i == 0) {
      out << "⎡" << row.str() << "⎤" << std::endl;
    } else if (i == numRows - 1) {
      out << "⎣" << row.str() << "⎦" << std::endl;
    } else {
      out << "⎢" << row.str() << "⎥" << std::endl;
    }
    row.str("");
  }
}

void Matrix::prettyPrint() const {
  this->prettyPrintTo(std::cout);
}

//Rounding:

void Matrix::roundZero() {
  for (size_t i = 0; i < this->rows(); ++i) 
    for (size_t j = 0; j < this->cols(); ++j) 
      if (((*this)(i,j) < ROUND_CUTOFF) && (-(*this)(i,j) < ROUND_CUTOFF)) 
        (*this)(i,j) = 0;
}

//Silly Wrappers: TODO: Make these not needed. 
size_t Matrix::size1() const {
  return this->rows();
}

size_t Matrix::size2() const {
  return this->cols();
}

//Private:
void Matrix::rowEchelonFormRIP(const size_t colNum, const size_t startRow) {
  if (colNum == this->numCols_) 
    return ;

  size_t row = startRow;
  while (row < this->numRows_ && (*this)(row,colNum) == 0) 
    ++row;
  if (row == this->numRows_) {
    return this->rowEchelonFormRIP(colNum + 1, startRow);
  } else {
    vec leadRow = this->getRow(row);
    vec scaledLeadRow = (1/(leadRow[colNum]))*leadRow;
    for (size_t otherRow=row+1; otherRow < this->numRows_; ++otherRow) {
      vec oldRow = this->getRow(otherRow);
      scalar s = (oldRow[colNum])/(leadRow[colNum]);
      vec newRow = oldRow - s*leadRow;
      this->setRow(otherRow,newRow);
    }
    if (row != startRow) {
      vec oldRow = this->getRow(startRow);
      this->setRow(startRow,scaledLeadRow);
      this->setRow(row,oldRow);
    }
    this->roundZero();
    return this->rowEchelonFormRIP(colNum + 1, startRow + 1);
  }
}

/*
 * IDENTITY MATRIX ============================================================
 *
 * This class provides a sparse representation of identity matricies, which is
 * made possible by the fact that we do not need to explicitly store any entries
 * in the matrix as we can know all of them simply by their index pair. 
 *

IdentityMatrix::IdentityMatrix() 
  : dim_(1) 
{}

IdentityMatrix::IdentityMatrix(const size_t n)
  : dim_(n)
{}

IdentityMatrix& operator=(const IdentityMatrix& copy) {
  this->dim_ = copy.dim_;
}

Matrix operator+(const Matrix& rhs) const;
*/



Matrix operator*(const scalar lhs, const Matrix& rhs) {
  Matrix toReturn = rhs;
  for (size_t row = 0; row < rhs.rows(); row++) 
    for (size_t col = 0; col < rhs.cols(); col++)
      toReturn(row, col) *= lhs;
  return toReturn;
}



void writeTo(std::ostream& out, const Matrix& m) {
  out << m.rows() << " " << m.cols() << "\n";
  for (size_t row = 0; row < m.rows(); ++row)
    for (size_t col = 0; col < m.cols(); ++col)
      out << m(row,col) << " ";
}


Matrix readFrom(std::istream& in) {
  std::string line;
  getline(in,line);
  std::istringstream dims(line);
  size_t rows, cols;
  dims >> rows;
  dims >> cols;
  //size_t rows = atoi(line);
  //getline(in,line);
  //size_t cols = atoi(line);
  Matrix m(rows,cols);
  for (size_t row = 0; row < rows; ++row) {
    for (size_t col = 0; col < cols; ++col) {
      in >> m(row,col);
    }
  }
  return m;
}

void writeToFile(const std::string& fileName, const Matrix& m) {
  std::ofstream file;
  file.open(fileName);
  writeTo(file, m);
}

Matrix readFromFile(const std::string& fileName) {
  std::ifstream file;
  file.open(fileName);
  return readFrom(file);
}

scalar norm(const vec& v) {
  scalar sum = 0;
  for (vec::const_iterator it=v.begin(); it != v.end(); ++it)
    sum += (*it) * (*it);
  return sqrt(sum);
}

void roundZero(vec& v) {
  for (size_t i = 0; i < v.size(); ++i) 
    if (v[i] < ROUND_CUTOFF && -v[i] < ROUND_CUTOFF) 
      v[i] = 0;
}

scalar eucInnerProd(const vec& v1, const vec& v2) {
  scalar sum = 0;
  vec::const_iterator it1 = v1.begin();
  vec::const_iterator it2 = v2.begin();
  for (; it1 != v1.end(); ++it1, ++it2)
    sum += (*it1) * (*it2);
  if (abs(sum) < ROUND_CUTOFF)
    sum = 0;
  return sum;
}

bool orthogonal(const vec& v1, const vec& v2) {
  return (eucInnerProd(v1, v2) == 0);
}

vec operator*(const scalar s, const vec& v) {
  vec temp = {};
  for (vec::const_iterator it=v.begin(); it != v.end(); ++it)
    temp.push_back(s * (*it));
  return temp;
}

vec operator+(const vec& lhs, const vec& rhs) {
  vec::const_iterator itL = lhs.begin();
  vec::const_iterator itR = rhs.begin();
  vec temp;
  for (; itL != lhs.end(); ++itL, ++itR)
    temp.push_back(*itL + *itR);
  roundZero(temp);
  return temp;
}

vec operator-(const vec& lhs, const vec& rhs) {
  return lhs + ((-1)*rhs);
}

void operator+=(vec& lhs, const vec& rhs) {
  lhs = lhs + rhs;
}

void operator-=(vec& lhs, const vec& rhs) {
  lhs = lhs - rhs;
}

bool isZero(const vec& v) {
  return (norm(v) == 0);
}

vec proj(const vec& base, const vec& v) {
  if (isZero(base))
    return base;
  else 
    return ((eucInnerProd(base, v)/pow(norm(base),2)) * base);
}

void normalize(vec& v) {
  scalar mag = norm(v);
  for (vec::iterator it=v.begin(); it != v.end(); ++it) 
    (*it) = (*it)/mag;
}
