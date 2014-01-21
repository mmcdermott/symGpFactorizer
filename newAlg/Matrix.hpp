#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED 1

#include <vector>
#include <math.h>
#include <string>
#include <iostream>
#include "Vector.hpp"

//TODO: Make this a template. 
class Matrix {
  public:
    //Constructors: 
    Matrix();
    Matrix(const Matrix& copy);
    Matrix(const size_t rows, const size_t cols);
    Matrix(const size_t dim);
    //Matrix(const IdentityMatrix& copy);
    //Destructor:
    ~Matrix();

    //Assignment Operators:
    Matrix& operator=(const Matrix& copy);

    //Math: 
    // Linear Algebra: 
    //  Transpose:
    Matrix transpose() const;
    void transposeIP();
    //  Trace:
    scalar trace() const;
    //  Row Echelon Form:
    Matrix rowEchelonForm() const;
    void rowEchelonFormIP(); // In place
    //  Gram Scmidt Column Orthagonalization:
    Matrix gramSchmidtCols() const;
    void gramSchmidtColsIP();
    //  Inverse:
    Matrix inverse() const;
    //Orthogonality: 
    bool orthogonal(bool debug) const;

    // Arithmetic:
    //  Multiplication:
    Matrix  operator* (const Matrix& rhs) const;
    Matrix& operator*=(const Matrix& rhs);
    //  Addition:
    Matrix  operator+ (const Matrix& rhs) const;
    Matrix  operator- (const Matrix& rhs) const;
    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator-=(const Matrix& rhs);
    //  Scalar Multiplication: 
    Matrix  operator* (const scalar rhs) const;
    Matrix& operator*=(const scalar rhs);
    //  Vector Operations:
    vec     operator* (const vec& rhs) const;
    //  # of nonzero entries:
    size_t numNonzeroEntries() const;

    //Element Access: 
    const scalar operator()(const size_t row, const size_t col) const;
    scalar& operator()(const size_t row, const size_t col);
    vec getCol(const size_t colNum) const;
    vec getRow(const size_t rowNum) const;
    void setCol(const size_t colNum, const vec& newCol);
    void setRow(const size_t rowNum, const vec& newRow);
    
    //Getters:
    size_t rows() const;
    size_t cols() const;

    //Printing: 
    void prettyPrintTo(std::ostream& out) const;
    void prettyPrint() const;

    //Rounding: 
    void roundZero();

    //Silly Wrappers: TODO: Make these not needed.
    size_t size1() const;
    size_t size2() const;
  private:
    size_t numRows_;
    size_t numCols_;
    scalar* data_;
    
    void rowEchelonFormRIP(const size_t colNum, const size_t startRow);
    Matrix orthMatInv() const;
};

Matrix operator*(const scalar lhs, const Matrix& rhs);

void writeTo(std::ostream& out, const Matrix& m);

Matrix readFrom(std::istream& in);

void writeToFile(const std::string& fileName, const Matrix& m);

Matrix readFromFile(const std::string& fileName);
#endif
