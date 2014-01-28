#ifndef TABLEAU_HPP_INCLUDED
#define TABLEAU_HPP_INCLUDED 1

#include <iostream>
#include <vector>
#include "symGp.hpp"

class Xlambda {
  public:
    Xlambda();

    Xlambda(std::vector<int> lambda);
  
    std::vector<int> lambda();
  private: 
    std::vector<int> lambda_;
};


class CXlambdaVec {
  public:
    CXlambdaVec();

    CXlambdaVec(std::vector<int> lambda);

    std::vector<int> lambda();
  private:
    std::vector<int> lambda_;

    std::vector<int> stdBasisDecomp_;
};

class LambdaTableau {
  public:
    LambdaTableau();

    LambdaTableau(const LambdaTableau& copy);

    //LambdaTableau(const std::vector<int> lambda, bool lambdaUsed);

    LambdaTableau(std::vector<int> pos);
    
    std::string toString() const;

    std::string toStringExact() const;
    
    void print(std::ostream& out) const;

    void printExact(std::ostream& out) const;

    bool sameRow(int i, int j) const;
    
    bool sameCol(int i, int j) const;

    int axialDistance(int i, int j) const;

    std::vector<int> pos() const;

    bool operator<(const LambdaTableau& rhs) const;

    bool operator>(const LambdaTableau& rhs) const;

    bool operator<=(const LambdaTableau& rhs) const;

    bool operator>=(const LambdaTableau& rhs) const;

    bool operator==(const LambdaTableau& rhs) const;

    bool operator!=(const LambdaTableau& rhs) const;
  private: 
    std::vector<int> pos_;
};

std::ostream& operator<<(std::ostream& out, const LambdaTableau& lt);

LambdaTableau operator*(const SymGpElm& lhs, const LambdaTableau& rhs);

std::vector<LambdaTableau> CXlambdaBasisPermuted(std::vector<int> lambda, const SymGpElm& sigma);

std::vector<LambdaTableau> CXlambdaBasis(std::vector<int> lambda);

class CXlambda {
  public: 
    CXlambda();

    CXlambda(std::vector<int> lambda);

    std::vector<int> lambda();

    std::vector<LambdaTableau> basis();
  private: 
    std::vector<int> lambda_;

    std::vector<LambdaTableau> stdBasis_;
};

//void printVec(std::vector<LambdaTableau> Vec);
//
//void printVec(std::vector<int> Vec);
//void printVec(std::vector<double> Vec);
//
//void printVec(std::vector<std::string> Vec);

void printVecOfVecs(std::vector<std::vector<int>> Vec);
void printVecOfVecs(std::vector<std::vector<double>> Vec);

std::vector<LambdaTableau> youngTabloids(const std::vector<int> &lambda);
#endif
