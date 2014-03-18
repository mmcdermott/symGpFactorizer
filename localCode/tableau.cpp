#include <iostream>
#include <sstream>
#include <vector>
#include "tableau.hpp"
#include "string"


void printVec(std::vector<LambdaTableau> Vec) {
  std::cout << "VecBegin: " << std::endl;
  std::vector<LambdaTableau>::iterator it = Vec.begin();
  std::cout << *it;
  ++it;
  for (; it != Vec.end(); ++it) {
    std::cout << "\n" << *it;
  }
  std::cout << "VecEnd." << std::endl;
}

void printVec(std::vector<double> Vec) {
  std::cout << "[";
  std::vector<double>::iterator it = Vec.begin();
  std::cout << *it;
  ++it;
  for (; it != Vec.end(); ++it) {
    std::cout << ", " << *it;
  }
  std::cout << "]" << std::endl;
}

void printVec(std::vector<int> Vec) {
  std::cout << "[";
  std::vector<int>::iterator it = Vec.begin();
  std::cout << *it;
  ++it;
  for (; it != Vec.end(); ++it) {
    std::cout << ", " << *it;
  }
  std::cout << "]" << std::endl;
}

void printVec(std::vector<std::string> Vec) {
  std::cout << "[";
  std::vector<std::string>::iterator it = Vec.begin();
  std::cout << *it;
  ++it;
  for (; it != Vec.end(); ++it) {
    std::cout << ", " << *it;
  }
  std::cout << "]" << std::endl;
}

void printVecOfVecs(std::vector<std::vector<int>> Vec) {
  for (std::vector<std::vector<int>>::iterator it = Vec.begin(); it != Vec.end(); ++it) {
    printVec(*it);
  }
}

void printVecOfVecs(std::vector<std::vector<double>> Vec) {
  for (std::vector<std::vector<double>>::iterator it = Vec.begin(); it != Vec.end(); ++it) {
    printVec(*it);
  }
}

//Xlambda
Xlambda::Xlambda() 
  :lambda_({1})
{}

//CXlambda
CXlambda::CXlambda()
  :lambda_({1}), stdBasis_({LambdaTableau({1})})
{}

CXlambda::CXlambda(std::vector<int> lambda)
  :lambda_(lambda)
{
  stdBasis_ = CXlambdaBasis(lambda);
}

//LambdaTableau
LambdaTableau::LambdaTableau()
  :pos_({1, 1, 1, 2, 2})
{}

LambdaTableau::LambdaTableau(const LambdaTableau& copy)
  :pos_(copy.pos())
{}

//LambdaTableau::LambdaTableau(std::vector<int> lambda, bool lambdaUsed) {
//  std::vector<int> pos;
//  int count = 1;
//  for (std::vector<int>::iterator it = lambda.begin(); it != lambda.end(); ++it) {
//    for (int i = 0; i < *it; ++i) {
//      pos.push_back(count);
//    }
//    ++count;
//  }
//  this->pos_ = pos;
//}

LambdaTableau::LambdaTableau(std::vector<int> pos)
  :pos_(pos)
{}

std::string LambdaTableau::toString() const
{
  std::stringstream ss;
  size_t count = 0;
  int miniC = 1;
  int level = 1;
  while (count < pos_.size()) {
    ss << "|";
    for (std::vector<int>::const_iterator it = pos_.begin(); it != pos_.end(); it++) {
      if (*it == level) {
        ss << miniC << "|";
        count++;
      } 
      miniC++;
    }
    miniC = 1;
    level++;
    ss << "\n";
  }
  return ss.str();
}

std::string LambdaTableau::toStringExact() const
{
  std::stringstream ss;
  ss << "|";
  for (std::vector<int>::const_iterator it = pos_.begin(); it != pos_.end(); it++) {
    ss << (*it) << "|";
  }
  return ss.str();
}

void LambdaTableau::print(std::ostream& out) const {
  out << this->toString() << std::endl;
}

void LambdaTableau::printExact(std::ostream& out) const {
  out << this->toStringExact() << std::endl;
}

bool LambdaTableau::sameRow(int i, int j) const {
  int n = this->pos().size();
  if (i > n || j > n) {
    std::cout << "Error: one of i=";
    std::cout << i <<" or j=" << j << " is not in this Tableau.";
    std::cout << std::endl;
    return false;
  } else {
    return ((this->pos())[i-1] == (this->pos())[j-1]);
  }
}

bool LambdaTableau::sameCol(int i, int j) const {
  int n = this->pos().size();
  if (i > n || j > n) {
    std::cout << "Error: one of i=";
    std::cout << i <<" or j=" << j << " is not in this Tableau.";
    std::cout << std::endl;
    return false;
  } else {
    int iRow = (this->pos())[i-1];
    int jRow = (this->pos())[j-1];
    int iCol = 0;
    int jCol = 0;
    for (int k = 0; k < std::max(i,j); ++k) {
      if ((this->pos())[k] == iRow && k < i-1) {
        ++iCol;
      } else if ((this->pos())[k] == jRow && k < j-1) {
        ++jCol;
      }
    }
    return (iCol == jCol);
  }
}

int LambdaTableau::axialDistance(int i, int j) const {
  int n = this->pos().size();
  if (i > n || j > n) {
    std::cout << "Error: one of i=";
    std::cout << i <<" or j=" << j << " is not in this Tableau.";
    std::cout << std::endl;
    return false;
  } else {
    int iRow = (this->pos())[i-1];
    int jRow = (this->pos())[j-1];
    int iCol = 0;
    int jCol = 0;
    for (int k = 0; k < std::max(i,j); ++k) {
      if ((this->pos())[k] == iRow && k < i-1) {
        ++iCol;
      } else if ((this->pos())[k] == jRow && k < j-1) {
        ++jCol;
      }
    }
    return (iCol - iRow) + (jRow - jCol);
  }
}

std::vector<int> LambdaTableau::pos() const
{
  return pos_;
}

bool LambdaTableau::operator<(const LambdaTableau& rhs) const {
  return ((*this) <= rhs) && ((*this) != rhs);
}

bool LambdaTableau::operator>(const LambdaTableau& rhs) const {
  return !((*this) <= rhs);
}

bool LambdaTableau::operator<=(const LambdaTableau& rhs) const {
  size_t n = this->pos_.size();
  if (n != rhs.pos_.size())
    return false;
  for (size_t row = 1; row <= n; row++) {
    size_t thisCount = 0;
    size_t rhsCount = 0;
    for (size_t i = 0; i < n; i++) {
      if (this->pos_[i] == row)
        thisCount++;
      if (rhs.pos_[i] == row) 
        rhsCount++;
    }
    if (thisCount > rhsCount) {
      return false;
    }
  }
  return true;
}

bool LambdaTableau::operator>=(const LambdaTableau& rhs) const {
  return !((*this) < rhs);
}

bool LambdaTableau::operator==(const LambdaTableau& rhs) const {
  return (this->pos_ == rhs.pos_);
}

bool LambdaTableau::operator!=(const LambdaTableau& rhs) const {
  return (this->pos_ != rhs.pos_);
}

LambdaTableau operator*(const SymGpElm& lhs, const LambdaTableau& rhs)
{
  std::vector<int> newPos = rhs.pos();
  for (size_t i = 1; i <= rhs.pos().size(); ++i) {
    newPos[i-1] = rhs.pos()[lhs(i)-1];
  }
  return LambdaTableau(newPos);
}

std::ostream& operator<<(std::ostream& out, const LambdaTableau& lt) {
  lt.print(out);
  return out;
}

std::vector<LambdaTableau> CXlambdaBasis(std::vector<int> lambda) {
  int n = 0;
  for (std::vector<int>::iterator it=lambda.begin(); it != lambda.end(); ++it) {
    n = n + *it;
  }

  int numRows = lambda.size();
  std::vector<std::vector<int>> positions;
  std::vector<std::vector<int>> lambdas;

  std::vector<int> imdPos;
  std::vector<int> imdLambda = lambda;
  for (int i = 1; i <= numRows; ++i) {
    imdPos.clear();
    imdLambda = lambda;
    if (imdLambda[i-1] != 0) {
      imdLambda[i-1] = imdLambda[i-1] - 1;
      imdPos.push_back(i);
      positions.push_back(imdPos);
      lambdas.push_back(imdLambda);
    }
  }
  std::vector<int> curPos;
  std::vector<int> curLambda;
  std::vector<std::vector<int>> imdPositions;
  std::vector<std::vector<int>> imdLambdas;
  std::vector<bool> seen;
  for (int i = 2; i <= n; ++i) {
    imdPositions.clear();
    imdLambdas.clear();
    for (size_t i = 0; i < positions.size(); ++i) {
      curPos    = positions[i];
      curLambda = lambdas[i];
      for (int i = 1; i <= numRows; ++i) {
        imdPos = curPos;
        imdLambda = curLambda;
        if (imdLambda[i-1] != 0) {
          imdLambda[i-1] = imdLambda[i-1] - 1;
          imdPos.push_back(i);
          imdPositions.push_back(imdPos);
          imdLambdas.push_back(imdLambda);
        }
      }
    }

    positions.clear();
    lambdas.clear();
    seen.clear();
    seen.assign(imdLambdas.size(),false);
    for (size_t j = 0; j < imdLambdas.size(); ++j) {
      if (!seen[j]) {
        curLambda = imdLambdas[j];
        for (size_t i = j; i < imdLambdas.size(); ++i) {
          if (imdLambdas[i] == curLambda) {
            seen[i] = true;
            lambdas.push_back(curLambda);
            positions.push_back(imdPositions[i]);
          }
        }
      }
    }
  }

  std::vector<LambdaTableau> basis;
  for (size_t i = 0; i < positions.size(); ++i) {
    basis.push_back(LambdaTableau(positions[i]));
  }
  return basis;
}


std::vector<LambdaTableau> youngTabloids(const std::vector<int> &lambda) {
  if (lambda.size() == 1) {
    std::vector<int> pos;
    for (int i = 1; i <= lambda[0]; ++i) {
      pos.push_back(1);
    }
    return {LambdaTableau(pos)};
  }
  std::vector<int> possibleEnds;
  for (size_t i = 1; i < lambda.size(); ++i) {
    if (lambda[i] < lambda[i-1]) {
      possibleEnds.push_back(i-1);
    }
  }
  possibleEnds.push_back(lambda.size()-1);
  std::vector<LambdaTableau> tabloids;
  std::vector<int> newLambda;
  for (size_t i = 0; i < possibleEnds.size(); ++i) {
    newLambda = lambda;
    int row = possibleEnds[i];
    newLambda[row] -= 1;
    if (newLambda[row] == 0) {
      newLambda.pop_back();
    }
    std::vector<LambdaTableau> result = youngTabloids(newLambda);
    for (std::vector<LambdaTableau>::iterator it=result.begin(); it!=result.end(); ++it) {
      std::vector<int> resPos = it->pos();
      resPos.push_back(possibleEnds[i]+1);
      (*it) = LambdaTableau(resPos);
    }
    tabloids.insert(tabloids.end(),result.begin(),result.end());
  }
  return tabloids;
}

