#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <array>
#include <algorithm>

#include "symGp.hpp"

using namespace std;

template <class T> 
ostream& operator<<(ostream& out, const list<T>& l) {
  out << "{ | ";
  for (T val : l) {
    out << val << " | ";
  }
  out << "}";
  return out;
}

int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

vector<list<SymGpElm>> adjTransDecRec(int n) {
  if (n == 1) {
    vector<int> identity = {1};
    SymGpElm e = SymGpElm(1,identity);
    list<SymGpElm> l = {e};
    return {l};
  }
  vector<list<SymGpElm>> Sn1 = adjTransDecRec(n-1);
  vector<list<SymGpElm>> Sn = Sn1;
  //stringstream ss;
  //ss << "(" << n-1 << n << ")";
  vector<int> adjTransMap;
  for (int i = 1; i <= n-2; i++) 
    adjTransMap.push_back(i);
  adjTransMap.push_back(n);
  adjTransMap.push_back(n-1);
  SymGpElm adjTransN = SymGpElm(n, adjTransMap);//ss.str());

  list<SymGpElm> multiplier = {adjTransN};
  list<SymGpElm> Sn1ElmDecomp;
  Sn.push_back(multiplier);
  vector<list<SymGpElm>>::iterator it = Sn1.begin();
  ++it;
  for (; it != Sn1.end(); ++it) {
    Sn1ElmDecomp = *it;
    Sn1ElmDecomp.insert(Sn1ElmDecomp.end(),multiplier.begin(),multiplier.end());
    Sn.push_back(Sn1ElmDecomp);
  }
  //stringstream temp;
  for (int i = 2; i < n; ++i) {
    //ss.str("");
    //ss << "(" << n-i+1 << n-i << ")";
    
    adjTransMap[n-i+1] = n-i+2;
    adjTransMap[n-i]   = n-i;
    adjTransMap[n-i-1] = n-i+1;

    SymGpElm adjTransI = SymGpElm(n, adjTransMap);//ss.str());
    multiplier.push_front(adjTransI);
    multiplier.push_back(adjTransI);
    Sn.push_back(multiplier);
    it = Sn1.begin();
    it++;
    for (; it != Sn1.end(); ++it) {
      Sn1ElmDecomp = *it;
      Sn1ElmDecomp.insert(Sn1ElmDecomp.end(),multiplier.begin(),multiplier.end());
      Sn.push_back(Sn1ElmDecomp);
    }
  }
  return Sn;
}

vector<string> adjTransDec(int n) {
  if (n == 1) {
    return {"(1)"};
  }
  vector<string> Sn1 = adjTransDec(n-1);
  vector<string> Sn = Sn1;
  stringstream ss;
  stringstream ssSide;
  ss << "(" << n-1 << n << ")";
  stringstream result;
  vector<string>::iterator it = Sn1.begin();
  Sn.push_back(ss.str());
  it++;
  for (; it != Sn1.end(); ++it) {
    result.str("");
    result << *it << (ss.str());
    Sn.push_back(result.str());
  }
  stringstream temp;
  for (int i = 2; i < n; ++i) {
    temp.str("");
    temp << "(" << n-i+1 << n-i << ")" << ss.str() << "(" << n-i << n-i+1 << ")"; 
    ss.str(temp.str());
    it = Sn1.begin();
    Sn.push_back(ss.str());
    it++;
    for (; it != Sn1.end(); ++it) {
      result.str("");
      result << *it << ss.str();
      Sn.push_back(result.str());
    }
  }
  return Sn;
}

int hookProduct(const vector<int>& lambda) {
  int hookP = 1;
  int numRows = lambda.size();
  for (int row=1; row <= numRows; ++row) {
    int rowL = lambda[row-1];
    for (int i = 1; i <= rowL; ++i) {
      int hookL = 1 + (rowL - i);
      int height = 0;
      for (int nextRow=row; nextRow < numRows; ++nextRow) {
        if (lambda[nextRow] >= i) {
          ++height;
        } else {
          break;
        }
      }
      hookL += height;
      hookP *= hookL;
    }
  }
  return hookP;
}

int fLambda(const vector<int>& lambda, int n) {
  //Calculates dimension of irreducible representation of S_n associated with lambda. 
  return factorial(n)/hookProduct(lambda);
}

vector<vector<int>> nPartitions(int n) {
  return nPartitionsV(n, {});
}

vector<vector<int>> nPartitionsV(int n, vector<int> curr)
{
  if (n == 0) {
    return {curr};
  } else if (n == 1) {
    curr.push_back(1);
    return {curr};
  } else {
    int bound;
    if (curr.size() == 0) {
      bound = n;
    } else {
      bound = min(n, curr.back());
    }
    vector<vector<int>> result;
    vector<vector<int>> temp;
    vector<int>         option;
    for (int i = 1; i <= bound; ++i) {
      option = curr;
      option.push_back(i);
      temp = nPartitionsV(n-i, option);
      result.insert(result.end(), temp.begin(), temp.end());
    }
    return result;
  }
}

int numNPartitions(int n) {
  return numNPartitionsL(n, n);
}

int numNPartitionsL(int n, int limit) {
  if (n <= 1) {
    return 1;
  } else {
    int bound = min(n,limit);
    int sum = 0;
    for (int i = 1; i <= bound; ++i) {
      sum += numNPartitionsL(n-i, i);
    }
    return sum;
  }
}

SymGp::SymGp(int n) 
  : n_(n) {
  //Nothing to do here
}

SymGpElm::SymGpElm()
  : n_(1) {
  map_ = {1};
}

SymGpElm::SymGpElm(const SymGpElm& copy) 
  : n_(copy.n_), map_(copy.map_)
{
  //Nothing to do here
}

SymGpElm::SymGpElm(int n, vector<int> map) 
  : n_(n), map_(map)
{
  //Nothing to do here
}

SymGpElm::SymGpElm(int n, string cycle)
  : n_(n)
{
  cout << endl << endl <<  "AAAAAGGGHHGHGHG DONT CALL ME IM BAD" << endl << endl;
  /* Constructs a SymGpElm out of a more friendly, string cycle input */
  vector<int> curMap = vector<int>();
  vector<int> imdMap = vector<int>();
  for (int i = 1; i <= n; ++i) {
    curMap.push_back(i);
    imdMap.push_back(i);
  }

  bool inCycle         = false;
  int curCycleLen      = 0;
  char c               = '(';
  for (size_t i = 0; i < cycle.length(); ++i) {
    if (inCycle) {
      c = cycle[i+1];
      if (c == ')') {
        inCycle = false;
        char domain = cycle[i];
        char image  = cycle[i-curCycleLen];
        imdMap[domain - '0' - 1] = curMap[image - '0' - 1];
        curCycleLen = 0;
      } else {
        char domain = cycle[i];
        char image  = c;
        imdMap[domain - '0' - 1] = curMap[image - '0' - 1];
        ++curCycleLen;
      }
    } else if (cycle[i] == '(') {
      inCycle = true;
    } else {
      curMap = imdMap;
    }

  }
  map_ = curMap;
}

SymGpElm SymGpElm::operator*(const SymGpElm& rhs) const {
  /* Returns the left multiplication operation on calling object and rhs */
  vector<int> curMap = vector<int>();
  int n = max(this->n_, rhs.n_);
  for (int i = 1; i <= n; ++i) {
    curMap.push_back((*this)(rhs(i)));
  }
  SymGpElm result = SymGpElm(n, curMap);
  return result;
}

int SymGpElm::operator()(const int k) const {
  if (k > this->n_) {
    return k;
  } else {
    return this->map_[k-1];
  }
}

string SymGpElm::toString() const {
  list<array<int,2>> cycles = permRep();
  if (cycles.empty()) {
    return "(1)";
  }
  stringstream ss;
  for (list<array<int,2>>::iterator it=cycles.begin(); it != cycles.end(); ++it) {
    ss << "(" << (*it)[0] << (*it)[1] << ")";
  }
  return ss.str();
}

string SymGpElm::toStringExact() const {
  stringstream ss;
  ss << "[" << map_[0];
  for (size_t i = 1; i < map_.size(); ++i) {
    ss << "," << map_[i];
  }
  ss << "]";
  return ss.str();
}

void SymGpElm::print(ostream& out) const {
  out << toString();
}

void SymGpElm::printExact(ostream& out) const {
  out << toStringExact();
}

bool SymGpElm::even() const {
  list<array<int,2>> cycles = permRep();
  return (cycles.size() % 2 == 0);
}

bool SymGpElm::identity() const {
  vector<int> temp;
  for (int i=1; i <= this->n_; ++i) {
    temp.push_back(i);
  }
  return (this->map_ == temp);
}

bool SymGpElm::operator==(const SymGpElm& rhs) const {
  return (this->map_ == rhs.map_);
}

SymGpElm& SymGpElm::operator=(const SymGpElm& rhs) {
  swap(rhs);
  return *this;
}

void SymGpElm::swap(const SymGpElm& rhs) {
  this->n_   = rhs.n_;
  this->map_ = rhs.map_;
}

void SymGpElm::extend(int newN) {
  int n = this->n_;
  if (n < newN) {
    for (int i = n+1; i <= newN; ++i) {
      this->map_.push_back(i);
    }
  } else if (n > newN) {
    for (int i = n; i > newN; ++i) {
      this->map_.pop_back();
    }
  }
  this->n_ = newN;
}

int SymGpElm::n() const {
  return this->n_;
}

vector<int> SymGpElm::map() const {
  return this->map_;
}

list<array<int,2>> SymGpElm::permRep() const {
  if (this->n_ == 1) {
    return {};
  }
  //Setup
  list<array<int,2>> cycles;
  array<int,2> permutation; 
  int count = 0;
  std::vector<bool> seen;
  for (int i = 0; i < this->n_; ++i) {
    seen.push_back(false);
  }

  //Initial Case
  int base   = 1;
  int index  = 1;
  int result = (*this)(index);
  int seqCt  = 1;

  //Loop
  while (count < this->n_) {
    seen[index-1]  = true;
    count++;
    if (result != base) {
      permutation[0] = base;
      permutation[1] = result;
      index          = result;
      result         = (*this)(index);
      cycles.push_front(permutation);
    } else {
      for (int i = seqCt; i < this->n_; ++i) {
        seqCt++;
        if (!seen[i]) {
          index  = i+1;
          base   = i+1;
          result = (*this)(index);
          break;
        }
      }
    }
  }
  return cycles;
}

ostream& operator<<(ostream& out, const SymGpElm& elm) {
  elm.print(out);
  return out;
}

vector<list<SymGpElm>> adjTransDecList(int n) {
  vector<list<SymGpElm>> tempV = adjTransDecRec(n);
  for (list<SymGpElm> l : tempV) {
    for (SymGpElm sigma : l) {
      sigma.extend(n);
    }
  }
  return tempV;
}

SymGpElm prod(const list<SymGpElm>& decomp) {
  SymGpElm sigma = decomp.front();
  for (list<SymGpElm>::const_iterator it = ++(decomp.begin()); it != decomp.end(); ++it) {
    sigma = sigma * (*it);
  }
  return sigma;
}

list<SymGpElm> decompose(const SymGpElm& s) {
  int n = s.n();
  vector<list<SymGpElm>> Sn = adjTransDecList(n);
  for (list<SymGpElm> decomp : Sn){
    SymGpElm option = prod(decomp);
    option.extend(n);
    if (s == option) {
      return decomp;
    }
  }
}
