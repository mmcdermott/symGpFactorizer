#ifndef SYM_GP_HPP_INCLUDED
#define SYM_GP_HPP_INCLUDED 1

#include <iostream>
#include <vector>
#include <list>
#include <array>

int factorial(int n);

std::vector<std::string> adjTransDec(int n);

int hookProduct(const std::vector<int>& lambda);

int fLambda(const std::vector<int> &lambda, int n);

std::vector<std::vector<int>> nPartitions(int n);

std::vector<std::vector<int>> nPartitionsV(int n, std::vector<int> curr);

int numNPartitions(int n);

int numNPartitionsL(int n, int limit);

class SymGp {
  public:
    SymGp(int n);

  private:
    int n_;
};

class SymGpElm { 
  public: 
    SymGpElm();

    SymGpElm(const SymGpElm& copy);

    SymGpElm(int n, std::vector<int> map);

    SymGpElm(int n, std::string cycle);

    SymGpElm operator * (const SymGpElm& rhs) const;

    int operator()(const int k) const;

    std::string toString() const;

    std::string toStringExact() const;

    void extend(int newN);

    void print(std::ostream& out) const;
    
    void printExact(std::ostream& out) const;
    
    bool even() const;

    bool identity() const;

    bool operator==(const SymGpElm& rhs) const;

    SymGpElm& operator=(const SymGpElm& rhs);

    int n() const;

    std::vector<int> map() const;
  private: 
    int n_;
    std::vector<int> map_;

    void swap(const SymGpElm& rhs);
    std::list<std::array<int,2>> permRep() const;
};

std::ostream& operator<<(std::ostream& out, const SymGpElm& elm);

std::vector<std::list<SymGpElm>> adjTransDecList(int n);
std::vector<std::list<SymGpElm>> adjTransDecRec(int n);

std::list<SymGpElm> decompose(const SymGpElm& s);
SymGpElm prod(const std::list<SymGpElm>& decomp);
#endif
