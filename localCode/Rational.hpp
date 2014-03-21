#ifndef RATIONAL_HPP_INCLUDED
#define RATIONAL_HPP_INCLUDED 1

#include <iostream>

typedef unsigned long long int longUS;
typedef long long int longS;
typedef long double floatT;

class Rational {
  public:
    //Constructors:
    Rational();
    Rational(const Rational& copy);
    Rational(const longS& numerator, const longUS& denominator);
    Rational(const longS& wholeNum);
    
    //Destrurctor is unneeded; this is a simple class. 

    //Assignment Operator:
    Rational& operator=(const Rational& copy);

    //Getters & Setters:
    longS numerator() const;
    longUS denominator() const;
    void setNum(const longS& newNum);
    void setDenom(const longS& newDenom);

    //Math:
    // Specific:
    floatT evaluate() const;
    void reduce() const;
    // Operations:
    Rational  operator* (const Rational& rhs) const;
    Rational& operator*=(const Rational& rhs);
    //  Addition:
    Rational  operator+ (const Rational& rhs) const;
    Rational  operator- (const Rational& rhs) const;
    Rational& operator+=(const Rational& rhs);
    Rational& operator-=(const Rational& rhs);
    ////  Integral Multiplication: 
    //Rational  operator* (const longS rhs) const;
    //Rational& operator*=(const longS rhs);
    //  Comparisons:
    bool operator==(const Rational& rhs) const;
    bool operator!=(const Rational& rhs) const;
    bool operator< (const Rational& rhs) const;
    bool operator<=(const Rational& rhs) const;
    bool operator> (const Rational& rhs) const;
    bool operator>=(const Rational& rhs) const;

    //Printing: 
    std::string toString() const;
    void print (std::ostream& out) const;
    void printEvaled (std::ostream& out) const;

  private:
    mutable longS numerator_;
    mutable longUS denominator_;
    //bool sign_;
};

std::ostream& operator<<(std::ostream& out, const Rational& q);

Rational abs(const Rational& x);
Rational pow(const Rational& x, const size_t exp);
longUS gcd(longUS a, longUS b);
#endif
