#include "Rational.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
//GCD:
longUS gcd(longUS a, longUS b) {
  longUS c = a % b;
  while (c != 0) {
    a = b;
    b = c;
    c = a % b;
  }
  return b;
}

//Constructors:
Rational::Rational()
  :numerator_(1), denominator_(1)
{}

Rational::Rational(const Rational& copy)
  :numerator_(copy.numerator_),denominator_(copy.denominator_)
{}

Rational::Rational(const longS& numerator, const longUS&denominator)
  :numerator_(numerator),denominator_(denominator)
{
  this->reduce();
}

Rational::Rational(const longS& wholeNum)
  :numerator_(wholeNum),denominator_(1)
{}

Rational& Rational::operator=(const Rational& copy) {
  this->numerator_ = copy.numerator_;
  this->denominator_ = copy.denominator_;
  return *this;
}

//Getters & Setters
longS Rational::numerator() const {
  return this->numerator_;
}

longUS Rational::denominator() const {
  return this-> denominator_;
}

void Rational::setNum(const longS& newNum) {
  this->numerator_ = newNum;
}
void Rational::setDenom(const longS& newDenom) {
  this->denominator_ = newDenom;
}

//Math:
// Evaluation:
floatT Rational::evaluate() const {
  return ((floatT) this->numerator_)/this->denominator_;
}
void Rational::reduce() const {
  longS d = gcd(abs(this->numerator_),this->denominator_);
  if (d != 1) {
    this->numerator_ = this->numerator_ / d;
    this->denominator_ = this->denominator_ / d;
  }
}
// Operations:
Rational  Rational::operator* (const Rational& rhs) const {
  return Rational(this->numerator_*rhs.numerator_, 
                  this->denominator_*rhs.denominator_);
}

Rational& Rational::operator*=(const Rational& rhs) {
  this->numerator_ *= rhs.numerator_;
  this->denominator_ *= rhs.denominator_;
  this->reduce();
  return (*this);
}
//  Addition:
Rational  Rational::operator+ (const Rational& rhs) const {
  std::cout << "In +: " << std::endl;
  std::cout << "this = "; 
  this->print(std::cout);
  std::cout << "rhs = "; 
  rhs.print(std::cout);
  longUS newDenom = this->denominator_ * rhs.denominator_;
  longS newNumerator = this->numerator_*rhs.denominator_ + 
                       rhs.numerator_*this->denominator_;
  std::cout << "this + rhs = " << newNumerator << "/" << newDenom << std::endl;
  std::cout << "Leaving +" << std::endl;
  return Rational(newNumerator,newDenom);
}
Rational  Rational::operator- (const Rational& rhs) const {
  return (*this) + (rhs * (-1));
}
Rational& Rational::operator+=(const Rational& rhs) {
  Rational temp = (*this) + rhs;
  (*this) = temp;
  return (*this);
}
Rational& Rational::operator-=(const Rational& rhs) {
  Rational temp = (*this) - rhs;
  (*this) = temp;
  return (*this);
}
////  Integral Multiplication: 
//Rational  Rational::operator* (const longS rhs) const {
//  return Rational(this->numerator_*rhs,this->denominator_);
//}
//Rational& Rational::operator*=(const longS rhs) {
//  this->numerator_ *= rhs;
//  this->reduce();
//  return (*this);
//}
//  Comparisons
bool Rational::operator==(const Rational& rhs) const {
  return ((this->numerator_==rhs.numerator_) && 
          (this->denominator_==rhs.denominator_));
}
bool Rational::operator!=(const Rational& rhs) const {
  return !((*this) == rhs);
}
bool Rational::operator< (const Rational& rhs) const {
  Rational diff = rhs - (*this);
  return (diff.numerator_ > 0);
}
bool Rational::operator<=(const Rational& rhs) const {
  return ((*this) < rhs || (*this) == rhs);
}
bool Rational::operator> (const Rational& rhs) const {
  return !((*this) <= rhs);
}
bool Rational::operator>=(const Rational& rhs) const {
  return !((*this) < rhs);
}

//Printing:
std::string Rational::toString() const {
  std::stringstream ss;
  ss << this->numerator_ << "/" << this->denominator_;
  return (ss.str());
}
void Rational::print(std::ostream& out) const {
  out << this->toString() << std::endl;
}
void Rational::printEvaled(std::ostream& out) const {
  out << this->evaluate() << std::endl;
}


Rational abs(const Rational& x) {
  return Rational(abs(x.numerator()),x.denominator());
}

Rational fabs(const Rational& x) {
  return abs(x);
}

//Crappy helper function: eliminate me!
longS pow(const longS& base, size_t exp) {
  longS result = base;
  while (exp > 1) {
    result *= base;
    exp--;
  }
  return result;
}
longUS pow(const longUS& base, size_t exp) {
  longUS result = base;
  while (exp > 1) {
    result *= base;
    exp--;
  }
  return result;
}

Rational pow(const Rational& x, size_t exp) {
  return Rational(pow(x.numerator(),exp),pow(x.denominator(),exp));
}

std::ostream& operator<<(std::ostream& out, const Rational& q) {
  out << q.evaluate();
  return out;
}

//Integral Comparison
bool operator==(const longUS& lhs, const Rational& rhs) {
  return Rational(lhs) == rhs;
}
bool operator!=(const longUS& lhs, const Rational& rhs) {
  return Rational(lhs) != rhs;
}
bool operator< (const longUS& lhs, const Rational& rhs) {
  return Rational(lhs) < rhs;
}
bool operator<=(const longUS& lhs, const Rational& rhs) {
  return Rational(lhs) <= rhs;
}
bool operator> (const longUS& lhs, const Rational& rhs) {
  return Rational(lhs) > rhs;
}
bool operator>=(const longUS& lhs, const Rational& rhs) {
  return Rational(lhs) >= rhs;
}

bool operator==(const longS& lhs, const Rational& rhs) {
  return Rational(lhs) == rhs;
}
bool operator!=(const longS& lhs, const Rational& rhs) {
  return Rational(lhs) != rhs;
}
bool operator< (const longS& lhs, const Rational& rhs) {
  return Rational(lhs) < rhs;
}
bool operator<=(const longS& lhs, const Rational& rhs) {
  return Rational(lhs) <= rhs;
}
bool operator> (const longS& lhs, const Rational& rhs) {
  return Rational(lhs) > rhs;
}
bool operator>=(const longS& lhs, const Rational& rhs) {
  return Rational(lhs) >= rhs;
}

bool operator==(const Rational& lhs, const longUS& rhs) { 
  return lhs == rhs;
}
bool operator!=(const Rational& lhs, const longUS& rhs) { 
  return lhs != rhs;
}
bool operator< (const Rational& lhs, const longUS& rhs) { 
  return lhs > rhs;
}
bool operator<=(const Rational& lhs, const longUS& rhs) { 
  return lhs >= rhs;
}
bool operator> (const Rational& lhs, const longUS& rhs) { 
  return lhs < rhs;
}
bool operator>=(const Rational& lhs, const longUS& rhs) { 
  return lhs <= rhs;
}

bool operator==(const Rational& lhs, const longS& rhs) {
  return lhs == rhs;
}
bool operator!=(const Rational& lhs, const longS& rhs) {
  return lhs != rhs;
}
bool operator< (const Rational& lhs, const longS& rhs) {
  return lhs > rhs;
}
bool operator<=(const Rational& lhs, const longS& rhs) {
  return lhs >= rhs;
}
bool operator> (const Rational& lhs, const longS& rhs) {
  return lhs < rhs;
}
bool operator>=(const Rational& lhs, const longS& rhs) {
  return lhs <= rhs;
}

// Float Comparisons:
bool operator==(const floatT lhs, const Rational& rhs) {
  return lhs == rhs.evaluate();
}
bool operator!=(const floatT lhs, const Rational& rhs) {
  return lhs != rhs.evaluate();
}
bool operator< (const floatT lhs, const Rational& rhs) {
  return lhs < rhs.evaluate();
}
bool operator<=(const floatT lhs, const Rational& rhs) {
  return lhs <= rhs.evaluate();
}
bool operator> (const floatT lhs, const Rational& rhs) {
  return lhs > rhs.evaluate();
}
bool operator>=(const floatT lhs, const Rational& rhs) {
  return lhs >= rhs.evaluate();
}
// Reversed: 
bool operator==(const Rational& rhs, const floatT lhs) {
  return lhs == rhs;
}
bool operator!=(const Rational& rhs, const floatT lhs) {
  return lhs != rhs;
}
bool operator< (const Rational& rhs, const floatT lhs) {
  return lhs > rhs;
}
bool operator<=(const Rational& rhs, const floatT lhs) {
  return lhs >= rhs;
}
bool operator> (const Rational& rhs, const floatT lhs) {
  return lhs < rhs;
}
bool operator>=(const Rational& rhs, const floatT lhs) {
  return lhs <= rhs;
}
