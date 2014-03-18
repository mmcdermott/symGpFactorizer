#include <math.h>
#include <vector>
#include <algorithm>
#include "Vector.hpp"
#define ROUND_CUTOFF 0.001

typedef float scalar;
typedef std::vector<scalar> vec;

scalar norm(const vec& v) {
  scalar sum = 0;
  for (vec::const_iterator it=v.begin(); it != v.end(); ++it)
    sum += (*it) * (*it);
  return sqrt(sum);
}

void roundZero(vec& v) {
  for (size_t i = 0; i < v.size(); ++i) 
    if (fabs(v[i]) < ROUND_CUTOFF) 
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
