#ifndef VECTOR_HPP_INCLUDED
#define VECTOR_HPP_INCLUDED 1

#include <vector>
#include <math.h>

typedef long double scalar;
typedef std::vector<scalar> vec;

scalar norm(const vec& v);

void round(vec& v);

void roundZero(vec& v);

scalar eucInnerProd(const vec& v1, const vec& v2);

bool orthogonal(const vec& v1, const vec& v2);

vec operator*(const scalar s, const vec& v);

vec operator+(const vec& lhs, const vec& rhs);

vec operator-(const vec& lhs, const vec& rhs);

void operator+=(vec& lhs, const vec& rhs);

void operator-=(vec& lhs, const vec& rhs);

bool isZero(const vec& v);

vec proj(const vec& base, const vec& v);

void normalize(vec& vec);

#endif
