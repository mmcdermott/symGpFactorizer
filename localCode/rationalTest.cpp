#include "Rational.hpp"

using namespace std;

int main(int argc, const char* argv[]) {
  Rational one = Rational();
  //cout << "One: " << one << endl;
  Rational oneHalf = Rational(1,2);
  //cout << "One Half: " << oneHalf << endl;
  //oneHalf.print(cout);
  //cout << "One Half * 2 = " << oneHalf * 2 << endl;
  //Rational q = Rational(136984,153986944);
  //q.print(cout);
  //cout << "q = " << q << endl;
  //cout << "q - 3 = " << q - 3 << endl;
  //cout << "abs(q - 3) = " << abs(q-3) << endl;
  Rational tf = Rational(3,4);
  //cout << "tf^3 = " << pow(tf, 3) << endl;
  cout << "3/4 < 1/2: " << (tf < oneHalf) << endl;
  cout << "(3/4)*(-1) = "<< (tf * (-1)) << " or, in exact form, -3/4 = ";
  (tf * (-1)).print(cout);
  cout << "1/2-3/4 = "; 
  (oneHalf - tf).print(cout);
  cout << "N[1/2-3/4] = " << (oneHalf - tf) << endl;
  //cout << "2 > 1999/1000: " << (Rational(2) > Rational(1999,1000)) << endl;
  //cout << "3/4 == 3/4: " << (tf == tf) << endl;
  //cout << "3/4 == 1: " << (tf == 1) << endl;
  cout << "gcd(2,8) = " << gcd(2,8) << endl;
  cout << "(-2)/gcd(2,8) = " << (-2/gcd(2,8)) << endl;
  return 0;
}
