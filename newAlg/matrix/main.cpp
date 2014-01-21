#include "myMatrix.hpp"

int main(const int argc, const char* argv[]) {
  Matrix m(5,5);
  m(1,2) = 3;
  Matrix m2 = m;
  m *= 4;
  m2(0,0) = 1;
  m2(0,1) = 1;
  m2(0,2) = 1;
  m2(0,3) = 1;
  m2(0,4) = 1;
  m2(2,2) = 5;
  m(0,0) = 18;
  m(1,0) = 4;
  m(2,0) = 4;
  m(3,0) = 8;
  Matrix m3 = m * m2;
  std::cout << "m: " << std::endl;
  m.prettyPrint();
  std::cout << "m2: " << std::endl;
  m2.prettyPrint();
  std::cout << "m3: " << std::endl;
  m3.prettyPrint();
  //const tests: 
  //
  const Matrix constMatrix(5,5);
  std::cout << "constMatrix: " << std::endl;
  constMatrix.prettyPrint();
  std::cout << "Const Matrix entry 3,3: " << constMatrix(3,3) << std::endl;
  return 0;
}
