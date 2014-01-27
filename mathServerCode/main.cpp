#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <utility>
#include "symGp.hpp"
#include "tableau.hpp"
#include "Matrix.hpp"
#include "printUtils.cpp"
using namespace std;

bool permuteCols(Matrix& m, const SymGpElm& sigma) {
  /*sigma.n() <= m.cols();
   */
  Matrix mCopy = m;
  vector<vec> cols;
  for (size_t i = 0; i < m.cols(); ++i) {
    cols.push_back(m.getCol(i));
  }
  for (size_t i = 0; i < sigma.n(); ++i) {
    m.setCol(sigma(i+1)-1, cols[i]);
  }
  return (mCopy != m);
}

bool canRead(ifstream& ifile) {
  if (ifile) {
    if ((ifile.peek() != ifstream::traits_type::eof())) {
      return true;
    } else {
      //cout << "File is empty!" << endl;
    }
  } else {
    //cout << "File doesn't exist!" << endl;
  }
  return false;
  //return (ifile && (ifile.peek() != ifstream::traits_type::eof()));
}

//void createDir(const path& dir) {
//  if (!boost::filesystem::exists(dir)) boost::filesystem::create_directory(dir);
//}

void createDir(const string& dir) {
  mkdir(dir.c_str(), 0755);
}

Matrix computeMatrixAdj(const SymGpElm& adjTrans, 
    const vector<LambdaTableau>& standardTableaux, 
    const vector<int>& lambdaRep) 
{
  stringstream fileName;
  fileName << "irrReps/";
  printVec(lambdaRep, fileName);
  createDir(fileName.str());
  fileName << "/" << adjTrans.toStringExact() << ".matrix";
  ifstream ifile(fileName.str());
  if (canRead(ifile)) {
    return readFrom(ifile);
  }
  unsigned fLam = standardTableaux.size();
  if (adjTrans.identity()) {
    Matrix I(fLam);
    return I;
  }
  int k = 1;
  int n = adjTrans.n();
  for (; k < n; ++k) {
    if (adjTrans(k) != k) {
      break;
    }
  }
  //Now we know that adjTrans: k-> k+1 -> k and that's it. 
  Matrix sigmaImg(fLam, fLam);
  vector<bool> seen;
  for (int i=0; i < fLam; ++i) {
    seen.push_back(false);
  }
  for (int i=0; i < fLam; ++i) {
    const LambdaTableau tableau = standardTableaux[i];
    if (tableau.sameRow(k,k+1)) {
      sigmaImg(i,i) = 1;
    } else if (tableau.sameCol(k,k+1)) {
      sigmaImg(i,i) = -1;
    } else if (!seen[i]) {
      //Now we must find its image. 
      LambdaTableau img = adjTrans*tableau;
      int j = i+1;
      for (; j < fLam; ++j) {
        if (img == standardTableaux[j]) {
          break;
        }
      }
      //Compute Axial Distance.
      int d = tableau.axialDistance(k,k+1);
      sigmaImg(i,i) = -(1.0/d);
      sigmaImg(j,i) = sqrt(1-pow(d,-2.0));
      sigmaImg(i,j) = sqrt(1-pow(d,-2.0)); 
      sigmaImg(j,j) = (1.0/d);
      seen[j] = true;
    }
  }
  sigmaImg.roundZero();
  writeToFile(fileName.str(), sigmaImg);
  return sigmaImg;
}

Matrix computeMatrix(const list<SymGpElm>& adjDecomp, const vector<int>& lambdaRep) {
  SymGpElm sigma = prod(adjDecomp);
  stringstream fileName;
  fileName << "irrReps/";
  printVec(lambdaRep, fileName);
  createDir(fileName.str());
  fileName << "/" << sigma.toStringExact() << ".matrix";
  ifstream ifile(fileName.str());
  if (canRead(ifile)) {
    return readFrom(ifile);
  }

  vector<LambdaTableau> standardTableaux = youngTabloids(lambdaRep); 
  unsigned fLam = standardTableaux.size();
  Matrix product(fLam);
  for (const SymGpElm adjTrans : adjDecomp) {
    product *= computeMatrixAdj(adjTrans, standardTableaux, lambdaRep);
    product.roundZero();
  }
  writeToFile(fileName.str(), product);
  return product;
}

int characterIrrDecomp(const list<SymGpElm>& s, const vector<int>& lambda) {
  return computeMatrix(s,lambda).trace();
}

int irrChar(const SymGpElm& s, const vector<int>& lambda) {
  /*TODO: Make this not incredibly stupid. 
   */
  return characterIrrDecomp(decompose(s), lambda);
}

vector<vector<Matrix>> irrReps(int n) {
  vector<vector<Matrix>> reps;
  vector<Matrix> irrRep;
  vector<list<SymGpElm>> Sn = adjTransDecList(n);
  vector<vector<int>> repTypes = nPartitions(n);
  for (vector<int> lambda : repTypes) {
    irrRep.clear();
    for (list<SymGpElm> adjDecomp : Sn) {
      Matrix gpElm = computeMatrix(adjDecomp,lambda);
      irrRep.push_back(gpElm);
    }
    reps.push_back(irrRep);
  }
  return reps;
}

Matrix D(const SymGpElm& s, const vector<int>& lambdaSpace, vector<LambdaTableau>& basis) {
  stringstream fileName;
  fileName << "tabloidReps/";
  printVec(lambdaSpace, fileName);
  createDir(fileName.str());
  fileName << "/" << s.toStringExact() << ".matrix";
  ifstream ifile(fileName.str());
  if (canRead(ifile)) {
    return readFrom(ifile);
  }
  int n = basis.size(); 
  Matrix Ds(n,n);
  //TODO: Maybe not use linear search here? 
  for (size_t i = 0; i < n; ++i) {
    LambdaTableau img = s*basis[i];
    for (size_t j = 0; j < n; ++j) 
      if (img == basis[j])
        Ds(i,j) = 1;
  }
  Ds.roundZero();
  writeToFile(fileName.str(), Ds);
  return Ds;
}

Matrix Pd(const vector<int>& lambdaRep, const vector<int>& lambdaSpace, int n, int mu) {
  vector<list<SymGpElm>> Sn = adjTransDecList(n);
  size_t total = Sn.size();
  size_t count = 1;
  //Takes nontrivial time \/\/
  vector<LambdaTableau> basis = CXlambdaBasis(lambdaSpace);
  unsigned dim = basis.size();
  Matrix sum(dim, dim);
  for (list<SymGpElm> s : Sn) {
    //Progress Bar: 
    cout << "SymGpElm " << count << " of " << total << "|______________\r";
    cout.flush();
    SymGpElm sigma = prod(s);
    Matrix Ds = D(sigma, lambdaSpace, basis);
    cout << "SymGpElm " << count << " of " << total << "|Ds computed___\r";
    cout.flush();
    reverse(s.begin(), s.end()); //Computing s^{-1}
    Matrix dLam = computeMatrix(s, lambdaRep);
    cout << "SymGpElm " << count << " of " << total << "|dLam computed_\r";
    cout.flush();
    sum += dLam(0,mu-1)*Ds;
    count++;
  }
  sum.roundZero();
  return sum;
}

Matrix P(const vector<int>& lambdaRep, const vector<int>& lambdaSpace, int n, int mu) {
  // if g = fac(n) and nj = fLambda(lambdaRep, n), then nj/g =
  // 1/hookProduct(lambdaRep)
  double hookP = hookProduct(lambdaRep);
  return (1/hookP)*Pd(lambdaRep, lambdaSpace, n, mu);
}

Matrix pi(const vector<int>& lambdaRep, const vector<int>& lambdaSpace, int n) {
  stringstream fileName;
  //fileName << "pijs/" << n << "/";
  //createDir(fileName.str());
  //printVec(lambdaSpace, fileName);
  //createDir(fileName.str());
  //fileName << "/" ;
  //printVec(lambdaRep, fileName);
  //fileName << ".matrix";
  //ifstream ifile(fileName.str());
  //if (canRead(ifile))
  //  return readFrom(ifile);
  Matrix pd = Pd(lambdaRep, lambdaSpace, n, 1);
  //writeToFile(fileName.str(),pd);
  return pd;
}

scalar character(const SymGpElm& s, const vector<int>& lambda) {
  //Presumes using representation $D(s,lambda)$!
  vector<LambdaTableau> basis = CXlambdaBasis(lambda);
  return D(s,lambda, basis).trace();
}

int numCycleTypes(int n, const vector<int>& cycleType) {
  int denom = 1;
  vector<int> c;
  for (int i = 0; i < n; ++i) {
    c.push_back(0);
  }
  for (int j : cycleType) {
    c[j] += 1;
    denom *= j;
  }
  for (int i : c) {
    denom *= factorial(i);
  }
  return (factorial(n)/denom);
}

SymGpElm cycType(int n, const vector<int>& cycleType) {
  stringstream sigma;
  vector<int> map;
  for (size_t i = 2; i <= n+1; ++i)
    map.push_back(i);
  size_t sum = 0;
  for (size_t i = 0; i < cycleType.size(); ++i) {
    sum += cycleType[i];
    map[sum-1] = i+1;
  }
  cout << "map = " << map << endl;
  int count = 1;
  for (int i : cycleType) {
    sigma << "(";
    for (int j = 1; j <= i; ++j, ++count) {
      sigma << count;
    }
    sigma << ")";
  }
  cout << "sigma = " << sigma.str() << endl;
  return (SymGpElm(n, sigma.str()));
}

vector<int> charRow(int n, const vector<int>& lambda) {
  vector<vector<int>> cycleTypes = nPartitions(n);
  vector<int> characters;
  for (vector<int> cycleType : cycleTypes) {
    SymGpElm sigma = cycType(n,cycleType);
    characters.push_back(character(sigma,lambda));
  }
  return characters;
}

vector<int> irrCharRow(int n, const vector<int>& lambda) {
  vector<vector<int>> cycleTypes = nPartitions(n);
  vector<int> characters;
  for (vector<int> cycleType : cycleTypes) {
    SymGpElm sigma = cycType(n, cycleType);
    characters.push_back(irrChar(sigma, lambda));
  }
  return characters;
}

vector<int> repDecomp(int n, vector<int> lambda) {
  vector<int> decomp;
  vector<vector<int>> partitions = nPartitions(n);
  vector<vector<int>> cycleTypes = partitions;
  vector<int> repCharRow = charRow(n, lambda);
  int orG = factorial(n);
  for (vector<int> irrRepLambda : partitions) {
    vector<int> irrRepCharRow = irrCharRow(n, irrRepLambda);
    //Inner product: 
    int sum = 0;
    for (size_t i = 0; i < cycleTypes.size(); ++i) {
      vector<int> cycleType = cycleTypes[i];
      int repChar = repCharRow[i];
      int irrRepChar = irrRepCharRow[i];
      sum += repChar*irrRepChar*numCycleTypes(n, cycleType);
    }
    decomp.push_back(sum/orG);
  }
  return decomp;
}

vector<vec> cjSetGS(const Matrix& pij) {
  vector<vec> cjSet;
  Matrix grammedCols = pij.gramSchmidtCols();
  for (int col = 0; col < grammedCols.cols(); ++col) {
    vec column = grammedCols.getCol(col);
    if (!isZero(column)) {
      cjSet.push_back(column);
    }
  }
  return cjSet;
}

vector<vec> cjSetRR(const Matrix& pij) {
  vector<vec> cjSet;
  Matrix reducedPij = pij.rowEchelonForm();
  size_t row = 0;
  for (size_t col = 0; col < reducedPij.cols(); ++col) {
    if (reducedPij(row,col) != 0) {
      cjSet.push_back(pij.getCol(col));
      row++;
    }
  }
  return cjSet;
}

vector<vec> cjSet(const Matrix& pij) {
  return cjSetGS(pij);
}

pair<vector<vector<vec>>, vector<SymGpElm>> cjSetPermuted(Matrix pij, const vector<int>& repType, int nj) {
  int n = pij.cols();
  int count = 1;
  Matrix pijCopy = pij;
  vector<vector<vec>> cjSets;
  vector<SymGpElm> permutations;
  vector<list<SymGpElm>> Sn = adjTransDecList(n);
  for (list<SymGpElm> sigmaDec : Sn) {
    SymGpElm sigma = prod(sigmaDec);

    //string padding = "______________________________________";
    //cout << padding << "repType " << repType << "(" << count << " of " << total;
    //cout << "). Computed pij (1 of "<< nj<<") "; 
    //cout << "Computing Permutation set (" << count << " of " << factorial(n);
    //cout << ")" << padding << "\r";
    //cout.flush();

    bool meaningful = permuteCols(pij, sigma);
    if (meaningful) {
      permutations.push_back(sigma);
      vector<vec> cjs = cjSet(pij);
      cjSets.push_back(cjs);
    }
    ++count;
  }
  return {cjSets, permutations};
}

int sum(const vector<int>& lambda) {
  int sum = 0;
  for (int i : lambda) 
    sum += i;
  return sum;
}

bool contributes(const vector<int>& irrRep, const vector<int>& lambdaSpace) {
  //Note: This presumes that irrRep is on at most \sum(lambdaSpace) elements.
  //Note: This is certainly necessary, but not sufficient. However, its fast
  //      to code and compute and has some positive benefit. 
  //TODO: Make version of nPartitions(n) that takes lambdaSpace and only 
  //      returns those partitions that will contribute. 
  return (irrRep.size() <= lambdaSpace.size());
}

pair<vector<vector<vec>>, vector<vector<SymGpElm>>> allPermutationsBases(int n, const vector<int>& lambdaSpace) {
  vector<vector<vec>> bases;
  vector<vector<SymGpElm>> permutations;
  vector<vector<int>> partitions = nPartitions(n);
  size_t count = 1;
  size_t total = partitions.size();
  //string padding = "______________________________________";
  for (vector<int> repType : partitions) {
    if (!contributes(repType, lambdaSpace)) {
      //cout << "lambdaSpace = " << lambdaSpace << endl;
      //cout << "repType = " << repType << endl;
      continue;
    }
    int nj = fLambda(repType, n);
    //cout << padding << "repType " << repType << "(" << count << " of " << total;
    //cout << "). Computing pij (1 of "<< nj<<")" << padding;
    //cout << "\r";
    //cout.flush();
    Matrix pij = pi(repType, lambdaSpace, n);
    pair<vector<vector<vec>>, vector<SymGpElm>> cjPermuted = cjSetPermuted(pij, repType, nj);
    permutations.push_back(cjPermuted.second);
    vector<Matrix> Pmats;
    for (size_t k = 2; k <= nj; ++k) {
      //cout << padding << "repType " << repType << "(" << count << " of " << total;
      //cout << "). Computing Pmat (" << k << " of "<< nj<<")" << padding;
      //cout << "\r";
      //cout.flush();
      Pmats.push_back(P(repType, lambdaSpace, n, k));
    }
    vector<vector<vec>> cjs = cjPermuted.first;
    for (vector<vec> cj : cjs) {
      vector<vec> Bfinal;
      for (vec vi : cj) {
        Bfinal.push_back(vi);
        for (size_t k = 0; k <= nj-2; ++k) {
          Bfinal.push_back(Pmats[k]*vi);
        }
      }
      bases.push_back(Bfinal);
    }
    count++;
  }
  return {bases,permutations};
}

vector<vec> finalBasis(int n, const vector<int>& lambdaSpace) {
  vector<vec> Bfinal;
  vector<vector<int>> partitions = nPartitions(n);
  size_t count = 1;
  size_t total = partitions.size();
  string padding = "______________________________________";
  for (vector<int> repType : partitions) {
    if (!contributes(repType, lambdaSpace)) {
      //cout << "lambdaSpace = " << lambdaSpace << endl;
      //cout << "repType = " << repType << endl;
      continue;
    }
    int nj = fLambda(repType, n);
    cout << padding << "repType " << repType << "(" << count << " of " << total;
    cout << "). Computing pij (1 of "<< nj<<")" << padding;
    cout << "\r";
    cout.flush();
    Matrix pij = pi(repType, lambdaSpace, n);
    vector<vec> cj = cjSet(pij);
    if (!cj.empty()) {
      for (vec vi : cj) {
        Bfinal.push_back(vi);
      }
      for (size_t k = 2; k <= nj; ++k) {
        cout << padding << "repType " << repType << "(" << count << " of " << total;
        cout << "). Computing Pmat (" << k << " of "<< nj<<")" << padding;
        cout << "\r";
        cout.flush();
        Matrix Pmat = P(repType, lambdaSpace, n, k);
        for (vec vi: cj) {
          Bfinal.push_back(Pmat*vi);
        }
      }
    } else {
      cout << endl << endl << "repType: " << repType << endl;
      cout << "lambdaSpace: " << lambdaSpace << endl;
      cout << "contributes(repType,lambdaSpace): ";
      cout << contributes(repType, lambdaSpace) << endl << endl;
    }
    count++;
  }
  return Bfinal;
}

Matrix COBmatrix(const vector<vec>& Bf) {
  size_t dim = Bf.size();
  Matrix BfTB1(dim,dim);
  for (size_t i = 0; i < dim; ++i) {
    BfTB1.setCol(i, Bf[i]);
  }
  BfTB1.roundZero();
  return BfTB1;
}

Matrix COBmatrix(const vector<vec>& BStart, const vector<vec>& BEnd) {
  //B1TBStart
  Matrix BStartTB1 = COBmatrix(BStart);
  //B1TBEnd
  Matrix BEndTB1 = COBmatrix(BEnd);
  Matrix B1TBEnd = BEndTB1.inverse();
  return B1TBEnd*BStartTB1;
}

void findBasisDecomps(const string& filePath, const string& fileName, const int n, const vector<int>& lambdaSpace) {
  stringstream totalFileName;
  totalFileName << filePath << "/" << fileName;
  ofstream file(totalFileName.str());
  file << "n: " << n << " lambdaSpace: " << lambdaSpace << endl;
  Matrix cobMatrix;
  vector<vec> Bprev, Bcurr;
  Bcurr = finalBasis(1,lambdaSpace);
  for (int i = 2; i <= n; ++i) {
    stringstream cobFileName;
    cobFileName << "S" << i-1 << "->S" << i << ".matrix";
    stringstream cobPath;
    cobPath << filePath << "/" << cobFileName.str();
    cout << "cobPath.string: " << cobPath.str() << endl;
    ifstream iCobFile(cobPath.str());
    if (iCobFile) {
      cout << "iCobFile found!" << endl;
      cobMatrix = readFrom(iCobFile);
    } else {
      Bprev = Bcurr;
      Bcurr = finalBasis(i, lambdaSpace);

      cobMatrix = COBmatrix(Bprev, Bcurr);
      ofstream oCobFile(cobPath.str());
      writeTo(oCobFile, cobMatrix);
    }

    file << endl << "B" << i-1 << " -> B" << i << ": " << endl;
    cobMatrix.prettyPrintTo(file);
    file << "(# non-zero entries)/dimension = " << cobMatrix.numNonzeroEntries();
    file << "/" << cobMatrix.rows() << " = ";
    file << (1.0*cobMatrix.numNonzeroEntries())/cobMatrix.rows() << endl;
    //TODO: Maybe doing extra work here, don't need to write twice?
  }
}

void findAllPermutedBases(const string& filePath, const int n, const vector<int>& lambdaSpace) {
  /* In this function, we find all possible finalBases that arise via changing
   * the ordering of Gram Schmidt during computation for all $1 < k \le n$. 
   */
  stringstream totalFileName;
  totalFileName << filePath << "/";
  vector<vector<vec>> bases;
  vector<vector<SymGpElm>> permutations;
  pair<vector<vector<vec>>,vector<vector<SymGpElm>>> res;
  for (int i = 2; i <= n; ++i) {
    res = allPermutationsBases(i, lambdaSpace);
    bases = res.first;
    permutations = res.second;
    totalFileName << "S" << i << "/";
    createDir(totalFileName.str());
    for (size_t j = 0; j < bases.size(); ++j) {
      totalFileName << "Basis" << j+1;
      vector<vec> basis = bases[j];
      vector<SymGpElm> permSet = permutations[j];
      Matrix cobMatrix = COBmatrix(basis);
      ofstream oBasisFile(totalFileName.str());
      printVec(permSet, oBasisFile);
      writeTo(oBasisFile, cobMatrix);
    }
  }
}

void test() {
  int n = 3;
  Matrix testM(n);
  testM.prettyPrint();
  SymGpElm sigma(3, "(123)");
  bool hi = permuteCols(testM, sigma);
  testM.prettyPrint();
}

int main(int argc, const char* argv[]) {
  int n;
  vector<int> lambdaSpace;
  if (argc == 2) {
    n = atoi(argv[1]);
    lambdaSpace = {n-2, 1, 1};
  } else if (argc > 2) {
    n = atoi(argv[1]);
    lambdaSpace.push_back(n);
    int count = 0;
    for (int i = 2; i < argc; ++i) {
      int row = atoi(argv[i]);
      if (row > lambdaSpace[i-2]) {
        cout << "Tabloid Shape Not Allowed!" << endl;
        return 1;
      }
      count += row;
      lambdaSpace.push_back(row);
    }
    lambdaSpace[0] = n - count;
    if (lambdaSpace[0] < lambdaSpace[1]) {
      cout << "Tabloid Shape Not Allowed!" << endl;
      printVec(lambdaSpace);
      cout << endl;
      return 1;
    }
  } else {
    n = 5;
    lambdaSpace = {3, 1, 1};
  }
  stringstream ss;
  ss << "factorizations/";
  printVec(lambdaSpace, ss);
  const string filePath = ss.str();
  createDir(filePath);
  ss.str("");
  ss << "matricies-S" << n << ".txt";
  const string fileName = ss.str();
  //test();
  //findBasisDecomps(filePath, fileName, n, lambdaSpace);
  findAllPermutedBases(filePath,  n, lambdaSpace);
  cout << endl;
  return 0;
}
