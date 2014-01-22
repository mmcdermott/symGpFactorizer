#include <iostream>
#include <vector>
#include <list> 

template <class T>
void printVec(const std::vector<T>& v, std::ostream& out = std::cout) {
  if (v.size() == 0) {
    out << "[]";
    return;
  }
  typename std::vector<T>::const_iterator it = v.begin();
  out << "[" << *it;
  ++it;
  for (; it != v.end(); ++it)
    out << ", " << *it;
  out << "]";
}

template <class T> 
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  printVec(v,out);
  return out;
}

template <class T> 
//typename std::list<T>::iterator lIterT;
void printList(const std::list<T>& L, std::ostream& out = std::cout) {
  typename std::list<T>::const_iterator it = L.begin();
  out << "(" << *it;
  for (; it != L.end(); ++it) {
    out << ", " << *it;
  }
  out << ")";
}

template <class T>
std::ostream& operator<<(std::ostream& out, const std::list<T>& L) {
  printList(L,out);
  return out;
}
