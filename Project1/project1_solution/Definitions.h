// -*- C++ -*-

// IFNDEF-DEFINE-ENDIF-Construct allows to prevent that the same library is loaded twice
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// Basic math, I/O and other fundamental includes
#include <cmath>
#include <cstdio>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <set>
#include <cstdarg>

// Linear Algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#pragma GCC diagnostic pop
#define EIGEN_SUPERLU_SUPPORT

// Postprocessing
#include <vtkCellType.h>

// Default error statement definition
// from here: http://stackoverflow.com/questions/6420194/how-to-print-both-to-stdout-and-file-in-c
#define errorStatement(s, ...)                                  \
  do {                                                          \
    fprintf (stderr, "(%30s:%40s:%4d) -- " s,                   \
             __FILE__, __func__, __LINE__, ##__VA_ARGS__);      \
    fflush (stderr);                                            \
  } while (0)

#define FILENAMEWITHEXTENSION (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

// Default exception statement
char exceptionBuffer[10000];
#define throwException(s, ...)                  \
  sprintf(exceptionBuffer, s, ##__VA_ARGS__);   \
  throw std::runtime_error(exceptionBuffer);

// Switch Output On/Off
enum VerboseFlag {Quiet, Verbose};
  
// Allows to use array, vector, etc. without having to prepend by std:
using std::array;
using std::vector;
using std::string;
using std::ifstream;
using std::cout;
using std::endl;

// Allows to use Matrix<.,.,.> without having to prepend by Eigen:
using Eigen::Matrix;

typedef Matrix<double, 1, 1> Vector1d;
typedef Matrix<double, 2, 1> Vector2d;
typedef Matrix<double, 3, 1> Vector3d;
typedef Matrix<double, 4, 1> Vector4d;
typedef Matrix<double, 1, 1> Matrix1d;
typedef Matrix<double, 2, 2> Matrix2d;
typedef Matrix<double, 3, 3> Matrix3d;
typedef Matrix<double, 4, 4> Matrix4d;

// Avoid unused variable warnings from the compiler with -Wall on
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template <typename Variable>
void ignoreUnusedVariable(Variable dummy) {
}
template <typename T>
void ignoreUnusedVariables(const T & t) {
}
template <typename T, typename U>
void ignoreUnusedVariables(const T & t, const U & u) {
}
template <typename T, typename U, typename V>
void ignoreUnusedVariables(const T & t, const U & u, const V & v) {
}
template <typename T, typename U, typename V, typename W>
void ignoreUnusedVariables(const T & t, const U & u, const V & v, const W & w) {
}
#pragma GCC diagnostic pop









class Node {
  
  // Public members
  public:
  
  // Constructor
  Node (int ID, Vector2d position){
    _ID       = ID;
    _position = position;
  }
  
  int       _ID;
  Vector2d  _position;
  
};

class Mesh {

  // Public members
  public:
  vector<Node>         _nodes;
  vector<array<int,3>> _connectivity;
  
};



#endif // DEFINITIONS_H
