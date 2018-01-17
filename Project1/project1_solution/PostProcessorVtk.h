// -*- C++ -*-
#ifndef POST_PROCESSOR_VTK_H
#define POST_PROCESSOR_VTK_H

#include "Definitions.h"

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCubeSource.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkVersion.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkProperty.h>
#include <vtkTriangle.h>

// #include "Utilities.h"

namespace PostProcessors {

enum Configuration {UndeformedConfiguration, DeformedConfiguration};

namespace Vtk {


class EmptyElement{
  
  public:
  
  typedef Matrix2d Stress;
  typedef Matrix2d Strain;
  typedef Matrix2d Vector;
  
  static const unsigned int NumberOfNodes    = 3;
  static const unsigned int SpatialDimension = 2;
  static const unsigned int VtkCellType      = VTK_TRIANGLE;
  
  typedef Vector2d Point;
  
  EmptyElement(){}
  
};

// =============================================================================
// *********************< Machinery for NamedArrays >***************************
// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

enum ElementWiseOrNodeWise {ElementWise, NodeWise};

template <class T>
struct NamedArray {
  vector<T>             _array;
  string                _title;
  ElementWiseOrNodeWise _elementWiseOrNodeWise;
  vector<string>        _componentTitles;
};

// default implementation is for eigen objects
template <class T>
unsigned int
getNumberOfComponents(const T & t) {
  return t.size();
}

// special implementations for standard types
template <> unsigned int getNumberOfComponents<int>(const int & t) {
  ignoreUnusedVariable<int>(t);
  return 1;
}
template <> unsigned int getNumberOfComponents<double>(const double & t) {
  ignoreUnusedVariable<double>(t);
  return 1;
}
template <> unsigned int getNumberOfComponents<float>(const float & t) {
  ignoreUnusedVariable<float>(t);
  return 1;
}
template <> unsigned int getNumberOfComponents<unsigned int>(const unsigned int & t) {
  ignoreUnusedVariable<unsigned int>(t);
  return 1;
}
template <> unsigned int getNumberOfComponents<long unsigned int>(const long unsigned int & t) {
  ignoreUnusedVariable<long unsigned int>(t);
  return 1;
}

// default implementation is for eigen objects
template <class T>
void
copyValueToDock(const T & t, vector<double> * dock) {
  for (unsigned int i = 0; i < t.size(); ++i) {
    dock->at(i) = t(i);
  }
}

// special implementations for standard types
template <> void copyValueToDock<double>(const double & t, vector<double> * dock) {dock->at(0) = t;}
template <> void copyValueToDock<float>(const float & t, vector<double> * dock) {dock->at(0) = t;}
template <> void copyValueToDock<int>(const int & t, vector<double> * dock) {dock->at(0) = t;}
template <> void copyValueToDock<unsigned int>(const unsigned int & t, vector<double> * dock) {dock->at(0) = t;}
template <> void copyValueToDock<long unsigned int>(const long unsigned int & t, vector<double> * dock) {dock->at(0) = t;}

template <class T, class VtkDataPointer>
void
addNamedArrayToVtkData(const NamedArray<T> & namedArray,
                       VtkDataPointer vtkData) {

  if (namedArray._array.size() == 0) {
    fprintf(stderr, "cannot add vtkNamedArray named %s to vtk data, it "
            "has %zu entries\n", namedArray._title.c_str(),
            namedArray._array.size());
    exit(1);
  }
  if (namedArray._elementWiseOrNodeWise == ElementWise &&
      namedArray._array.size() != vtkData->GetNumberOfCells()) {
    fprintf(stderr, "cannot add vtkNamedArray named %s to vtk data cells, it "
            "has %zu entries and there are %zd cells\n",
            namedArray._title.c_str(), namedArray._array.size(),
            size_t(vtkData->GetNumberOfCells()));
    exit(1);
  }
  if (namedArray._elementWiseOrNodeWise == NodeWise &&
      namedArray._array.size() != vtkData->GetNumberOfPoints()) {
    fprintf(stderr, "cannot add vtkNamedArray named %s to vtk data points, it "
            "has %zu entries and there are %zd points\n",
            namedArray._title.c_str(), namedArray._array.size(),
            size_t(vtkData->GetNumberOfPoints()));
    exit(1);
  }

  const unsigned int numberOfComponents =
    getNumberOfComponents<T>(namedArray._array[0]);

  if (namedArray._componentTitles.size() > 0 &&
      namedArray._componentTitles.size() != numberOfComponents) {
    fprintf(stderr, "cannot add vtkNamedArray named %s to vtk data, the data "
            "type has %u components but the _componentTitles.size() = %zu\n",namedArray._title.c_str(),
            numberOfComponents, namedArray._componentTitles.size());
  }

  vtkSmartPointer<vtkDoubleArray> vtkArray =
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkArray->SetName(namedArray._title.c_str());
  vtkArray->SetNumberOfComponents(numberOfComponents);
  if (namedArray._componentTitles.size() > 0) {
    for (unsigned int componentIndex = 0; componentIndex < numberOfComponents;
         ++componentIndex) {
      vtkArray->SetComponentName(componentIndex,
                                 namedArray._componentTitles[componentIndex].c_str());
    }
  }
  vector<double> dock;
  dock.resize(numberOfComponents);
  for (unsigned int index = 0; index < namedArray._array.size(); ++index) {
    copyValueToDock<T>(namedArray._array[index], &dock);
    vtkArray->InsertNextTuple(&dock[0]);
  }

  if (namedArray._elementWiseOrNodeWise == ElementWise) {
    vtkData->GetCellData()->AddArray(vtkArray);
  } else {
    vtkData->GetPointData()->AddArray(vtkArray);
  }
}

template <class ... NamedArrayTypes_ButIReallyCanPutAnythingIWantHereSoWhyDontIDeclareMyLoveForMyWifeAndKidsThatsNiceIsntIt_Jeff>
class NamedArrays {
};

template <>
struct NamedArrays <> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    ignoreUnusedVariable<VtkDataType>(data);
  }

};

template <class T0>
struct NamedArrays <T0> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
};

template <class T0, class T1>
struct NamedArrays <T0, T1> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
};

template <class T0, class T1, class T2>
struct NamedArrays <T0, T1, T2> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
};

template <class T0, class T1, class T2, class T3>
struct NamedArrays <T0, T1, T2, T3> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T3> & namedArray : _namedArrays3) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
  vector<NamedArray<T3> > _namedArrays3;
};

template <class T0, class T1, class T2, class T3, class T4>
struct NamedArrays <T0, T1, T2, T3, T4> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T3> & namedArray : _namedArrays3) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T4> & namedArray : _namedArrays4) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T3> & namedArray) {
    _namedArrays3.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T4> & namedArray) {
    _namedArrays4.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
  vector<NamedArray<T3> > _namedArrays3;
  vector<NamedArray<T4> > _namedArrays4;
};

template <class T0, class T1, class T2, class T3, class T4, class T5>
struct NamedArrays <T0, T1, T2, T3, T4, T5> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T3> & namedArray : _namedArrays3) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T4> & namedArray : _namedArrays4) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T5> & namedArray : _namedArrays5) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T3> & namedArray) {
    _namedArrays3.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T4> & namedArray) {
    _namedArrays4.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T5> & namedArray) {
    _namedArrays5.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
  vector<NamedArray<T3> > _namedArrays3;
  vector<NamedArray<T4> > _namedArrays4;
  vector<NamedArray<T5> > _namedArrays5;
};

template <class T0, class T1, class T2, class T3, class T4, class T5, class T6>
struct NamedArrays <T0, T1, T2, T3, T4, T5, T6> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T3> & namedArray : _namedArrays3) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T4> & namedArray : _namedArrays4) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T5> & namedArray : _namedArrays5) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T6> & namedArray : _namedArrays6) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T3> & namedArray) {
    _namedArrays3.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T4> & namedArray) {
    _namedArrays4.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T5> & namedArray) {
    _namedArrays5.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T6> & namedArray) {
    _namedArrays6.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
  vector<NamedArray<T3> > _namedArrays3;
  vector<NamedArray<T4> > _namedArrays4;
  vector<NamedArray<T5> > _namedArrays5;
  vector<NamedArray<T6> > _namedArrays6;
};

template <class T0, class T1, class T2, class T3, class T4, class T5, class T6, class T7>
struct NamedArrays <T0, T1, T2, T3, T4, T5, T6, T7> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T3> & namedArray : _namedArrays3) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T4> & namedArray : _namedArrays4) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T5> & namedArray : _namedArrays5) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T6> & namedArray : _namedArrays6) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T7> & namedArray : _namedArrays7) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T3> & namedArray) {
    _namedArrays3.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T4> & namedArray) {
    _namedArrays4.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T5> & namedArray) {
    _namedArrays5.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T6> & namedArray) {
    _namedArrays6.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T7> & namedArray) {
    _namedArrays7.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
  vector<NamedArray<T3> > _namedArrays3;
  vector<NamedArray<T4> > _namedArrays4;
  vector<NamedArray<T5> > _namedArrays5;
  vector<NamedArray<T6> > _namedArrays6;
  vector<NamedArray<T7> > _namedArrays7;
};

template <class T0, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8>
struct NamedArrays <T0, T1, T2, T3, T4, T5, T6, T7, T8> {

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T3> & namedArray : _namedArrays3) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T4> & namedArray : _namedArrays4) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T5> & namedArray : _namedArrays5) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T6> & namedArray : _namedArrays6) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T7> & namedArray : _namedArrays7) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T8> & namedArray : _namedArrays8) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T3> & namedArray) {
    _namedArrays3.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T4> & namedArray) {
    _namedArrays4.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T5> & namedArray) {
    _namedArrays5.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T6> & namedArray) {
    _namedArrays6.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T7> & namedArray) {
    _namedArrays7.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T8> & namedArray) {
    _namedArrays8.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
  vector<NamedArray<T3> > _namedArrays3;
  vector<NamedArray<T4> > _namedArrays4;
  vector<NamedArray<T5> > _namedArrays5;
  vector<NamedArray<T6> > _namedArrays6;
  vector<NamedArray<T7> > _namedArrays7;
  vector<NamedArray<T8> > _namedArrays8;
};

template <class T0, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9>
struct NamedArrays <T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>{

  NamedArrays() {
  }

  template <class VtkDataType>
  void
  addArraysToVtkData(VtkDataType data) const {
    for (const NamedArray<T0> & namedArray : _namedArrays0) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T1> & namedArray : _namedArrays1) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T2> & namedArray : _namedArrays2) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T3> & namedArray : _namedArrays3) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T4> & namedArray : _namedArrays4) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T5> & namedArray : _namedArrays5) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T6> & namedArray : _namedArrays6) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T7> & namedArray : _namedArrays7) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T8> & namedArray : _namedArrays8) addNamedArrayToVtkData(namedArray, data);
    for (const NamedArray<T9> & namedArray : _namedArrays9) addNamedArrayToVtkData(namedArray, data);
  }

  void
  addArray(const NamedArray<T0> & namedArray) {
    _namedArrays0.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T1> & namedArray) {
    _namedArrays1.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T2> & namedArray) {
    _namedArrays2.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T3> & namedArray) {
    _namedArrays3.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T4> & namedArray) {
    _namedArrays4.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T5> & namedArray) {
    _namedArrays5.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T6> & namedArray) {
    _namedArrays6.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T7> & namedArray) {
    _namedArrays7.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T8> & namedArray) {
    _namedArrays8.push_back(namedArray);
  }

  void
  addArray(const NamedArray<T9> & namedArray) {
    _namedArrays9.push_back(namedArray);
  }

private:
  vector<NamedArray<T0> > _namedArrays0;
  vector<NamedArray<T1> > _namedArrays1;
  vector<NamedArray<T2> > _namedArrays2;
  vector<NamedArray<T3> > _namedArrays3;
  vector<NamedArray<T4> > _namedArrays4;
  vector<NamedArray<T5> > _namedArrays5;
  vector<NamedArray<T6> > _namedArrays6;
  vector<NamedArray<T7> > _namedArrays7;
  vector<NamedArray<T8> > _namedArrays8;
  vector<NamedArray<T9> > _namedArrays9;
};

// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// *********************< Machinery for NamedArrays >***************************
// =============================================================================




// =============================================================================
// ***************************< General functions >*****************************
// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

template <class Element, class NamedArrays>
vtkSmartPointer<vtkUnstructuredGrid>
makeUnstructuredGrid(const Mesh &      mesh        ,
                     const NamedArrays &  namedArrays ) {

  typedef typename Element::Point Point;

  vtkSmartPointer<vtkPoints> points
    = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkIntArray> elementIds
    = vtkSmartPointer<vtkIntArray>::New();
  elementIds->SetName("ElementId");

  std::array<double, 3> pointStaticArray;
  pointStaticArray[0] = 0;
  pointStaticArray[1] = 0;
  pointStaticArray[2] = 0;
  for (unsigned int nodeIndex = 0; nodeIndex < mesh._nodes.size(); ++nodeIndex) {
    const Point & undeformedPosition = mesh._nodes[nodeIndex]._position;

    const Point displayPoint = undeformedPosition;
    for (unsigned int coordinate = 0; coordinate < Element::SpatialDimension; ++coordinate) {
      pointStaticArray[coordinate] = displayPoint[coordinate];
    }
    points->InsertNextPoint(&pointStaticArray[0]);
  }

  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  unstructuredGrid->SetPoints(points);

  std::array<vtkIdType, Element::NumberOfNodes> pointIds;
  for (size_t elementIndex = 0; elementIndex < mesh._connectivity.size(); ++elementIndex) {
    const array<int, Element::NumberOfNodes> & nodeIds
      = mesh._connectivity[elementIndex];
    for (size_t nodeNumber = 0; nodeNumber < Element::NumberOfNodes; ++nodeNumber) {
      pointIds[nodeNumber] = nodeIds[nodeNumber];
    }
    unstructuredGrid->InsertNextCell(Element::VtkCellType,
                                     Element::NumberOfNodes,
                                     &pointIds[0]);
    elementIds->InsertNextValue(elementIndex);
  }

  unstructuredGrid->GetCellData()->AddArray(elementIds);

  namedArrays.addArraysToVtkData(unstructuredGrid);

  return unstructuredGrid;
}


template <class NamedArrays>
void
makeMeshFile(const Mesh &         mesh        ,
             const string &       filename    ,
             const NamedArrays &  namedArrays ) {

  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    makeUnstructuredGrid<EmptyElement>(mesh, namedArrays);

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  string appendedOutputFilename = filename;
  appendedOutputFilename += ".vtu";
  writer->SetFileName(appendedOutputFilename.c_str());
  writer->SetInputData(unstructuredGrid);

  writer->SetDataModeToBinary();
  int result = writer->Write();
  if (result != 1) {
    printf("Failed to write vtk file with filename %s",
           appendedOutputFilename.c_str());
  }
}

// undeformed, no displacements, no stresses, no boundary conditions (0 0 0)
template <class NamedArrays = PostProcessors::Vtk::NamedArrays<> >
void
makeUndeformedMeshFile(const Mesh &    mesh                        ,
                       const string &     filename                    ,
                       const NamedArrays  namedArrays = NamedArrays() ) {
  makeMeshFile<NamedArrays>(mesh, filename, namedArrays);
}

} // Vtk

} // PostProcessors


#endif  // POST_PROCESSOR_VTK_H
