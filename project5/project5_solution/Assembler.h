#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "Definitions.h"
#include "Utilities.h"

template <class E>
class Assembler {

public:

  typedef E                                              Element;
  typedef typename Element::Vector                       ElementVector;
  typedef typename Element::Forces                       ElementForces;
  typedef typename Element::StiffnessMatrix              ElementStiffnessMatrix;
  typedef typename Element::Stress                       ElementStress;
  static const unsigned int NumberOfNodesPerElement =    Element::NumberOfNodes;
  static const unsigned int SpatialDimension =           Element::SpatialDimension;
  static const unsigned int DegreesOfFreedom =           Element::DegreesOfFreedom;


  Assembler () :_numberOfNodes(0) {
  }
  
  Assembler ( const vector<Element> & elements      ,
              const size_t            numberOfNodes ) :
      _elements(elements),
      _numberOfNodes(numberOfNodes) {
  }
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // All the following functions are given, nothing needs to be changed
  // They're the solution to the last problem set. Feel free to use your
  // own, if you think you did a better job in implementing everything :)
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                  Assemble energy               %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  double
  assembleEnergy(const vector<ElementVector> & displacements) const {
    double energy = 0.0;
    // normal elements
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) {
      const Element & element = _elements[elementIndex];
      array<size_t, NumberOfNodesPerElement> elementNodeIds =
        element.getNodeIds();
      array<ElementVector, NumberOfNodesPerElement> elementDisplacements =
        Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);
      energy += element.computeEnergy(elementDisplacements);
    }
    return energy;
  }

  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                  Assemble forces               %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  
  Eigen::VectorXd
  assembleForceVector(const vector<ElementVector> & displacements) const {
    size_t numberOfDofs = displacements.size() * DegreesOfFreedom;
    Eigen::VectorXd forceVector(numberOfDofs);
    forceVector.fill(0);

    // normal elements
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) {
      const Element & element = _elements[elementIndex];
      array<size_t, NumberOfNodesPerElement> elementNodeIds =
        element.getNodeIds();
      array<ElementVector, NumberOfNodesPerElement> elementDisplacements =
        Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);
      ElementForces elementForces = element.computeForces(elementDisplacements);
      for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size();++nodeIndex) {
        size_t nodeId = elementNodeIds[nodeIndex];
        for (size_t i = 0; i < DegreesOfFreedom; ++i) {
          forceVector(nodeId * DegreesOfFreedom + i) +=
            elementForces[nodeIndex](i);
        }
      }
    }

    return forceVector;
  }
  
  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%               Assemble stiffness               %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%  
  
  Eigen::MatrixXd
  assembleStiffnessMatrix(const vector<ElementVector> & displacements) const {
    size_t numberOfDofs = displacements.size() * DegreesOfFreedom;
    Eigen::MatrixXd stiffnessMatrix(numberOfDofs, numberOfDofs);
    stiffnessMatrix.fill(0);

    // normal elements
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) {
      const Element & element = _elements[elementIndex];
      array<size_t, NumberOfNodesPerElement> elementNodeIds =
        element.getNodeIds();
      array<ElementVector, NumberOfNodesPerElement> elementDisplacements =
        Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);
      ElementStiffnessMatrix elementStiffnessMatrix =
        element.computeStiffnessMatrix(elementDisplacements);
      for (size_t nodeIndex1 = 0; nodeIndex1 < elementNodeIds.size(); ++nodeIndex1) {
        size_t nodeId1 = elementNodeIds[nodeIndex1];
        for (size_t nodeIndex2 = 0; nodeIndex2 < elementNodeIds.size(); ++nodeIndex2) {
          size_t nodeId2 = elementNodeIds[nodeIndex2];
          for (size_t i = 0; i < DegreesOfFreedom; ++i) {
            for (size_t j = 0; j < DegreesOfFreedom; ++j) {
              // heaven help me if this is wrong
              stiffnessMatrix(nodeId1 * DegreesOfFreedom + i,
                              nodeId2 * DegreesOfFreedom + j) +=
                elementStiffnessMatrix(nodeIndex1 * DegreesOfFreedom + i,
                                       nodeIndex2 * DegreesOfFreedom + j);
            }
          }
        }
      }
    }
    return stiffnessMatrix;
  };

  
  
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // All the following functions are given, nothing needs to be changed
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  

  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%             Compute nodal stresses             %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%  
  
  vector<ElementStress>
  computeNodalStresses(const vector<ElementVector> & displacements) const {

    const size_t numberOfNodes = displacements.size();

    vector<ElementStress> nodalStresses(numberOfNodes, ElementStress::Zero());
    vector<double> volumeSums(numberOfNodes, 0);

    vector<ElementStress> elementStresses = computeElementStresses(displacements);

    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) {
      const Element & element = _elements[elementIndex];
      array<size_t, NumberOfNodesPerElement> elementNodeIds = element.getNodeIds();
      array<double, NumberOfNodesPerElement> elementWeights =
        element.computeNodalWeights();
      const ElementStress & elementStress = elementStresses[elementIndex];

      for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size(); nodeIndex++) {
        nodalStresses[elementNodeIds[nodeIndex]] += elementStress * elementWeights[nodeIndex];
        volumeSums[elementNodeIds[nodeIndex]] += elementWeights[nodeIndex];
      }
    }

    for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++) {
      nodalStresses[nodeIndex] /= volumeSums[nodeIndex];
    }
    return nodalStresses;
  }

  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%             Compute nodal strains              %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  
  vector<typename Element::Strain>
  computeNodalStrains(const vector<typename Element::Vector> & displacements) const {

    const size_t numberOfNodes = displacements.size();

    vector<typename Element::Strain> nodalStrains(numberOfNodes,
                                                  Element::Strain::Zero());
    vector<double> volumeSums(numberOfNodes, 0);

    vector<typename Element::Strain> elementStrains =
      computeElementStrains(displacements);

    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) {
      const Element & element = _elements[elementIndex];
      array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
      array<double, Element::NumberOfNodes> elementWeights =
        element.computeNodalWeights();
      const typename Element::Strain & elementStrain = elementStrains[elementIndex];

      for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size(); nodeIndex++) {
        nodalStrains[elementNodeIds[nodeIndex]] +=
          elementStrain / elementWeights[nodeIndex];
        volumeSums[elementNodeIds[nodeIndex]] += 1./elementWeights[nodeIndex];
      }
    }

    for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
      nodalStrains[nodeIndex] /= volumeSums[nodeIndex];
    }
    return nodalStrains;
  }

  
  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%            Compute element stresses            %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%  

  vector<ElementStress>
  computeElementStresses(const vector<ElementVector> & displacements) const {

    vector<ElementStress> allElementStresses(_elements.size());
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) {
      const Element & element = _elements[elementIndex];
      array<size_t, NumberOfNodesPerElement> elementNodeIds = element.getNodeIds();
      array<ElementVector, NumberOfNodesPerElement> elementDisplacements =
        Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,displacements);
      array<ElementStress, Element::QuadPoints> elementStresses =
        element.computeStressesAtGaussPoints(elementDisplacements);
      ElementStress average = ElementStress::Zero();
      for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
        average += elementStresses[qpIndex];
      }
      average /= Element::QuadPoints;
      allElementStresses[elementIndex] = average;
    }
    return allElementStresses;
  }

  
  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%            Compute element strains             %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%  
  
  vector<typename Element::Strain>
  computeElementStrains(const vector<typename Element::Vector> & displacements) const {

    vector<typename Element::Strain> allElementStrains(_elements.size());
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) {
      const Element & element = _elements[elementIndex];
      array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
      array<typename Element::Vector, Element::NumberOfNodes> elementDisplacements =
        ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                    displacements);
      array<typename Element::Strain, Element::QuadPoints> elementStrains =
        element.computeStrainsAtGaussPoints(elementDisplacements);
      typename Element::Strain average = Element::Strain::Zero();
      for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
        average += elementStrains[qpIndex];
      }
      average /= Element::QuadPoints;
      allElementStrains[elementIndex] = average;
    }
    return allElementStrains;
  }

  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%            Get-Functions             %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                      %%%%%%%%%%%%%%%%%%%%  
  
  size_t
  getNumberOfElements() const {
    return _elements.size();
  }

  size_t
  getNumberOfNodes() const {
    return _numberOfNodes;
  }

private:

  vector<Element> _elements;
  const size_t _numberOfNodes;
};

#endif  // ASSEMBLER_H
