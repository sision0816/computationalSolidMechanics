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

  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                  Assemble energy               %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  
  
  double
  assembleEnergy(const vector<ElementVector> & displacements) const {
    
    double energy = 0.0;
    
    // ToDo:
    // 1) Loop through all elements (Note: there are _elements.size() elements)
    // 2) Create an Element called element corresponding to the elementIndexth element from _elements
    // 3) Get the Node-IDs of that element - Note: the function getNodeIds() is provided by the element
    //    and returns an array<size_t, NumberOfNodesPerElement>
    // 4) Find the element displacements of the element using the Utility-function:
    //    Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds, displacements);
    //    Why? displacements is a vector which has a displacement for every node. We want an array
    //    with NumberOfNodesPerElement ElementVector's that corresponds to the displacement at the
    //    nodes with elementNodeIds.
    //    This also tells you the which array<...,...> you have to use ... :)
    // 5) Use the element displacement to calculate the energy in every element and add them together
    //    (Remind yourself of the computeEnergy function on the Element level)
    
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) { // (Step 1 - already given)
		Element element = _elements[elementIndex];
	    array<size_t, NumberOfNodesPerElement> elementNodeIds = element.getNodeIds();
		array<ElementVector, NumberOfNodesPerElement> elementDisplacements = Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds, displacements);
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

    // ToDo:
    // 1) Loop through all elements (Note: there are _elements.size() elements)
    // 2) Create an Element called element corresponding to the elementIndexth element from _elements
    // 3) Get the Node-IDs of that element - Note: the function getNodeIds() is provided by the element
    //    and returns an array<size_t, NumberOfNodesPerElement>
    // 4) Find the element displacements of the element using the Utility-function:
    //    Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds, displacements);
    //    Why? displacements is a vector which has a displacement for every node. We want an array
    //    with NumberOfNodesPerElement ElementVector's that corresponds to the displacement at the
    //    nodes with elementNodeIds.
    //    This also tells you the which array<...,...> you have to use ... :)
    // 5) Evaluate the element forces (of type ElementForces)
    //    (Remind yourself of the computeForce functionality on the Element level)
    // 6) Create the force vector by adding all the elementforces together. 
    //    Pay attention on where each contribution is located in the final global force vector!!!
    //
    // NOTE: The for-loops outlined simply reflect our recommendations. You may be able to solve it with
    //       less or with more. So long as the result is not hard-coded, please feel free to append or
    //       erase for-loops to your liking :)
    
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) { // (Step 1 - given)
      
        Element element = _elements[elementIndex];
	    array<size_t, NumberOfNodesPerElement> elementNodeIds = element.getNodeIds();
		array<ElementVector, NumberOfNodesPerElement> elementDisplacements = Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds, displacements);
		ElementForces elementForces = element.computeForces(elementDisplacements);
      
      for (size_t nodeIndex = 0; nodeIndex < NumberOfNodesPerElement /* TODO ... */ ; ++nodeIndex) {
             
        size_t nodeId = elementNodeIds[nodeIndex]; //... save the node ID
        for (size_t i = 0; i < DegreesOfFreedom; ++i) {
          forceVector(nodeId * DegreesOfFreedom + i) += elementForces[nodeIndex](i); 
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
    
    // Todo:
    // 1) Loop through all elements (Note: there are _elements.size() elements)
    // 2) Create an Element called element corresponding to the elementIndexth element from _elements
    // 3) Get the Node-IDs of that element - Note: the function getNodeIds() is provided by the element
    //    and returns an array<size_t, NumberOfNodesPerElement>
    // 4) Find the element displacements of the element using the Utility-function:
    //    Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds, displacements);
    //    Why? displacements is a vector which has a displacement for every node. We want an array
    //    with NumberOfNodesPerElement ElementVector's that corresponds to the displacement at the
    //    nodes with elementNodeIds.
    //    This also tells you the which array<...,...> you have to use ... :)
    // 5) Evaluate the stiffness tensor corresponding to the displacements of that element
    // 6) So we have the stiffness tensor >for our element<. Great! We now need to include it in the
    //    appropriate location of our >global< stiffness matrix. This corresponds to the line written
    //    in the innermost for-loop. Please finish it :)
    //
    // NOTE: The for-loops outlined simply reflect our recommendations. You may be able to solve it with
    //       less or with more. So long as the result is not hard-coded, please feel free to append or
    //       erase for-loops to your liking :)
    
    for (size_t elementIndex = 0; elementIndex < _elements.size(); ++elementIndex) { // (Step 1 - given)
           
        Element element = _elements[elementIndex];
	    array<size_t, NumberOfNodesPerElement> elementNodeIds = element.getNodeIds();
		array<ElementVector, NumberOfNodesPerElement> elementDisplacements = Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds, displacements);
        ElementStiffnessMatrix elementStiffness = element.computeStiffnessMatrix(elementDisplacements);
      
      for (size_t nodeIndex1 = 0; nodeIndex1 < NumberOfNodesPerElement/* TODO ... */ ; ++nodeIndex1) {
             
        size_t nodeId1 = elementNodeIds[nodeIndex1]; // ...
        
        for (size_t nodeIndex2 = 0; nodeIndex2 < NumberOfNodesPerElement /* TODO ... */; ++nodeIndex2) {
               
          size_t nodeId2 = elementNodeIds[nodeIndex2]; // ...
          
          for (size_t i = 0; i < DegreesOfFreedom; ++i) {
            for (size_t j = 0; j < DegreesOfFreedom; ++j) {
               
              stiffnessMatrix(nodeId1 * DegreesOfFreedom + i,
                              nodeId2 * DegreesOfFreedom + j) += elementStiffness(nodeIndex1 * DegreesOfFreedom + i, nodeIndex2 * DegreesOfFreedom + j); // (Step 6) add the stiffness compon
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
