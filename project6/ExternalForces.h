#ifndef ELEMENT_EXTERNAL_FORCE
#define ELEMENT_EXTERNAL_FORCE

#include "Definitions.h"

namespace Elements {
namespace ExternalForce {

template <class Element, template<class> class NodeType = NodeWithId> 
class ConstantBodyForce {

public:

  static const unsigned int NumberOfNodes    = Element::NumberOfNodes;
  static const unsigned int SpatialDimension = Element::SpatialDimension;
  static const unsigned int DegreesOfFreedom = Element::SpatialDimension;

  typedef typename Element::Point               Point;
  typedef typename Element::Vector              Vector;
  typedef typename Element::NodalDisplacements  NodalDisplacements;
  typedef typename Element::NodeStiffnessMatrix NodeStiffnessMatrix;
  typedef typename Element::Forces              Forces;
  typedef typename Element::Strain              Strain;
  typedef typename Element::Stress              Stress;
  typedef typename Element::TangentMatrix       TangentMatrix;
  typedef typename Element::StiffnessMatrix     StiffnessMatrix;
  typedef typename Element::MassMatrix          MassMatrix;
  
  typedef NodeType<Point>                       Node;

  ConstantBodyForce(const Element & element, const Vector & bodyForceVector) {
    
    _nodeIds = element._nodeIds;
    
    array<double, NumberOfNodes> nodalWeights; 
    
    nodalWeights.fill( element._properties._density * element._volume / (double)NumberOfNodes );
    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++){
      _nodalForces[nodeIndex] = bodyForceVector*nodalWeights[nodeIndex];
    }
    
  }
  
  double
  computeEnergy(const NodalDisplacements & displacements) const {
    
    double energy = 0.;
    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++) {
      energy -= _nodalForces[nodeIndex].dot(displacements[nodeIndex]);
    }
    
    return energy;
  }

  Forces
  computeForces(const NodalDisplacements &  displacements) const {
    
    ignoreUnusedVariables(displacements);
    
    Forces nodalForces;
    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++) {
      nodalForces[nodeIndex] = -_nodalForces[nodeIndex];
    }
    
    return nodalForces;
  }
  
  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements & displacements) const {
    ignoreUnusedVariables(displacements, time);
    return StiffnessMatrix::Zero();
  }

  
  MassMatrix
  computeConsistentMassMatrix() const {
    return MassMatrix::Zero();
  }

  array<size_t, NumberOfNodes>
  getNodeIds() const {
    return _nodeIds;
  }

  
private:

  array<Vector, NumberOfNodes> _nodalForces;
  array<size_t, NumberOfNodes> _nodeIds;
};

}
}

#endif //ELEMENT_EXTERNAL_FORCE
