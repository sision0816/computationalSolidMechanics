#ifndef ELEMENT_WALL_H
#define ELEMENT_WALL_H

#include "Definitions.h"

namespace Elements {
namespace ExternalForce {

template <unsigned SpatialDimensionSize,
          unsigned DegreesOfFreedomSize,
          template<class> class NodeType = NodeWithId>
class Wall {

public:

  static const unsigned int NumberOfNodes                   = 1;
  static const unsigned int SpatialDimension                = SpatialDimensionSize;
  static const unsigned int DegreesOfFreedom                = DegreesOfFreedomSize;

  typedef Matrix<double, SpatialDimension, 1>                 Point;
  typedef NodeType<Point>                                     Node;
  typedef Matrix<double, DegreesOfFreedom, 1>                 Vector;
  typedef array<Vector, 1>                                    NodalDisplacements;
  typedef Matrix<double, DegreesOfFreedom, DegreesOfFreedom>  StiffnessMatrix;
  typedef Matrix<double,
                 NumberOfNodes*DegreesOfFreedom,
                 NumberOfNodes*DegreesOfFreedom>              MassMatrix;
  typedef array<Vector, 1>                                    Forces;
  typedef Matrix<double, 1, 1>                                Stress;
  typedef Matrix<double, 1, 1>                                Strain;

  Wall(const Node &   node          ,
       const double & forceConstant ,
       const Vector & base          ,
       const Point &  direction     ) :
    _nodeId(node._id),
    _startingPosition(node._position),
    _forceConstant(forceConstant),
    _direction(direction / direction.norm()),
    _indenterBase(base) {
  }

  double
  computeEnergy(const NodalDisplacements & displacements) const {
    
    Point shortDisplacements;
    for (size_t spatialIndex = 0; spatialIndex < SpatialDimension;  spatialIndex++){
      shortDisplacements(spatialIndex) = displacements[0](spatialIndex);
    }
    const Point r = (_startingPosition + shortDisplacements) - _indenterBase;
    const double penetration = _direction.dot(r);
    if (penetration > 0) {
      return _forceConstant * penetration * penetration * penetration;
    } else {
      return 0.;
    }
  }

  Forces
  computeForces(const NodalDisplacements & displacements) const {
    
    Point shortDisplacements;
    for (size_t spatialIndex = 0; spatialIndex < SpatialDimension;  spatialIndex++){
      shortDisplacements(spatialIndex) = displacements[0](spatialIndex);
    }
    const Point r = (_startingPosition + shortDisplacements) - _indenterBase;
    const double penetration = _direction.dot(r);
    Vector longDirection = Vector::Zero();
    for (size_t spatialIndex = 0; spatialIndex < SpatialDimension;  spatialIndex++){
      longDirection(spatialIndex) = _direction(spatialIndex);
    }
    if (penetration > 0) {
      const Vector force =
        3. * _forceConstant * penetration * penetration * longDirection;
      return (array<Vector, 1>) {{force}};
    } else {
      return (array<Vector, 1>) {{Vector::Zero()}};
    }
  }

  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements & displacements) const {
    
    Point shortDisplacements;
    for (size_t spatialIndex = 0; spatialIndex < SpatialDimension;  spatialIndex++){
      shortDisplacements(spatialIndex) = displacements[0](spatialIndex);
    }
    const Point r = (_startingPosition + shortDisplacements) - _indenterBase;
    Vector longDirection = Vector::Zero();
    for (size_t spatialIndex = 0; spatialIndex < SpatialDimension;  spatialIndex++){
      longDirection(spatialIndex) = _direction(spatialIndex);
    }
    const double penetration = _direction.dot(r);
    if (penetration > 0) {
      return 6 * _forceConstant * penetration * (longDirection * longDirection.transpose());
    } else {
      return StiffnessMatrix::Zero();
    }
  }

  MassMatrix
  computeConsistentMassMatrix() const {
    return MassMatrix::Zero();
  }

  array<size_t, 1>
  getNodeIds() const {
    return (array<size_t, 1>) {{_nodeId}};
  }

private:

  const size_t _nodeId;
  const Point _startingPosition;
  const double _forceConstant;
  Point _direction;
  Point _indenterBase;
};
}// namespace Elements
}// namespace ExternalForce

#endif // ELEMENT_WALL_H
